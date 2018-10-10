

struct Stoch_Scenario
    # p probability, d right side, q obj, a coeffs, comps either "E"(quality), "G"(reater or equal) or "L"(ess or equal)
    p::Float64
    h::Vector{Float64}
    q::Vector{Float64}
    T:: Array{Float64,2}
    W:: Array{Float64,2}
    comps::Array{String, 1}
end


function add_cut!(master, SCENS, v_xs, π_hat, x_hat, numScens, x, θ, names1, is_integer)
    
    if !is_integer
        β = sum(SCENS[k].p * (v_xs[k][1].objVal - π_hat[k]'x_hat) for k in 1:numScens)
        α = sum(SCENS[k].p * π_hat[k] for k in 1:numScens)
        @constraint(master, θ >= sum(α[i] * x[names1[i]] for i in 1:length(names1)) + β)
    # else
    #     # Integer L-Shaped
    #     # theta >= (L - v(x'))[SUM_{i en S(x')}(1 - x_i) + SUM_{i no en S(x')}(x_{i})] + v(x')
    #     S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
    #     @constraint(master, θ >= (L - v_x_hat) * (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))) + v_x_hat)
    end
end

function get_duals_constr(v_xs, numScens, k)
    λ1 = getdual(v_xs[k][3])
    λ2 = getdual(v_xs[k][2])
    λ1, λ2
end

function add_feas_cut!(master, x, names1, SCENS, v_xs, numScens, ys, names2, k)
    # v_xs[i][3] tiene los constr_λ de ese problema, y v_xs[i][2] los constr_π.
    # primero conseguir valores duales, después multiplicarlo por las matrices, y armar la restricción
    # conseguir duales
    λ1, λ2 = get_duals_constr(v_xs, numScens, k)
    println(size(λ1), " ", size(λ2))
    println(size(names2), " ", size(names1))
    #multiplicar por matrices (SCENS[k].T, SCENS[k].W, SCENS[k].h)
    #println(length(λ1[1]), " ", length(SCENS[1].h), " ", length(λ2[1]))
    
    # bar_a es λ·A, donde A en este caso es [T W]
    #                                       [I 0]
    
    A1 = hcat(SCENS[k].T, SCENS[k].W)
    # I de tamaño (length(x), length(names1) + length(names2))?
    A2 = eye(length(x), length(names1) + length(names2))
    A = vcat(A1, A2)
    λ = vcat(λ1, λ2)
    bar_a = λ'A
    upper_affine1 = 0
    for i in 1:length(names2)
        if getupperbound(ys[k][names2[i]]) !== Inf && bar_a[i] < 0
            println(!"Algo")
            upper_affine1 += bar_a[i] * getupperbound(ys[k][names2[i]])
        end
    end

    upper_affine2 = 0
    for i in 1:length(names1)
        if getupperbound(x[names1[i]]) !== Inf && bar_a[length(names2) + i] < 0
            upper_affine2 += bar_a[length(names2) + i] * getupperbound(x[names1[i]])
        end
    end
    # println(upper_affine1, " ", upper_affine2)
    upper_affine = upper_affine1 + upper_affine2
    
    lower_affine1 = 0
    for i in 1:length(names2)
        if getlowerbound(ys[k][names2[i]]) !== -Inf && bar_a[i] > 0
            # println(!"Algo")
            lower_affine1 += bar_a[i] * getlowerbound(ys[k][names2[i]])
        end
    end

    lower_affine2 = 0
    for i in 1:length(names1)
        if getlowerbound(x[names1[i]]) !== -Inf && bar_a[length(names2) + i] > 0
            lower_affine2 += bar_a[length(names2) + i] * getlowerbound(x[names1[i]])
        end
    end
    # println(lower_affine1, " ", lower_affine2)
    lower_affine = lower_affine1 + lower_affine2
    
    @constraint(master, λ1'SCENS[k].h + sum(λ2[i] * x[names1[i]] for i in 1:length(names1)) <= upper_affine + lower_affine)

end

function update_subprob_values(v_xs, numScens, names1, SCENS, is_integer)
    v_x_hat = 0.0
    π_hat = []
    if !is_integer
        for k = 1:numScens
            status = solve(v_xs[k][1], suppress_warnings=true)
            if status == :InfeasibleOrUnbounded || status == :Infeasible
                return nothing, k
            end
            π_k = [getdual(v_xs[k][2][i]) for i in 1:length(names1)]
            push!(π_hat, π_k)
            v_x_hat += SCENS[k].p * v_xs[k][1].objVal
        end
    else
        for k = 1:numScens
            status = solve(v_xs[k][1], suppress_warnings=true)
            if status == :InfeasibleOrUnbounded || status == :Infeasible
                return nothing, nothing
            end
            v_x_hat += SCENS[k].p * v_xs[k][1].objVal
        end
    end
    return v_x_hat, π_hat
end

function print_iter_info(iter_count, θ_hat, v_x_hat, objVal)
    println("Iteración ", iter_count)
    println("θ_hat: ", θ_hat)
    println("v_x_hat: ", v_x_hat)
    println("Master objective: ", objVal)
    nothing
end

function change_category!(vars, COLS, to_int)
    if to_int
        for (name, var_type) in COLS
            if var_type == "int"
                setcategory(vars[name], :Int)
            elseif var_type == "bin"
                setcategory(vars[name], :Bin)
            end
        end
    else 
        for (name, var_type) in COLS
            setcategory(vars[name], :Cont)
        end
    end

end
