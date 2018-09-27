

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


function add_feas_cut!(SCENS, v_xs, numScens, ys, names2)
    # max (b − D ˆy)Tu
    # s.t. ATu ≤ c
    # u ≥ 0
    # agregar (b − D ˆy)Tu ≤ 0
    for k = 1:numScens
        @constraint(v_xs[k][1], sum(SCENS[k].q[i] * ys[k][names2[i]] for i = 1:length(names2)) >= 0)
    end
end

function update_subprob_values(v_xs, numScens, names1, SCENS, is_integer)
    v_x_hat = 0.0
    π_hat = []
    if !is_integer
        for k = 1:numScens
            status = solve(v_xs[k][1])
            if status == :InfeasibleOrUnbounded
                return nothing, nothing
            end
            π_k = [getdual(v_xs[k][2][i]) for i in 1:length(names1)]
            push!(π_hat, π_k)
            v_x_hat += SCENS[k].p * v_xs[k][1].objVal
        end
    else
        for k = 1:numScens
            status = solve(v_xs[k][1])
            if status == :InfeasibleOrUnbounded
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

