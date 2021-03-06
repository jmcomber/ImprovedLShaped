
struct Subproblem
    model::JuMP.Model
    constr_π::Array{JuMP.ConstraintRef,1}
    constr_λ::Array{Any,1}
end

mutable struct StageData
    cols::Array{Array{String,1},1}
    rhs::Dict{Any,Any}
    obj::Dict{SubString{String},Float64}
    rows::Array{Any,1}
    constrs::Dict{SubString{String},Dict{Any,Any}}
    L::Float64
end

struct Stoch_Scenario
    # p probability, d right side, q obj, a coeffs, comps either "E"(quality), "G"(reater or equal) or "L"(ess or equal)
    p::Float64
    h::Vector{Float64}
    q::Vector{Float64}
    T:: Array{Float64,2}
    W:: Array{Float64,2}
    comps::Array{String, 1}
end

function not_suited(linking_vars, first_cols)
    for link_var in LINKING_VARS
        for (name, var_type) in first_cols
            if link_var == name && var_type != "bin"
                return true
            end
        end
    end
    false
end


function add_cont_optimality_cut!(master, SCENS, v_xs, π_hat, x_hat, x, θ, names1)
    β = sum(SCENS[k].p * (v_xs[k].model.objVal - π_hat[k]'x_hat) for k in 1:length(SCENS))
    α = sum(SCENS[k].p * π_hat[k] for k in 1:length(SCENS))
    @constraint(master, θ >= sum(α[i] * x[names1[i]] for i in 1:length(names1)) + β)
end

function get_duals_constr(v_x)
    λ1 = getdual(v_x.constr_λ)    # Tz + Wy ~ h
    λ2 = getdual(v_x.constr_π)    # z = x
    λ1, λ2
end

function add_cont_feas_cut!(master, x, names1, scen, v_x)
    λ1, λ2 = get_duals_constr(v_x)
    @constraint(master, λ1'scen.h + sum(λ2[i] * x[names1[i]] for i in 1:length(names1)) <= 0)
end

function update_subprob_values(v_xs, names1, SCENS, is_integer)
    v_x_hat = 0.0
    π_hat = []
    if !is_integer
        for k = 1:length(SCENS)
            status = solve(v_xs[k].model, suppress_warnings=true)
            if status == :InfeasibleOrUnbounded || status == :Infeasible
                return nothing, k
            end
            π_k = [getdual(v_xs[k].constr_π[i]) for i in 1:length(names1)]
            push!(π_hat, π_k)
            v_x_hat += SCENS[k].p * v_xs[k].model.objVal
        end
    else
        for k = 1:length(SCENS)
            status = solve(v_xs[k].model, suppress_warnings=true)
            if status == :InfeasibleOrUnbounded || status == :Infeasible
                return nothing, nothing
            end
            v_x_hat += SCENS[k].p * v_xs[k].model.objVal
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

