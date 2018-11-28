include("structs.jl")

function solve_improved(master::JuMP.Model, x::JuMP.JuMPArray{JuMP.Variable,1,Tuple{Array{String,1}}}, θ::JuMP.Variable, master_data::StageData, SCENS::Array{Stoch_Scenario,1}, CORE, SEC_STG_COLS::Array{Array{String,1},1}, names1::Array{String,1}, L::Float64, multi_thread::Bool)
    if multi_thread
        solve(master, suppress_warnings=true)
    else
        solve(master, suppress_warnings=true)
    end
    x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

    v_xs = create_v_xs(x_hat, SCENS, CORE[4], master_data.cols, SEC_STG_COLS, multi_thread)

    integer_v_xs = create_v_xs_integer(x_hat, SCENS, CORE[4], master_data.cols, SEC_STG_COLS, multi_thread)

    for v_x in integer_v_xs
        change_category!(v_x.y, SEC_STG_COLS, true)
    end

    # # pi_hat[k] es valor del pi para problema del escenario k en la actual iteración
    v_x_hat, π_hat = update_subprob_values(integer_v_xs, names1, SCENS, true)

    while v_x_hat === nothing || θ_hat < v_x_hat - τ
        solve(master, suppress_warnings=true)
        x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]
        update_subproblems!(v_xs, x_hat)
        v_x_hat, π_hat = update_subprob_values(v_xs, names1, SCENS, false)
        if v_x_hat !== nothing
            add_cont_optimality_cut!(master, SCENS, v_xs, π_hat, x_hat, x, θ, names1)
        else
            println("FEASIBILITY! \n")
            add_cont_feas_cut!(master, x, names1, SCENS[π_hat], v_xs[π_hat])
        end
    end

    println("\nOptimal value Master LP: ", master.objVal, "\n\n")

    
    change_category!(x, master_data.cols, true)
    
    init = time_ns()

    solve(master, suppress_warnings=true)
    x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

    V, W = Set([]), Set([])
    feas_cont, feas_int, opt_cont, opt_int = 0, 0, 0, 0
    function add_lazy_improved(cb)
        x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

        update_subproblems!(v_xs, x_hat)
        update_subproblems!(integer_v_xs, x_hat)

        v_x_hat, π_hat = update_subprob_values(v_xs, names1, SCENS, false)
        if v_x_hat !== nothing   #optimality
            push!(V, x_hat)
        
            if θ_hat < v_x_hat - τ
                β = sum(SCENS[k].p * (v_xs[k].model.objVal - π_hat[k]'x_hat) for k in 1:length(SCENS))
                α = sum(SCENS[k].p * π_hat[k] for k in 1:length(SCENS))
                opt_cont += 1
                @lazyconstraint(cb, θ >= sum(α[i] * x[names1[i]] for i in 1:length(names1)) + β)
                return
            end

            if x_hat ∈ V          
                v_x_hat, π_hat = update_subprob_values(integer_v_xs, names1, SCENS, true)
                S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
                if v_x_hat !== nothing
                    opt_int += 1
                    @lazyconstraint(cb, θ >= (L - v_x_hat) * (sum(1 - x[names1[i]] for i ∈ S) + sum(x[names1[i]] for i in 1:length(names1) if i ∉ S)) + v_x_hat)
                    push!(W, x_hat)
                else
                    S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
                    feas_int += 1
                    @lazyconstraint(cb, 1 <= (sum(1 - x[names1[i]] for i ∈ S) + sum(x[names1[i]] for i in 1:length(names1) if i ∉ S)))
                end
            end
        else    #feasibility
            println("Feasibility")
            S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
            # CREO QUE ESTA DEBERÍA SER BENDERS, NO INTEGER
            λ1, λ2 = get_duals_constr(v_xs[π_hat])
            feas_cont += 1
            @lazyconstraint(cb, λ1'SCENS[π_hat].h + sum(λ2[i] * x[names1[i]] for i in 1:length(names1)) <= 0)
            # @lazyconstraint(cb, 1 <= (sum(1 - x[names1[i]] for i ∈ S) + sum(x[names1[i]] for i in 1:length(names1) if i ∉ S)))
        end            
    end
    addlazycallback(master, add_lazy_improved)
    solve(master, suppress_warnings=true)
    println("Continuous Feasibility: $feas_cont\nInteger Feasibility: $feas_int\nContinuous Optimality: $opt_cont\nInteger Optimality: $opt_int")
    return (time_ns() - init) / (10^9)
end

function solve_not_improved(master::JuMP.Model, x::JuMP.JuMPArray{JuMP.Variable,1,Tuple{Array{String,1}}}, θ::JuMP.Variable, master_data::StageData, SCENS::Array{Stoch_Scenario,1}, CORE, SEC_STG_COLS, names1::Array{String,1}, names2::Array{String,1}, L::Float64,multi_thread::Bool)
    solve(master, suppress_warnings=true)
    x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

    # # CREAR v_x con referencia a constraints de las que necesito dual (z = x), para poder obtener duales (getdual(constr))
    v_xs = create_v_xs(x_hat, SCENS, CORE[4], master_data.cols, SEC_STG_COLS, multi_thread)

    # # v_xs = [[modelo1, constr_pi1], [...], ...]
    # # pi_hat[k] es valor del pi para problema del escenario k en la actual iteración

    v_x_hat, π_hat = update_subprob_values(v_xs, names1, SCENS, true)

    while v_x_hat === nothing || θ_hat < v_x_hat - τ
          
        solve(master, suppress_warnings=true)
        x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]
        update_subproblems!(v_xs, x_hat)

        v_x_hat, π_hat = update_subprob_values(v_xs, names1, SCENS, false)
        
        if v_x_hat !== nothing
            add_cont_optimality_cut!(master, SCENS, v_xs, π_hat, x_hat, x, θ, names1)
        else
            println("FEASIBILITY! \n")
            add_cont_feas_cut!(master, x, names1, SCENS[π_hat], v_xs[π_hat])
        end
        
    end

    println("\nOptimal value Master LP: ", master.objVal)
    init = time_ns()
    # Set category for vars
    # change_category!(x, master_data.cols, true)
    # for v_x in v_xs
    #     change_category!(v_x.y, SEC_STG_COLS, true)
    # end

    v_x_hat = 1e7

    println("Now MIP")
    function add_lazy_ilsm(cb)
        x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]
        x_hat = [round(i, 0) for i in x_hat]
        update_subproblems!(v_xs, x_hat)
        v_x_hat, π_hat = update_subprob_values(v_xs, names1, SCENS, true)
        if v_x_hat != nothing
            S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
            @lazyconstraint(cb, θ >= (L - v_x_hat) * (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))) + v_x_hat)
        else
            S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
            @lazyconstraint(cb, 1 <= (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))))     
        end
    end

    addlazycallback(master, add_lazy_ilsm)
    solve(master, suppress_warnings=true)
    return (time_ns() - init) / (10^9)
end