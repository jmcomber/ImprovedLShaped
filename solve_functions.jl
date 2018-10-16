function solve_improved(master, x, master_data, SCENS, CORE, SEC_STG_COLS, names1)
    change_category!(x, master_data.cols, true)
    solve(master, suppress_warnings=true)
    x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

    # # CREAR v_x con referencia a constraints de las que necesito dual (z = x), para poder obtener duales (getdual(constr))
    v_xs, ys = create_v_xs(x_hat, SCENS, CORE[4], master_data.cols, SEC_STG_COLS)

    for y in ys
        change_category!(y, SEC_STG_COLS, true)
    end

    V, W = Set([]), Set([])

    function add_lazy_improved(cb)
        x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]
        x_hat = [round(i, 0) for i in x_hat]

        update_subproblems!(v_xs, x_hat)

        for y in ys
            change_category!(y, SEC_STG_COLS, false)
        end

        v_x_hat, π_hat = update_subprob_values(v_xs, names1, SCENS, false)
        if v_x_hat !== nothing   #optimality
            push!(V, x_hat)
        
            if θ_hat < v_x_hat - τ
                β = sum(SCENS[k].p * (v_xs[k][1].objVal - π_hat[k]'x_hat) for k in 1:length(SCENS))
                α = sum(SCENS[k].p * π_hat[k] for k in 1:length(SCENS))
                @lazyconstraint(cb, θ >= sum(α[i] * x[names1[i]] for i in 1:length(names1)) + β)
                return
            end

            if x_hat ∈ V          
                for y ∈ ys
                    change_category!(y, SEC_STG_COLS, true)
                end
                v_x_hat, π_hat = update_subprob_values(v_xs, names1, SCENS, true)
                S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
                if v_x_hat !== nothing
                    @lazyconstraint(cb, θ >= (L - v_x_hat) * (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))) + v_x_hat)
                    if x_hat ∉ W
                        push!(W, x_hat)
                    end
                else
                    S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
                    @lazyconstraint(cb, 1 <= (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))))
                end
            end
        else    #feasibility
            println("Feasibility")
            S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
            # CREO QUE ESTA DEBERÍA SER BENDERS, NO INTEGER
            @lazyconstraint(cb, 1 <= (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))))
        end            
    end
    addlazycallback(master, add_lazy_improved)
    solve(master, suppress_warnings=true)
end