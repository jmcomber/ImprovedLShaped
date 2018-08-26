using JuMP, Gurobi
include("smps_parser.jl")
include("utils.jl")

function change_category_2!(vars, COLS, to_int)
    if to_int
        for (name, var_type) in COLS
            if var_type == "int"
                for k = 1:numScens
                    setcategory(vars[k, name], :Int)
                end
            elseif var_type == "bin"
                for k = 1:numScens
                    setcategory(vars[k, name], :Bin)
                end
            end
        end
    else 
        for (name, var_type) in COLS
            for k = 1:numScens
                setcategory(vars[k, name], :Cont)
            end
        end
    end

end



function create_v_x(x_hat, SCENS, numScens, BOUNDS, FIRST_STG_COLS, SEC_STG_COLS)
    v_x = Model(solver=GurobiSolver(OutputFlag=0))

    names1 = [var[1] for var in FIRST_STG_COLS]
    @variable(v_x, z[i = names1])
    names2 = [var[1] for var in SEC_STG_COLS]
    @variable(v_x, y[k = 1:numScens, j = names2])


    @objective(v_x, :Min, sum(SCENS[k].q[j] * y[k, names2[j]] for k=1:numScens for j in 1:length(names2)))

    # Tz + Wy ~ h
    # Its duals are λ: not necessary to get reference
    for k = 1:numScens
        for r = 1:size(SCENS[k].T, 1)
            if SCENS[k].comps[r] == "E"
                @constraint(v_x, sum(SCENS[k].T[r, c] * z[names1[c]] for c in 1:size(SCENS[k].T, 2)) + 
                sum(SCENS[k].W[r, c] * y[k, names2[c]] for c in 1:size(SCENS[k].W, 2)) == SCENS[k].h[r])
            elseif SCENS[k].comps[r] == "G"
                @constraint(v_x, sum(SCENS[k].T[r, c] * z[names1[c]] for c in 1:size(SCENS[k].T, 2)) + 
                sum(SCENS[k].W[r, c] * y[k, names2[c]] for c in 1:size(SCENS[k].W, 2)) >= SCENS[k].h[r])
            else 
                @constraint(v_x, sum(SCENS[k].T[r, c] * z[names1[c]] for c in 1:size(SCENS[k].T, 2)) + 
                sum(SCENS[k].W[r, c] * y[k, names2[c]] for c in 1:size(SCENS[k].W, 2)) <= SCENS[k].h[r])
            end
        end
    end

    # z = x_hat
    # Its duals are π: get reference to return v_x, constr_π
    @constraint(v_x, constr_π[i=1:length(names1)], z[names1[i]] == x_hat[i])
    add_bounds_2!(v_x, numScens, names2, y, BOUNDS, SEC_STG_COLS)

    return [v_x, constr_π], y

end

function add_bounds!(m, names, z, BOUNDS, COLS)
    non_neg = Dict(i[1] => true for i=COLS)
    for bound in BOUNDS
        if bound[3] in names
            if bound[1] == "UP" || bound[1] == "UI"
                @constraint(m, z[bound[3]] <= get(bound[4]))
            elseif bound[1] == "FX"
                @constraint(m, z[bound[3]] == get(bound[4]))
            elseif bound[1] == "LO"
                @constraint(m, z[bound[3]] >= get(bound[4]))
                non_neg[[0, bound[3]]] = false
            elseif bound[1] == "FR"
                non_neg[[0, bound[3]]] = false
            end
        end
    end
    for (key, val) in non_neg
        if val
            # println(key)
            @constraint(m, z[key] >= 0)
        end
    end
end


function add_bounds_2!(m, numScens, names, z, BOUNDS, COLS)
    non_neg2 = Dict(i[1] => true for i=COLS)
    for bound in BOUNDS
        if bound[3] in names
            for k = 1:numScens
                if bound[1] == "UP" || bound[1] == "UI"
                    @constraint(m, z[k, bound[3]] <= get(bound[4]))
                elseif bound[1] == "FX"
                    @constraint(m, z[k, bound[3]] == get(bound[4]))
                elseif bound[1] == "LO"
                    @constraint(m, z[k, bound[3]] >= get(bound[4]))
                    non_neg2[[k, bound[3]]] = false
                elseif bound[1] == "FR"
                    non_neg2[[k, bound[3]]] = false
                end
            end
        end
    end
    for (key, val) in non_neg2
        if val
            println(key)
            try
                @constraint(m, z[key[1], key[2]] >= 0)
            catch
                1+1;
            end
        end
    end
end


function update_subproblem!(x_hat, constr_π)
    for i in 1:length(x_hat)
        JuMP.setRHS(constr_π[i], x_hat[i])
    end
end

function update_subprob_values(v_x, numScens, names1, SCENS, is_integer)
    v_x_hat = 0.0
    π_hat = []
    if !is_integer
        solve(v_x[1])
        π_k = [getdual(v_x[2][i]) for i in 1:length(names1)]
        push!(π_hat, π_k)
        v_x_hat += v_x[1].objVal
    else
        for k = 1:numScens
            solve(v_x[1])
            v_x_hat +=  v_x[1].objVal
        end
    end
    return v_x_hat, π_hat
end

function add_cut!(master, SCENS, v_x, π_hat, x_hat, numScens, x, θ, names1, is_integer)
    
    if !is_integer
        β = v_x[1].objVal - π_hat[1]'x_hat
        α = π_hat
        @constraint(master, θ >= sum(α[i] * x[names1[i]] for i in 1:length(names1)) + β)
    # else
    #     # Integer L-Shaped
    #     # theta >= (L - v(x'))[SUM_{i en S(x')}(1 - x_i) + SUM_{i no en S(x')}(x_{i})] + v(x')
    #     S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
    #     @constraint(master, θ >= (L - v_x_hat) * (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))) + v_x_hat)
    end
end


# time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape\\shape-3-3_3-3-2_1.tim"
# core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape\\shape-3-3_3-3-2_1.mps"
# stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape\\shape-3-3_3-3-2_1.sto"

time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_100.tim"
core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_100.cor"
stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_100.sto"


# Utilidades
TIME = get_TIME(time)
STG2_C, STG2_R = TIME[2][1], TIME[2][2]
CORE = get_CORE(core)
STOCH = get_STOCH(stoch)
STOCH, perturb = STOCH
numScens = length(STOCH)
FIRST_STG_COLS = get_first_stg_cols(CORE, STG2_C)
FIRST_STG_COLS = binarize_first_stg_cols(FIRST_STG_COLS, CORE[4])
just_names = [i[1] for i in FIRST_STG_COLS]
FIRST_STG_OBJECT, object_name = get_first_stg_object(CORE)
AUX = [i[1] for i in FIRST_STG_OBJECT]
for var in FIRST_STG_COLS
    if !(var[1] in AUX)
        push!(FIRST_STG_OBJECT, [var[1], 0.0])
    end
end
AUX = Dict(i[1] => i[2] for i in FIRST_STG_OBJECT)
FIRST_STG_OBJECT = AUX
FIRST_STG_ROWS = get_first_stg_rows(CORE, STG2_R, object_name)
FIRST_STG_CONSTR = Dict(i[2] => Dict() for i in FIRST_STG_ROWS)
FIRST_STG_RHS = get_first_stg_rhs(CORE, FIRST_STG_ROWS)
SEC_STG_COLS = get_sec_stg_cols(CORE, STG2_C)
SEC_STG_OBJECT = get_sec_stg_object(CORE, object_name)
AUX2 = [i[1] for i in SEC_STG_OBJECT]
for var in SEC_STG_COLS
    if !(var[1] in AUX2)
        push!(SEC_STG_OBJECT, [var[1], 0.0])
    end
end
AUX2 = Dict((i[1], STOCH[s][1][2]) => i[2] for i in SEC_STG_OBJECT for s=1:numScens)
println(2)
SEC_STG_ROWS = get_sec_stg_rows(CORE, STG2_R)
SEC_STG_CONSTR = Dict(i[2] => Dict() for i in SEC_STG_ROWS)
SEC_STG_RHS = get_sec_stg_rhs(CORE, SEC_STG_ROWS)
for row in SEC_STG_ROWS
    for col in FIRST_STG_COLS
        if !haskey(SEC_STG_CONSTR[row[2]], col[1])
            SEC_STG_CONSTR[row[2]][col[1]] = 0.0
        end
    end
    for col in SEC_STG_COLS
        if !haskey(SEC_STG_CONSTR[row[2]], col[1])
            SEC_STG_CONSTR[row[2]][col[1]] = 0.0
        end
    end
end
println(3)
SCENS = create_scenarios(numScens, STOCH, SEC_STG_RHS, SEC_STG_CONSTR, SEC_STG_ROWS, SEC_STG_COLS, FIRST_STG_COLS)
non_neg = Dict(i[1] => true for i=FIRST_STG_COLS)
# Fin Utilidades

println("[INITIATING PREPROCESSING]")

# ver que linking_vars \subseteq de BINS_IDX
LINKING_VARS = linking_vars(SEC_STG_ROWS, FIRST_STG_COLS, SEC_STG_CONSTR)
println("LINKING VARS: ", LINKING_VARS)

for link_var in LINKING_VARS
    for (name, var_type) in FIRST_STG_COLS
        if link_var == name && var_type != "bin"
            throw("Problem not suited for Integer L-Shaped Method: linking variables are not all binary")
        end
    end
end

# Ahora crear los L_k y L
# A son los coeficientes de FIRST_STG_ROWS
# b es FIRST_STG_RHS

L = 0.0
for i in 1:numScens
    Lk = create_Lk(FIRST_STG_COLS, FIRST_STG_ROWS, FIRST_STG_RHS, SEC_STG_COLS, SCENS[i], CORE[4])
    solve(Lk)
    L += SCENS[i].p * Lk.objVal
end
println("[PREPROCESSING READY, L = ", L, "]")

# Crear Maestro
# min c^{T}x + θ
#       Ax ~ b
#       x ~ 0 (bounds)
#       θ >= L
master, x, θ = init_master(FIRST_STG_COLS, FIRST_STG_OBJECT, FIRST_STG_ROWS, FIRST_STG_RHS, GurobiSolver, CORE[4], L)

names1 = [var[1] for var in FIRST_STG_COLS]

function solve_decomposed(master, improved)
    if improved        
        change_category!(x, FIRST_STG_COLS, true)
        solve(master)
        x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

        # # CREAR v_x con referencia a constraints de las que necesito dual (z = x), para poder obtener duales (getdual(constr))
        v_x, y = create_v_x(x_hat, SCENS, numScens, CORE[4], FIRST_STG_COLS, SEC_STG_COLS)

        
        change_category_2!(y, SEC_STG_COLS, true)
        

        τ = 1e-4

        V = Set([])
        W = Set([])

        function add_lazy_improved(cb)
            x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]
            x_hat = [round(i, 0) for i in x_hat]

            update_subproblem!(x_hat, v_x[2])

            # for y in ys
            change_category_2!(y, SEC_STG_COLS, false)
            # end

            names1 = [var[1] for var in FIRST_STG_COLS]

            v_x_hat, π_hat = update_subprob_values(v_x, numScens, names1, SCENS, false)
            push!(V, x_hat)
            
            if θ_hat < v_x_hat - τ
                β = v_x[1].objVal - π_hat[1]'x_hat
                α = π_hat
                @lazyconstraint(cb, θ >= sum(α[i] * x[names1[i]] for i in 1:length(names1)) + β)
                return
            end

            if x_hat ∈ V          
                # for y ∈ ys
                change_category_2!(y, SEC_STG_COLS, true)
                # end
                v_x_hat, π_hat = update_subprob_values(v_x, numScens, names1, SCENS, true)
                S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
                    @lazyconstraint(cb, θ >= (L - v_x_hat) * (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))) + v_x_hat)
                if x_hat ∉ W
                    push!(W, x_hat)
                end
            end
        end
        addlazycallback(master, add_lazy_improved)
        solve(master)
    else
        solve(master)
        x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

        # # CREAR v_x con referencia a constraints de las que necesito dual (z = x), para poder obtener duales (getdual(constr))
        v_x, y = create_v_x(x_hat, SCENS, numScens, CORE[4], FIRST_STG_COLS, SEC_STG_COLS)

        # # v_xs = [[modelo1, constr_pi1], [...], ...]
        # # pi_hat[k] es valor del pi para problema del escenario k en la actual iteración

        v_x_hat = 0.0
        π_hat = []

        # names1 = [var[1] for var in FIRST_STG_COLS]

        
        solve(v_x[1])
        π_k = [getdual(v_x[2][i]) for i in 1:length(names1)]
        push!(π_hat, π_k)
        v_x_hat +=  v_x[1].objVal
    

        τ = 1e-4


        iter_count = 0
        while θ_hat < v_x_hat - τ
            # print_iter_info(iter_count, θ_hat, v_x_hat, master.objVal)
            iter_count += 1
              
            solve(master)
            x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]
            update_subproblem!(x_hat, v_x[2])

            v_x_hat = 0.0
            π_hat = []

            names1 = [var[1] for var in FIRST_STG_COLS]

            solve(v_x[1])
            π_k = [getdual(v_x[2][i]) for i in 1:length(names1)]
            push!(π_hat, π_k)
            v_x_hat += v_x[1].objVal

            add_cut!(master, SCENS, v_x, π_hat, x_hat, numScens, x, θ, names1, false)
        end

        println("\nOptimal value Master LP: ", master.objVal)

        # Set category for vars
        change_category!(x, FIRST_STG_COLS, true)
        
        change_category_2!(y, SEC_STG_COLS, true)
    

        v_x_hat = 1e7

        println("Now MIP")
        function add_lazy_ilsm(cb)
            x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]
            x_hat = [round(i, 0) for i in x_hat]
            names1 = [var[1] for var in FIRST_STG_COLS]
            update_subproblem!(x_hat, v_x[2])
            v_x_hat, π_hat = update_subprob_values(v_x, numScens, names1, SCENS, true)
            S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
            @lazyconstraint(cb, θ >= (L - v_x_hat) * (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))) + v_x_hat)

        end

        addlazycallback(master, add_lazy_ilsm)

        solve(master)

    end

end


solve_decomposed(master, false)

println("\nOptimal value Master MIP: ", master.objVal)

for i in 1:master.numCols
    if master.colVal[i] != 0.0
        println(master.colNames[i], ": ", master.colVal[i])
    end
end






