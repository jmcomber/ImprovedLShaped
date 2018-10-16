using JuMP, Gurobi
include("smps_parser.jl")
include("utils.jl")
include("solve_functions.jl")


time = "/Users/jmcomber/Universidad/MISTI/SMPS_Parser/sslp/sslp_15_45_5.tim"
core = "/Users/jmcomber/Universidad/MISTI/SMPS_Parser/sslp/sslp_15_45_5.cor"
stoch = "/Users/jmcomber/Universidad/MISTI/SMPS_Parser/sslp/sslp_15_45_5.sto"

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
FIRST_STG_RHS, FIRST_STG_CONSTR = get_first_stg_rhs(CORE, FIRST_STG_ROWS, FIRST_STG_CONSTR)
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
SEC_STG_RHS, SEC_STG_CONSTR = get_sec_stg_rhs(CORE, SEC_STG_ROWS, SEC_STG_CONSTR)
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
SCENS = create_scenarios(numScens, STOCH, SEC_STG_RHS, SEC_STG_CONSTR, SEC_STG_ROWS, SEC_STG_COLS, FIRST_STG_COLS, AUX2)
non_neg = Dict(i[1] => true for i=FIRST_STG_COLS)
# Fin Utilidades

println("[INITIATING PREPROCESSING]")

# ver que linking_vars \subseteq de BINS_IDX
LINKING_VARS = linking_vars(SEC_STG_ROWS, FIRST_STG_COLS, SEC_STG_CONSTR)

if not_suited(LINKING_VARS, FIRST_STG_COLS)
    throw("Problem not suited for Integer L-Shaped Method: linking variables are not all binary")
end

# Ahora crear los L_k y L
# A son los coeficientes de FIRST_STG_ROWS
# b es FIRST_STG_RHS
master_data = StageData(FIRST_STG_COLS, FIRST_STG_RHS, FIRST_STG_OBJECT, FIRST_STG_ROWS, FIRST_STG_CONSTR, -Inf)

L = 0.0
for i in 1:numScens
    Lk = create_Lk(master_data, SEC_STG_COLS, SCENS[i], CORE[4])
    solve(Lk)
    L += SCENS[i].p * Lk.objVal
end

master_data.L = L
println("[PREPROCESSING READY, L = ", master_data.L, "]")


# Crear Maestro
# min c^{T}x + θ
#       Ax ~ b
#       x ~ 0 (bounds)
#       θ >= L
master, x, θ = init_master(master_data, GurobiSolver, CORE[4])

names1 = [var[1] for var in master_data.cols]
names2 = [var[1] for var in SEC_STG_COLS]

τ = 1e-4



function solve_decomposed(master, improved)
    if improved
        solve_improved(master, x, master_data, SCENS, CORE, SEC_STG_COLS, names1)
    else
        solve(master, suppress_warnings=true)
        x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

        # # CREAR v_x con referencia a constraints de las que necesito dual (z = x), para poder obtener duales (getdual(constr))
        v_xs, ys = create_v_xs(x_hat, SCENS, CORE[4], master_data.cols, SEC_STG_COLS)

        # # v_xs = [[modelo1, constr_pi1], [...], ...]
        # # pi_hat[k] es valor del pi para problema del escenario k en la actual iteración

        v_x_hat, π_hat = update_subprob_values(v_xs, names1, SCENS, true)

        # iter_count = 0
        while v_x_hat === nothing || θ_hat < v_x_hat - τ
              
            solve(master, suppress_warnings=true)
            x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]
            update_subproblems!(v_xs, x_hat)

            v_x_hat, π_hat = update_subprob_values(v_xs, names1, SCENS, false)
            
            if v_x_hat !== nothing
                add_cut!(master, SCENS, v_xs, π_hat, x_hat, x, θ, names1, false)
            else
                println("FEASIBILITY! \n")
                k = π_hat
                add_feas_cut!(master, x, names1, SCENS, v_xs, ys, names2, k)
            end
            
        end

        println("\nOptimal value Master LP: ", master.objVal)

        # Set category for vars
        change_category!(x, master_data.cols, true)
        for y in ys
            change_category!(y, SEC_STG_COLS, true)
        end

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

    end
end


solve_decomposed(master, false)

println("\nOptimal value Master MIP: ", master.objVal)

for i in 1:master.numCols
    if master.colVal[i] != 0.0
        println(master.colNames[i], ": ", master.colVal[i])
    end
end






