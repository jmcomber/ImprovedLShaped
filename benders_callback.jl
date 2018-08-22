using JuMP, Gurobi
include("smps_parser.jl")


struct Stoch_Scenario
    # p probability, d right side, q obj, a coeffs, comps either "E"(quality), "G"(reater or equal) or "L"(ess or equal)
    p::Float64
    h::Vector{Float64}
    q::Vector{Float64}
    T:: Array{Float64,2}
    W:: Array{Float64,2}
    comps::Array{String, 1}
end


function update_subprob_values(v_xs, numScens, names1, SCENS, is_integer)
    v_x_hat = 0.0
    π_hat = []
    if !is_integer
        for k = 1:numScens
            solve(v_xs[k][1])
            π_k = [getdual(v_xs[k][2][i]) for i in 1:length(names1)]
            push!(π_hat, π_k)
            v_x_hat += SCENS[k].p * v_xs[k][1].objVal
        end
    else
        for k = 1:numScens
            solve(v_xs[k][1])
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


# time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape\\shape-3-3_3-3-2_1.tim"
# core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape\\shape-3-3_3-3-2_1.mps"
# stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape\\shape-3-3_3-3-2_1.sto"

time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_10_50_100.tim"
core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_10_50_100.cor"
stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_10_50_100.sto"

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

change_category!(x, FIRST_STG_COLS, true)


solve(master)

x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

# # CREAR v_x con referencia a constraints de las que necesito dual (z = x), para poder obtener duales (getdual(constr))
v_xs, ys = create_v_xs(x_hat, SCENS, numScens, CORE[4], FIRST_STG_COLS, SEC_STG_COLS)

for y in ys
    change_category!(y, SEC_STG_COLS, true)
end


τ = 1e-4

V = Set([])

function add_lazy_ilsm(cb)
    x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]
    x_hat = [round(i, 0) for i in x_hat]
    # println("x_hat: $x_hat, θ_hat: $θ_hat")


    update_subproblems!(v_xs, x_hat)

    for y in ys
        change_category!(y, SEC_STG_COLS, false)
    end

    names1 = [var[1] for var in FIRST_STG_COLS]

    v_x_hat, π_hat = update_subprob_values(v_xs, numScens, names1, SCENS, false)
    push!(V, x_hat)
    
    if θ_hat < v_x_hat - τ
        # println("Going to add subgradient cut")
        β = sum(SCENS[k].p * (v_xs[k][1].objVal - π_hat[k]'x_hat) for k in 1:numScens)
        α = sum(SCENS[k].p * π_hat[k] for k in 1:numScens)
        @lazyconstraint(cb, θ >= sum(α[i] * x[names1[i]] for i in 1:length(names1)) + β)
        return
    end

    if x_hat ∈ V
        
        for y ∈ ys
            change_category!(y, SEC_STG_COLS, true)
        end
    
        v_x_hat, π_hat = update_subprob_values(v_xs, numScens, names1, SCENS, true)
        S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
            @lazyconstraint(cb, θ >= (L - v_x_hat) * (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))) + v_x_hat)
        if x_hat ∉ W
            push!(W, x_hat)
            # println("Going to add integer cut")
            # S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
            # @lazyconstraint(cb, θ >= (L - v_x_hat) * (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))) + v_x_hat)
        end
    end


    # v_xs = [[modelo1, constr_pi1], [...], ...]
    # pi_hat[k] es valor del pi para problema del escenario k en la actual iteración

    # v_x_hat = 0.0
    # π_hat = []

    # names1 = [var[1] for var in FIRST_STG_COLS]

    # for k = 1:numScens
    #     solve(v_xs[k][1])
    #     π_k = [getdual(v_xs[k][2][i]) for i in 1:length(names1)]
    #     push!(π_hat, π_k)
    #     v_x_hat += SCENS[k].p * v_xs[k][1].objVal
    # end

    # println("v_x_hat: $v_x_hat")
    # if θ_hat < v_x_hat - τ
    #     S = [i for i in 1:length(names1) if x_hat[i] >= 0.9]
    #     @lazyconstraint(cb, θ >= (L - v_x_hat) * (sum(1 - x[names1[i]] for i in S) + sum(x[names1[i]] for i in 1:length(names1) if !(i in S))) + v_x_hat)
    # end

end


addlazycallback(master, add_lazy_ilsm)


solve(master)

println("\nOptimal value Master MIP: ", master.objVal)


for i in 1:master.numCols
    if master.colVal[i] != 0.0
        println(master.colNames[i], ": ", master.colVal[i])
    end
end






