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
    name::String
end


time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_100.tim"
core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_100.cor"
stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_100.sto"

TIME = get_TIME(time)

STG2_C, STG2_R = TIME[2][1], TIME[2][2]

CORE = get_CORE(core)

STOCH = get_STOCH(stoch)
STOCH, perturb = STOCH

numScens = length(STOCH)


FIRST_STG_COLS = get_first_stg_cols(CORE, STG2_C)


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

println(3)

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

println(5)

SCENS = create_scenarios(numScens, STOCH, SEC_STG_RHS, SEC_STG_CONSTR, SEC_STG_ROWS, SEC_STG_COLS, FIRST_STG_COLS)

non_neg = Dict(i[1] => true for i=FIRST_STG_COLS)

# CONTS_IDX = []
# INTS_IDX = []
# BINS_IDX = []


# for (name, type_var) in FIRST_STG_COLS
#     if type_var == "cont"
#         push!(CONTS_IDX, name)
#     elseif type_var == "int"
#         # revisar si tiene UP o UI 1. Si no, INTS_IDX. Si es que sí, BINS_IDX.
#         added = false
#         for bound in CORE[4]
#             if bound[3] == name && !added
#                 if (bound[1] == "UP" || bound[1] == "UI") && (get(bound[4]) == 1 || get(bound[4]) == 1.0)
#                     push!(BINS_IDX, name)
#                     added = true
#                 end
#             end
#         end
#         if !added
#             push!(INTS_IDX, name)
#         end
#     else
#         push!(BINS_IDX, name)
#     end
# end

# println("BINS_IDX: ", BINS_IDX)
println("[INITIATING PREPROCESSING]")

# ver que linking_vars \subseteq de BINS_IDX

LINKING_VARS = linking_vars(SEC_STG_ROWS, FIRST_STG_COLS, SEC_STG_CONSTR)
println("LINKING VARS: ", LINKING_VARS)

# Ahora crear los L_k y L
# A son los coeficientes de FIRST_STG_ROWS
# b es FIRST_STG_RHS
# bounds para z_k son las de primera etapa
# bounds para y_k son las de segunda etapa

L = 0.0
for i in 1:numScens
    Lk = create_Lk(FIRST_STG_COLS, FIRST_STG_ROWS, FIRST_STG_RHS, SEC_STG_COLS, SCENS[i], CORE[4])
    solve(Lk)
    L += SCENS[i].p * Lk.objVal
end
println("[PREPROCESSING READY, L = ", L, "]")
# L = -264.16

# Crear Maestro
# min c^{T}x + theta
#       Ax ~ b
#       x ~ 0 (bounds)
#       theta >= L

master, x, θ = init_master(FIRST_STG_COLS, FIRST_STG_OBJECT, FIRST_STG_ROWS, FIRST_STG_RHS, GurobiSolver, CORE[4])

solve(master)

x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

# CREAR v_x con referencia a constraints de las que necesito dual (z = x), para poder obtener duales (getdual(constr))
# JuMP.setRHS(mycon, 3)  # Now x + y <= 3
v_xs = create_v_xs(x_hat, SCENS, numScens, CORE[4], FIRST_STG_COLS, SEC_STG_COLS)

# v_xs = [[modelo1, constr_pi1], [...], ...]
# pi_hat[k] es valor del pi

v_x_hat = 0.0
π_hat = []

names1 = [var[1] for var in FIRST_STG_COLS]

for k = 1:numScens
    solve(v_xs[k][1])
    π_k = [getdual(v_xs[k][2][i]) for i in 1:length(names1)]
    push!(π_hat, π_k)
    v_x_hat += SCENS[k].p * v_xs[k][1].objVal
end

τ = 1e-4

iter_count = 0
while θ_hat < v_x_hat - τ
    println("Iteración ", iter_count)
    iter_count += 1
    println("θ_hat: ", θ_hat)
    println("v_x_hat: ", v_x_hat)
    # Create cut
    β = sum(SCENS[k].p * (v_xs[k][1].objVal - π_hat[k]'x_hat) for k in 1:numScens)
    α = sum(SCENS[k].p * π_hat[k] for k in 1:numScens)
    
    @constraint(master, θ >= sum(α[i] * x[names1[i]] for i in 1:length(names1)) + β)
    solve(master)
    x_hat, θ_hat = master.colVal[1:end-1], master.colVal[end]

    # Change v_xs, solve them and update π_hat
    update_subproblems!(v_xs, x_hat)
    v_x_hat = 0.0
    π_hat = []

    for k = 1:numScens
        solve(v_xs[k][1])
        π_k = [getdual(v_xs[k][2][i]) for i in 1:length(names1)]
        push!(π_hat, π_k)
        v_x_hat += SCENS[k].p * v_xs[k][1].objVal
    end

end

println("\nOptimal value: ", master.objVal)







