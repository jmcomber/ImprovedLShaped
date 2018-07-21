using StochasticPrograms, Gurobi, LShapedSolvers

# Backup en LShapedSolversCopy


struct Stoch_Scenario <: AbstractScenarioData
    # π probability, d right side, q obj, a coeffs, comps either "E"(quality), "G"(reater or equal) or "L"(ess or equal)
    π::Float64
    d::Vector{Float64}
    q::Vector{Float64}
    a:: Array{Float64,2}
    comps::Array{String, 1}
    name::String
end

function get_variable(name, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins)
    if name in CONTS_IDX
        return conts[name]
    elseif name in INTS_IDX
        return ints[name]
    else
        return bins[name]
    end
end

function StochasticPrograms.expected(sds::Vector{Stoch_Scenario})
    sd = Stoch_Scenario(1, sum([s.π * s.d for s in sds]), sum([s.π * s.q for s in sds]))
end


StochasticPrograms.probability(s::Stoch_Scenario) = s.π

include("smps_parser.jl")

# time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape\\shape-3-3_3-3-2_1.tim"
# core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape\\shape-3-3_3-3-2_1.mps"
# stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape\\shape-3-3_3-3-2_1.sto"

time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_50.tim"
core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_50.cor"
stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_50.sto"

TIME = get_TIME(time)

STG2_C, STG2_R = TIME[2][1], TIME[2][2]

CORE = get_CORE(core)

STOCH = get_STOCH(stoch)
STOCH, perturb = STOCH

numScens = length(STOCH)


FIRST_STG_COLS = []
var_type = "cont"
for i in CORE[2]
    # VA A SER SIEMPRE EL 2DO ELEMENTO LO QUE DIGA MARKER?
    if contains(i[2], "MARKER")
        var_type = "int"
    else
        if i[1] == STG2_C
            break
        end
        push!(FIRST_STG_COLS, [i[1], var_type])
    end
end

FIRST_STG_COLS = unique(FIRST_STG_COLS)

just_names = [i[1] for i in FIRST_STG_COLS]


FIRST_STG_OBJECT = []
object_name = ""

for i in CORE[1]
    if i[1] == "N"
        object_name = i[2]
    end
end
CORE[1] = CORE[1][2:end]

for i in CORE[2]
    if i[2] == object_name
        push!(FIRST_STG_OBJECT, [i[1], i[3]])
    elseif length(i) > 3 && i[4] == object_name
        push!(FIRST_STG_OBJECT, [i[1], i[5]])
    end
end

FIRST_STG_OBJECT = [[get(i[1]), get(i[2])] for i in FIRST_STG_OBJECT]

AUX = [i[1] for i in FIRST_STG_OBJECT]
for var in FIRST_STG_COLS
    if !(var[1] in AUX)
        push!(FIRST_STG_OBJECT, [var[1], 0.0])
    end
end
 
AUX = Dict(i[1] => i[2] for i in FIRST_STG_OBJECT)


FIRST_STG_ROWS = []
for i in CORE[1]
    if i[2] == STG2_R
        break
    end
    if i[2] != object_name
        push!(FIRST_STG_ROWS, i)
    end
end

FIRST_STG_CONSTR = Dict(i[2] => Dict() for i in FIRST_STG_ROWS)

FIRST_STG_RHS = Dict()

for row in FIRST_STG_ROWS
    comp, name = row

    for col in CORE[2]
        if name == col[2]
            FIRST_STG_CONSTR[name][col[1]] = get(col[3])
        end
        if length(col) > 3 && col[4] == name
            FIRST_STG_CONSTR[name][col[1]] = get(col[5])
        end
    end

    for rhs in CORE[3]
        if name == rhs[2] 
            FIRST_STG_RHS[name] = get(rhs[3])
        elseif length(rhs) > 3 && name == rhs[4]
            FIRST_STG_RHS[name] = get(rhs[5])
        end
    end

    
end


SEC_STG_COLS = []
flag = false
var_type = "cont"
for i in CORE[2]
    if contains(i[2], "MARKER")
        var_type = "int"
    else
        if i[1] == STG2_C
            flag = true
        end
        if flag
            push!(SEC_STG_COLS, [i[1], var_type])
        end
    end
end
println(1)
SEC_STG_COLS = unique(SEC_STG_COLS)

SEC_STG_OBJECT = []
for i in CORE[2]
    if i[2] == object_name
        push!(SEC_STG_OBJECT, [i[1], i[3]])
    elseif length(i) > 3 && i[4] == object_name
        push!(SEC_STG_OBJECT, [i[1], i[5]])
    end
end
println(2)
SEC_STG_OBJECT = [[get(i[1]), get(i[2])] for i in SEC_STG_OBJECT]


AUX2 = [i[1] for i in SEC_STG_OBJECT]
for var in SEC_STG_COLS
    if !(var[1] in AUX2)
        push!(SEC_STG_OBJECT, [var[1], 0.0])
    end
end
AUX2 = Dict((i[1], STOCH[s][1][2]) => i[2] for i in SEC_STG_OBJECT for s=1:numScens)

println(3)

SEC_STG_ROWS = [] 
flag = 0
for i in CORE[1]
    if i[2] == STG2_R
        flag = 1
    end
    if i[1] != "N" && flag != 0
        push!(SEC_STG_ROWS, i)
    end
end



SEC_STG_CONSTR = Dict(i[2] => Dict() for i in SEC_STG_ROWS)

SEC_STG_RHS = Dict()
println(4)

for row in SEC_STG_ROWS
    comp, name = row

    for col in CORE[2]
        if name == col[2]
            SEC_STG_CONSTR[name][col[1]] = get(col[3])
        end
            if length(col) > 3 && col[4] == name
            SEC_STG_CONSTR[name][col[1]] = get(col[5])
        end
    end

    for rhs in CORE[3]
        if name == rhs[2] 
            SEC_STG_RHS[name] = get(rhs[3])
        elseif length(rhs) > 3 && name == rhs[4]
            SEC_STG_RHS[name] = get(rhs[5])
        end
    end
end

for row in SEC_STG_ROWS
    comp, name = row
    if !(haskey(SEC_STG_RHS, name))
        SEC_STG_RHS[name] = 0.0
    end
end


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

SCENS = []

for s in 1:numScens
    # CREAR d, q, a
    copy_RHS = deepcopy(SEC_STG_RHS)
    copy_AUX2 = deepcopy(AUX2)
    copy_SEC_STG_CONSTR = deepcopy(SEC_STG_CONSTR)
    
    for change in STOCH[s][2]
        # ["RHS1", "R0000021", 8.1918]
        if contains(change[2], "RHS")
            copy_RHS[change[1]] = get(change[3])
        elseif contains(change[2], "obj")
            copy_AUX2[(change[1], STOCH[s][1][2])] = get(change[3])
        elseif contains(change[1], "RHS")
            copy_RHS[change[2]] = get(change[3])
        elseif contains(change[1], "obj")
            copy_AUX2[(change[2], STOCH[s][1][2])] = get(change[3])
        else
            println("Change not made: ", change)
        end
    end

    
    d = [copy_RHS[name[2]] for name in SEC_STG_ROWS]
    q = [copy_AUX2[(name[1], STOCH[s][1][2])] for name in SEC_STG_COLS]
    a = [[copy_SEC_STG_CONSTR[name[2]][col[1]] for col in vcat(FIRST_STG_COLS, SEC_STG_COLS)] for name in SEC_STG_ROWS]
    a = hcat(a...)'
    comps = [row[1] for row in SEC_STG_ROWS]
    s = Stoch_Scenario(get(STOCH[s][1][4]), d, q, a, comps, STOCH[s][1][2])
    push!(SCENS, s)
end

println(6)
sp = StochasticProgram([s for s in SCENS])


x = Dict()
non_neg = Dict(i[1] => true for i=FIRST_STG_COLS)

CONTS_IDX = []
INTS_IDX = []
BINS_IDX = []


# for bound in CORE[4]
#     if bound[3] in just_names
#         if bound[1] == "UP" || bound[1] == "UI"
#             @constraint(model, get_variable(bound[3], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) <= get(bound[4]))
#         elseif bound[1] == "FX"
#             @constraint(model, get_variable(bound[3], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) == get(bound[4]))
#         elseif bound[1] == "LO"
#             @constraint(model, get_variable(bound[3], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) >= get(bound[4]))
#             non_neg[bound[3]] = false
#         elseif bound[1] == "FR"
#             non_neg[bound[3]] = false
#         end
#     end
# end


for (name, type_var) in FIRST_STG_COLS
    if type_var == "cont"
        push!(CONTS_IDX, name)
    elseif type_var == "int"
        # revisar si tiene UP o UI 1. Si no, INTS_IDX. Si es que sí, BINS_IDX.
        added = false
        for bound in CORE[4]
            if bound[3] == name && !added
                # println(bound)
                if (bound[1] == "UP" || bound[1] == "UI") && (get(bound[4]) == 1 || get(bound[4]) == 1.0)
                    push!(BINS_IDX, name)
                    added = true
                end
            end
        end
        if !added
            push!(INTS_IDX, name)
        end
    else
        push!(BINS_IDX, name)
    end
end

println("BINS_IDX: ", BINS_IDX)
println("[INITIATING PREPROCESSING]")

# tenemos BINS_IDX: ver que eso sea igual a linking_vars
LINKING_VARS = []

for (comp, row) in SEC_STG_ROWS
    for (name, var_type) in FIRST_STG_COLS
        if haskey(SEC_STG_CONSTR[row], name) && !(name in LINKING_VARS)
            push!(LINKING_VARS, name)
        end
    end
end

# LINKING_VARS = unique(LINKING_VARS)
println("LINKING VARS: ", LINKING_VARS)

# for link_var in LINKING_VARS
#     if !(link_var in BINS_IDX)
#         throw("Problem not suited for Integer L-Shaped Method: linking variables are not all binary")
#     end
# end

# Fixed value for now
L = -1000000

println("[PREPROCESSING READY, L set to ", L, "]")



println(7)

@first_stage sp = begin

    @variable(model, conts[name=CONTS_IDX])
    # @variable(model, ints[name=INTS_IDX], category=:Int)
    # @variable(model, bins[name=BINS_IDX], category=:Bin)
    @variable(model, ints[name=INTS_IDX])
    @variable(model, bins[name=BINS_IDX])

    @objective(model, Min, sum(AUX[i[1]] * get_variable(i[1], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) for i=FIRST_STG_COLS))
    

    for (constr_type, name) in FIRST_STG_ROWS
        if constr_type == "E"
            @constraint(model, sum(val * get_variable(var, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) for (var, val) in FIRST_STG_CONSTR[name]) == FIRST_STG_RHS[name])
        elseif constr_type == "G"
            @constraint(model, sum(val * get_variable(var, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) for (var, val) in FIRST_STG_CONSTR[name]) >= FIRST_STG_RHS[name])
        else
            @constraint(model, sum(val * get_variable(var, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) for (var, val) in FIRST_STG_CONSTR[name]) <= FIRST_STG_RHS[name])
        end
    end
    
    just_names = [i[1] for i in FIRST_STG_COLS]

    for bound in CORE[4]
        if bound[3] in just_names
            if bound[1] == "UP" || bound[1] == "UI"
                @constraint(model, get_variable(bound[3], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) <= get(bound[4]))
            elseif bound[1] == "FX"
                @constraint(model, get_variable(bound[3], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) == get(bound[4]))
            elseif bound[1] == "LO"
                @constraint(model, get_variable(bound[3], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) >= get(bound[4]))
                non_neg[bound[3]] = false
            elseif bound[1] == "FR"
                non_neg[bound[3]] = false
            end
        end
    end

    for (key, val) in non_neg
        if val
            @constraint(model, get_variable(key, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins) >= 0)
        end
    end

    # println(model)
end


SEC_CONTS_IDX = []
SEC_INTS_IDX = []
SEC_BINS_IDX = []

for (name, type_var) in SEC_STG_COLS
    if type_var == "cont"
        push!(SEC_CONTS_IDX, name)
    elseif type_var == "int"
        push!(SEC_INTS_IDX, name)
    else
        push!(SEC_BINS_IDX, name)
    end
end


non_neg2 = Dict((i[1], STOCH[s][1][2]) => true for i=SEC_STG_COLS for s=1:numScens)


@second_stage sp = begin

    @decision conts ints bins
    s = scenario


    @variable(model, sec_conts[name=SEC_CONTS_IDX])
    # @variable(model, sec_ints[name=SEC_INTS_IDX], category=:Int)
    # @variable(model, sec_bins[name=SEC_BINS_IDX], category=:Bin)
    @variable(model, sec_ints[name=SEC_INTS_IDX])
    @variable(model, sec_bins[name=SEC_BINS_IDX])


    @objective(model, Min, sum(s.q[i] * sec_conts[SEC_STG_COLS[i][1]] for i=1:length(SEC_STG_COLS) if SEC_STG_COLS[i][1] in SEC_CONTS_IDX) + 
    sum(s.q[i] * sec_ints[SEC_STG_COLS[i][1]] for i=1:length(SEC_STG_COLS) if SEC_STG_COLS[i][1] in SEC_INTS_IDX) + 
    sum(s.q[i] * sec_bins[SEC_STG_COLS[i][1]] for i=1:length(SEC_STG_COLS) if SEC_STG_COLS[i][1] in SEC_BINS_IDX))


    for x in 1:length(SEC_STG_ROWS)

        if s.comps[x] == "E"
            @constraint(model, sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_conts[SEC_STG_COLS[y]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_CONTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_ints[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_INTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_bins[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_BINS_IDX) + 
            sum(s.a[x, y] * conts[FIRST_STG_COLS[y][1]] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in CONTS_IDX) +
            sum(s.a[x, y] * ints[FIRST_STG_COLS[y][1]] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in INTS_IDX) +
            sum(s.a[x, y] * bins[FIRST_STG_COLS[y][1]] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in BINS_IDX) == s.d[x])
        elseif s.comps[x] == "G"
            @constraint(model, sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_conts[SEC_STG_COLS[y]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_CONTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_ints[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_INTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_bins[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_BINS_IDX) +
            sum(s.a[x, y] * conts[FIRST_STG_COLS[y][1]] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in CONTS_IDX) +
            sum(s.a[x, y] * ints[FIRST_STG_COLS[y][1]] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in INTS_IDX) +
            sum(s.a[x, y] * bins[FIRST_STG_COLS[y][1]] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in BINS_IDX) >= s.d[x])
        else
            @constraint(model, sum(s.a[x, y] * sec_conts[SEC_STG_COLS[y]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_CONTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_ints[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_INTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_bins[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_BINS_IDX) +
            sum(s.a[x, y] * conts[FIRST_STG_COLS[y][1]] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in CONTS_IDX) +
            sum(s.a[x, y] * ints[FIRST_STG_COLS[y][1]] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in INTS_IDX) +
            sum(s.a[x, y] * bins[FIRST_STG_COLS[y][1]] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in BINS_IDX) <= s.d[x])
        end
    end

    just_names2 = [i[1] for i in SEC_STG_COLS]

    for bound in CORE[4]
        if bound[3] in just_names2
            if bound[1] == "UP" || bound[1] == "UI"
                if bound[3] in SEC_CONTS_IDX
                    @constraint(model, sec_conts[bound[3]] <= get(bound[4]))
                elseif bound[3] in SEC_INTS_IDX
                    @constraint(model, sec_ints[bound[3]] <= get(bound[4]))
                else
                    @constraint(model, sec_bins[bound[3]] <= get(bound[4]))
                end
            elseif bound[1] == "FX"
                if bound[3] in SEC_CONTS_IDX
                    @constraint(model, sec_conts[bound[3]] == get(bound[4]))
                elseif bound[3] in SEC_INTS_IDX
                    @constraint(model, sec_ints[bound[3]] == get(bound[4]))
                else
                    @constraint(model, sec_bins[bound[3]] == get(bound[4]))
                end
            elseif bound[1] == "LO"
                if bound[3] in SEC_CONTS_IDX
                    @constraint(model, sec_conts[bound[3]] >= get(bound[4]))
                elseif bound[3] in SEC_INTS_IDX
                    @constraint(model, sec_ints[bound[3]] >= get(bound[4]))
                else
                    @constraint(model, sec_bins[bound[3]] >= get(bound[4]))
                end
                non_neg2[(bound[3], s.name)] = false
            elseif bound[1] == "FR"
                non_neg2[(bound[3], s.name)] = false
            end
        end
    end

    for (key, val) in non_neg2
        if val
            if key[1] in SEC_CONTS_IDX
                @constraint(model, sec_conts[key[1]] >= 0)
            elseif key[1] in SEC_INTS_IDX
                @constraint(model, sec_ints[key[1]] >= 0)
            else
                @constraint(model, sec_bins[key[1]] >= 0)
            end
        end
    end    

    
end

# println(sp)
println(8)

# solve(sp, solver=GurobiSolver())

# multiple cuts :ls
# regularized decomposition :rd
# trust region :tr
# level sets :lv

# println(sp)
solve(sp, solver=LShapedSolver(:ls, GurobiSolver(OutputFlag=0)))


println("\nOptimal value: ", optimal_value(sp))

for i in 1:sp.numCols
    if sp.colVal[i] != 0.0
        println(i, " ", sp.colNames[i], ": ", sp.colVal[i])
    end
end
