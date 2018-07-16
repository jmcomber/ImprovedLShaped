using StochasticPrograms, Gurobi


# struct Stoch_Scenario <: AbstractScenarioData
#     π::Float64
#     d::Vector{Float64}
#     q::Vector{Float64}
#     a:: Array{Float64,2}
#     comps::Array{String, 1}
#     name::String
# end

function get_variable(name, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins, scen)
    if name in CONTS_IDX
        return conts[name, scen]
    elseif name in INTS_IDX
        return ints[name, scen]
    else
        return bins[name, scen]
    end
end

StochasticPrograms.probability(s::Stoch_Scenario) = s.π

include("smps_parser.jl")

# time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape-3-3_3-3-2_1.tim"
# core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape-3-3_3-3-2_1.mps"
# stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape-3-3_3-3-2_1.sto"

time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_50.tim"
core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_50.cor"
stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\sslp\\sslp_5_25_50.sto"

TIME = get_TIME(time)

STG2_C, STG2_R = TIME[2][1], TIME[2][2]

CORE = get_CORE(core)

STOCH = get_STOCH(stoch)
STOCH, perturb = STOCH

numScens = length(STOCH)



# PRIMERO CREEMOS EL CORE, PARA DE AHI EN UN FOR MODIFICAR Y CREAR CADA ESCENARIO?

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

# x = Dict()
# non_neg = Dict()

# for i in FIRST_STG_COLS
#     if i[2] == "int"
#         x[(i[1], "MASTER")] = @variable(m, category=:Int)
#     elseif i[2] == "bin"
#         x[(i[1], "MASTER")] = @variable(m, category=:Bin)
#     else
#         x[(i[1], "MASTER")] = @variable(m)
#     end
#     non_neg[(i[1], "MASTER")] = true
# end


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
 
AUX = Dict((i[1], "MASTER") => i[2] for i in FIRST_STG_OBJECT)


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

    # if comp == "E"
    #     @constraint(m, 
    #     sum(i[2] * x[i[1]] for i in FIRST_STG_CONSTR[name]) == FIRST_STG_RHS[name])
    # elseif comp == "L"
    #     @constraint(m,
    #     sum(i[2] * x[i[1]] for i in FIRST_STG_CONSTR[name]) <= FIRST_STG_RHS[name])
    # else
    #     @constraint(m,
    #     sum(i[2] * x[i[1]] for i in FIRST_STG_CONSTR[name]) >= FIRST_STG_RHS[name])
    # end
    
end

# BOUNDS

# just_names = [i[1] for i in FIRST_STG_COLS]

# for bound in CORE[4]
#     if bound[3] in just_names
#         if bound[1] == "UI"
#             @constraint(m, x[bound[3]] <= get(bound[4]))
#             for i in 1:length(FIRST_STG_COLS)
#                 if bound[3] == FIRST_STG_COLS[i][1]
#                     FIRST_STG_COLS[i] = [FIRST_STG_COLS[i][1], "int"]
#                 end
#             end
#         elseif bound[1] == "UP"
#             @constraint(m, x[(bound[3], "MASTER")] <= get(bound[4]))
#         elseif bound[1] == "FX"
#             @constraint(m, x[(bound[3], "MASTER")] == get(bound[4]))
#         elseif bound[1] == "LO"
#             @constraint(m, x[(bound[3], "MASTER")] >= get(bound[4]))
#             non_neg[(bound[3], "MASTER")] = false
#         elseif bound[1] == "FR"
#             non_neg[(bound[3], "MASTER")] = false
#         end
#     end
# end


# CREAR ESCENARIOS

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
# println(sum(AUX2[i] for i=SEC_STG_COLS))
println(3)

SEC_STG_ROWS = [] #Van todas las restricciones en los escenarios
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


SCENS = []


# for i in SEC_STG_CONSTR
#     println(i)
# end

# names_prob = Dict(scen[1][2] => get(scen[1][4]) for scen in STOCH)


# d lado derecho, q obj, a coefs

for s in 1:numScens
    # CREAR d, q, a
    copy_RHS = deepcopy(SEC_STG_RHS)
    copy_AUX2 = deepcopy(AUX2)
    # println(copy_AUX2)
    copy_SEC_STG_CONSTR = deepcopy(SEC_STG_CONSTR)
    

    # println("Escenario es ", STOCH[s][1][2])
    for change in STOCH[s][2]
        # ["RHS1", "R0000021", 8.1918]
        # println("change: ", change)
        # println("copy_AUX2: ", copy_AUX2)
        if contains(change[2], "RHS")
            copy_RHS[change[1]] = get(change[3])
        elseif contains(change[2], "obj")
            # println(change[1], " ", STOCH[s][1][2])
            # println((change[1], STOCH[s][1][2]))
            copy_AUX2[(change[1], STOCH[s][1][2])] = get(change[3])
        elseif contains(change[1], "RHS")
            copy_RHS[change[2]] = get(change[3])
        elseif contains(change[1], "obj")
            copy_AUX2[(change[2], STOCH[s][1][2])] = get(change[3])
        else
            println("Change not made: ", change)
        end
    end

    # println(SEC_STG_CONSTR)

    # none = [[println(name[2], " ", col[1], " ", col[2]) for col in vcat(FIRST_STG_COLS, SEC_STG_COLS)] for name in SEC_STG_ROWS]

    d = [copy_RHS[name[2]] for name in SEC_STG_ROWS]
    q = [copy_AUX2[(name[1], STOCH[s][1][2])] for name in SEC_STG_COLS]
    # a = [[col[2] for col in copy_SEC_STG_CONSTR[name[2]]] for name in SEC_STG_ROWS]
    # println(SEC_STG_ROWS)
    a = [[copy_SEC_STG_CONSTR[name[2]][col[1]] for col in vcat(FIRST_STG_COLS, SEC_STG_COLS)] for name in SEC_STG_ROWS]
    # println(size(a))
    a = hcat(a...)'
    # println(size(a))
    comps = [row[1] for row in SEC_STG_ROWS]
    println(STOCH[s][1])
    s = Stoch_Scenario(get(STOCH[s][1][4]), d, q, a, comps, STOCH[s][1][2])
    push!(SCENS, s)
end
# println(SCENS)




function StochasticPrograms.expected(sds::Vector{Stoch_Scenario})
    sd = Stoch_Scenario(1, sum([s.π * s.d for s in sds]), sum([s.π * s.q for s in sds]))
end



sp = StochasticProgram([s for s in SCENS], solver=GurobiSolver())

# for (name, constr) in FIRST_STG_CONSTR
#     println(name, ": ", constr)
# end

# println(AUX)
x = Dict()
non_neg = Dict(i[1] => true for i=FIRST_STG_COLS)

CONTS_IDX = []
INTS_IDX = []
BINS_IDX = []

for (name, type_var) in FIRST_STG_COLS
    if type_var == "cont"
        push!(CONTS_IDX, name)
    elseif type_var == "int"
        push!(INTS_IDX, name)
    else
        push!(BINS_IDX, name)
    end
end



@first_stage sp = begin
    # for (name, type_var) in FIRST_STG_COLS
    #     if type_var == "cont"
    #        x[(name, "MASTER")] = @variable(model)
    #     elseif type_var == "int"
    #        x[(name, "MASTER")] = @variable(model, category=:Int)
    #     else
    #        x[(name, "MASTER")] = @variable(model, category=:Bin)
    #     end
    #     non_neg[(name, "MASTER")] = true
    # end

    @variable(model, conts[name=CONTS_IDX, scen=["MASTER"]])
    @variable(model, ints[name=INTS_IDX, scen=["MASTER"]], category=:Int)
    @variable(model, bins[name=BINS_IDX, scen=["MASTER"]], category=:Bin)

    # println(ints)
    @objective(model, Min, sum(AUX[(i[1], "MASTER")] * get_variable(i[1], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins, "MASTER") for i=FIRST_STG_COLS))
    

    for (constr_type, name) in FIRST_STG_ROWS
        if constr_type == "E"
            @constraint(model, sum(val * get_variable(var, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins, "MASTER") for (var, val) in FIRST_STG_CONSTR[name]) == FIRST_STG_RHS[name])
        elseif constr_type == "G"
            @constraint(model, sum(val * get_variable(var, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins, "MASTER") for (var, val) in FIRST_STG_CONSTR[name]) >= FIRST_STG_RHS[name])
        else
            @constraint(model, sum(val * get_variable(var, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins, "MASTER") for (var, val) in FIRST_STG_CONSTR[name]) <= FIRST_STG_RHS[name])
        end
    end
    
    # println(model)

    
    just_names = [i[1] for i in FIRST_STG_COLS]

    for bound in CORE[4]
        if bound[3] in just_names
            if bound[1] == "UP" || bound[1] == "UI"
                @constraint(model, get_variable(bound[3], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins, "MASTER") <= get(bound[4]))
            elseif bound[1] == "FX"
                @constraint(model, get_variable(bound[3], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins, "MASTER") == get(bound[4]))
            elseif bound[1] == "LO"
                @constraint(model, get_variable(bound[3], CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins, "MASTER") >= get(bound[4]))
                non_neg[(bound[3], "MASTER")] = false
            elseif bound[1] == "FR"
                non_neg[(bound[3], "MASTER")] = false
            end
        end
    end

    for (key, val) in non_neg
        if val
            @constraint(model, get_variable(key, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins, "MASTER") >= 0)
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

# for s = 1:numScens
#     println(STOCH[s][1])
# end

non_neg2 = Dict((i[1], STOCH[s][1][2]) => true for i=SEC_STG_COLS for s=1:numScens)


@second_stage sp = begin

    @decision conts ints bins
    s = scenario

    # println(s.comps)

    @variable(model, sec_conts[name=SEC_CONTS_IDX])
    @variable(model, sec_ints[name=SEC_INTS_IDX], category=:Int)
    @variable(model, sec_bins[name=SEC_BINS_IDX], category=:Bin)

    # println(SEC_INTS_IDX)
    # println(sec_ints)

    # println("Scenario is ", s.name)
    @objective(model, Min, sum(s.q[i] * sec_conts[SEC_STG_COLS[i][1]] for i=1:length(SEC_STG_COLS) if SEC_STG_COLS[i][1] in SEC_CONTS_IDX) + 
    sum(s.q[i] * sec_ints[SEC_STG_COLS[i][1]] for i=1:length(SEC_STG_COLS) if SEC_STG_COLS[i][1] in SEC_INTS_IDX) + 
    sum(s.q[i] * sec_bins[SEC_STG_COLS[i][1]] for i=1:length(SEC_STG_COLS) if SEC_STG_COLS[i][1] in SEC_BINS_IDX))


    for x in 1:length(SEC_STG_ROWS)
        # println(length(s.a))
        # println(length(FIRST_STG_COLS), " ", length(SEC_STG_COLS), " ", length(SEC_STG_ROWS))
        if s.comps[x] == "E"
            @constraint(model, sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_conts[SEC_STG_COLS[y]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_CONTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_ints[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_INTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_bins[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_BINS_IDX) + 
            sum(s.a[x, y] * conts[FIRST_STG_COLS[y][1], "MASTER"] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in CONTS_IDX) +
            sum(s.a[x, y] * ints[FIRST_STG_COLS[y][1], "MASTER"] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in INTS_IDX) +
            sum(s.a[x, y] * bins[FIRST_STG_COLS[y][1], "MASTER"] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in BINS_IDX) == s.d[x])
        elseif s.comps[x] == "G"
            @constraint(model, sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_conts[SEC_STG_COLS[y]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_CONTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_ints[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_INTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_bins[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_BINS_IDX) +
            sum(s.a[x, y] * conts[FIRST_STG_COLS[y][1], "MASTER"] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in CONTS_IDX) +
            sum(s.a[x, y] * ints[FIRST_STG_COLS[y][1], "MASTER"] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in INTS_IDX) +
            sum(s.a[x, y] * bins[FIRST_STG_COLS[y][1], "MASTER"] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in BINS_IDX) >= s.d[x])
        else
            @constraint(model, sum(s.a[x, y] * sec_conts[SEC_STG_COLS[y]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_CONTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_ints[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_INTS_IDX) + 
            sum(s.a[x, length(FIRST_STG_COLS) + y] * sec_bins[SEC_STG_COLS[y][1]] for y in 1:length(SEC_STG_COLS) if SEC_STG_COLS[y][1] in SEC_BINS_IDX) +
            sum(s.a[x, y] * conts[FIRST_STG_COLS[y][1], "MASTER"] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in CONTS_IDX) +
            sum(s.a[x, y] * ints[FIRST_STG_COLS[y][1], "MASTER"] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in INTS_IDX) +
            sum(s.a[x, y] * bins[FIRST_STG_COLS[y][1], "MASTER"] for y in 1:length(FIRST_STG_COLS) if FIRST_STG_COLS[y][1] in BINS_IDX) <= s.d[x])
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
            # @constraint(model, get_variable(key, CONTS_IDX, INTS_IDX, BINS_IDX, conts, ints, bins, "MASTER") >= 0)
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


# solve(sp, solver=LShapedSolver(:ls, GurobiSolver()))

solve(sp, solver=GurobiSolver())



