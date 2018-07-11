
using JuMP, Gurobi


include("smps_parser.jl")

time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape-3-3_3-3-2_1.tim"
core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape-3-3_3-3-2_1.mps"
stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape-3-3_3-3-2_1.sto"

TIME = get_TIME(time)

CORE = get_CORE(core)

STOCH = get_STOCH(stoch)
STOCH, perturb = STOCH

numScens = length(STOCH)

m = Model(solver=GurobiSolver())

STG2_C, STG2_R = TIME[2][1], TIME[2][2]

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

x = Dict()
non_neg = Dict()

for i in FIRST_STG_COLS
    if i[2] == "int"
        x[(i[1], "MASTER")] = @variable(m, category=:Int)
    elseif i[2] == "bin"
        x[(i[1], "MASTER")] = @variable(m, category=:Bin)
    else
        x[(i[1], "MASTER")] = @variable(m)
    end
    non_neg[(i[1], "MASTER")] = true
end


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


# "SCEN0001" => 0.75
names_prob = Dict(scen[1][2] => get(scen[1][4]) for scen in STOCH)
names_prob["MASTER"] = 1.0

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
# println(1)
SEC_STG_COLS = unique(SEC_STG_COLS)

for (scen, prob) in names_prob
    for i in SEC_STG_COLS
        # println(i[1], " ", scen)
        if i[2] == "int"
            x[(i[1], scen)] = @variable(m, category=:Int)
        elseif i[2] == "bin"
            x[(i[1], scen)] = @variable(m, category=:Bin)
        else
            x[(i[1], scen)] = @variable(m)
        end
        non_neg[(i[1], scen)] = true
    end
end


SEC_STG_OBJECT = []
for i in CORE[2]
    if i[2] == object_name
        push!(SEC_STG_OBJECT, [i[1], i[3]])
    elseif length(i) > 3 && i[4] == object_name
        push!(SEC_STG_OBJECT, [i[1], i[5]])
    end
end
# println(2)
SEC_STG_OBJECT = Dict(get(i[1]) => get(i[2]) for i in SEC_STG_OBJECT)

# [["SC", "SCEN0001", "ROOT", 0.75, "STG00002"], [CHANGES]]
# CHANGES = [["q_00", "obj", 28], [...]]
for scen in STOCH
    COPY_SEC_STG_OBJECT = deepcopy(SEC_STG_OBJECT)
    for change in scen[2]
        if change[2] == "obj" && perturb == "replace"
           COPY_SEC_STG_OBJECT[change[1]] = get(change[3])
        elseif change[2] == "obj" && perturb == "add"
            COPY_SEC_STG_OBJECT[change[1]] += get(change[3])
        elseif change[2] == "obj" && perturb == "multiply"
            COPY_SEC_STG_OBJECT[change[1]] *= get(change[3])
        end
    end
    for (key, value) in COPY_SEC_STG_OBJECT
        AUX[(key, scen[1][2])] = value
    end
end



obj = sum(AUX[(i[1], "MASTER")] * x[(i[1], "MASTER")] for i=FIRST_STG_COLS) + sum(prob * AUX[(i[1], scen)] * x[(i[1], scen)] for i=SEC_STG_COLS for (scen, prob) in names_prob)

@objective(m, :Min, obj)

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

    if comp == "E"
        @constraint(m, 
        sum(i[2] * x[(i[1], "MASTER")] for i in FIRST_STG_CONSTR[name]) == FIRST_STG_RHS[name])
    elseif comp == "L"
        @constraint(m,
        sum(i[2] * x[(i[1], "MASTER")] for i in FIRST_STG_CONSTR[name]) <= FIRST_STG_RHS[name])
    else
        @constraint(m,
        sum(i[2] * x[(i[1], "MASTER")] for i in FIRST_STG_CONSTR[name]) >= FIRST_STG_RHS[name])
    end
end


SEC_STG_ROWS = [] #Van todas las restricciones en los escenarios
flag1 = 0
for i in CORE[1]
    if i[2] == STG2_R
        flag = 1
    end
    if i[2] != "OBJECTRW" && flag != 0
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
    comp, name = row

    for col in CORE[2]
        if name == col[2]
            CONSTR[name][col[1]] = get(col[3])
        end
        if length(col) > 3 && col[4] == name
            CONSTR[name][col[1]] = get(col[5])
        end
    end

    for rhs in CORE[3]
        if name == rhs[2] 
            RHS[name] = get(rhs[3])
        elseif length(rhs) > 3 && name == rhs[4]
            RHS[name] = get(rhs[5])
        end
    end

    if comp == "E"
        for (scen, prob) in names_prob
            @constraint(m, 
            sum(i[2] * x[(i[1], scen)] for i in CONSTR[name] if haskey(x[(i[1], scen)])) + 
            sum(i[2] * x[(i[1], "MASTER")] for i in CONSTR[name] if !haskey(x[(i[1], scen)])) == RHS[name])
        end
    elseif comp == "L"
        for (scen, prob) in names_prob
            @constraint(m, 
            sum(i[2] * x[(i[1], scen)] for i in CONSTR[name] if haskey(x[(i[1], scen)])) + 
            sum(i[2] * x[(i[1], "MASTER")] for i in CONSTR[name] if !haskey(x[(i[1], scen)])) <= RHS[name])
        end
    else
        for (scen, prob) in names_prob
            @constraint(m, 
            sum(i[2] * x[(i[1], scen)] for i in CONSTR[name] if haskey(x, (i[1], scen))) + 
            sum(i[2] * x[(i[1], "MASTER")] for i in CONSTR[name] if !haskey(x, (i[1], scen))) >= RHS[name])
        end
    end
    
end


# FALTAN BOUNDS. PERO ANTES CORREGIR QUE VARIABLES SEAN POR ESCENARIO

just_names = [i[1] for i in FIRST_STG_COLS]

for bound in CORE[4]
    if bound[3] in just_names
        if bound[1] == "UI"
            @constraint(m, x[bound[3]] <= get(bound[4]))
            for i in 1:length(FIRST_STG_COLS)
                if bound[3] == FIRST_STG_COLS[i][1]
                    FIRST_STG_COLS[i] = [FIRST_STG_COLS[i][1], "int"]
                end
            end
        elseif bound[1] == "UP"
            @constraint(m, x[(bound[3], "MASTER")] <= get(bound[4]))
        elseif bound[1] == "FX"
            @constraint(m, x[(bound[3], "MASTER")] == get(bound[4]))
        elseif bound[1] == "LO"
            @constraint(m, x[(bound[3], "MASTER")] >= get(bound[4]))
            non_neg[(bound[3], "MASTER")] = false
        elseif bound[1] == "FR"
            non_neg[(bound[3], "MASTER")] = false
        end
    end
end



just_names = [i[1] for i in SEC_STG_COLS]

for bound in CORE[4]
    if bound[3] in just_names
        if bound[1] == "UI"
            @constraint(m, x[bound[3]] <= get(bound[4]))
            for i in 1:length(SEC_STG_COLS)
                if bound[3] == SEC_STG_COLS[i][1]
                    SEC_STG_COLS[i] = [SEC_STG_COLS[i][1], "int"]
                end
            end
        elseif bound[1] == "UP"
            for (scen, prob) in names_prob
                @constraint(m, x[(bound[3], scen)] <= get(bound[4]))
            end
        elseif bound[1] == "FX"
            for (scen, prob) in names_prob
                @constraint(m, x[(bound[3], scen)] == get(bound[4]))
            end
        elseif bound[1] == "LO"
            for (scen, prob) in names_prob
                @constraint(m, x[(bound[3], scen)] >= get(bound[4]))
                non_neg[(bound[3], scen)] = false
            end
        elseif bound[1] == "FR"
            for (scen, prob) in names_prob
                non_neg[(bound[3], scen)] = false
            end
        end
    end
end


for var in FIRST_STG_COLS
    if non_neg[(var[1], "MASTER")]
        @constraint(m, x[(var[1], "MASTER")] >= 0)
    end
end

for var in SEC_STG_COLS
    for (scen, prob) in names_prob
        if non_neg[(var[1], scen)]
            @constraint(m, x[(var[1], scen)] >= 0)
        end
    end
end

# println(m)


solve(m)

for (key, value) in x
    println(key, ": ", getvalue(value))
end


