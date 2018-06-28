# using DSPsolver
using StructJuMP, Gurobi


function var_getter(STG2_C, x, y, i, FIRST_STG_COLS)
    just_names = [i[1] for i in FIRST_STG_COLS]
    if i in just_names
        return x[i]
    else
        return y[i]
    end
end


include("smps_parser.jl")

time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape-3-3_3-3-2_1.tim"
core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape-3-3_3-3-2_1.mps"
stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\shape-3-3_3-3-2_1.sto"

TIME = get_TIME(time)

CORE = get_CORE(core)

STOCH = get_STOCH(stoch)
STOCH, perturb = STOCH

numScens = length(STOCH)

m = StructuredModel(num_scenarios=numScens);

# CORE = [ROWS, COLUMNS, RHS, BOUNDS]

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

# @variable(m, x[i = FIRST_STG_COLS] >= 0, Bin)
just_names = [i[1] for i in FIRST_STG_COLS]

x = Dict()

for i in FIRST_STG_COLS
    if i[2] == "int"
        x[i[1]] = @variable(m, category=:Int)
    elseif i[2] == "bin"
        x[i[1]] = @variable(m, category=:Bin)
    else
        x[i[1]] = @variable(m)
    end
end


# if all(i[2] == "int" for i in FIRST_STG_COLS)
#     @variable(m, x[i = just_names], Int)
# else
#     @variable(m, x[i = just_names])
# end
# HACERLO BIEN: VER CADA UNA SU TIPO Y AGREGARLA. ADEMÁS, AGREGAR BOUNDS >= 0 POR DEFECTO EN PARTE DE BOUNDS

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


@objective(m, :Min, sum(AUX[i[1]] * x[i[1]] for i=FIRST_STG_COLS));


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
        sum(i[2] * x[i[1]] for i in FIRST_STG_CONSTR[name]) == FIRST_STG_RHS[name])
    elseif comp == "L"
        @constraint(m,
        sum(i[2] * x[i[1]] for i in FIRST_STG_CONSTR[name]) <= FIRST_STG_RHS[name])
    else
        @constraint(m,
        sum(i[2] * x[i[1]] for i in FIRST_STG_CONSTR[name]) >= FIRST_STG_RHS[name])
    end
    
end

# BOUNDS

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
            @constraint(m, x[bound[3]] <= get(bound[4]))
        elseif bound[1] == "FX"
            @constraint(m, x[bound[3]] == get(bound[4]))
        elseif bound[1] == "LO"
            @constraint(m, x[bound[3]] >= get(bound[4]))
        end
    end
end

println("Root scenario ready")

##### EMPIEZA SEGUNDA ETAPA #####


# [["SC", "SCEN0001", "ROOT", 0.75, "STG00002"], [CHANGES]]

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
AUX2 = Dict(i[1] => i[2] for i in SEC_STG_OBJECT)
# println(sum(AUX2[i] for i=SEC_STG_COLS))
println(3)

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

println("Before numScens")
for s in 1:numScens
    sb = StructuredModel(parent=m, id = s, prob = get(STOCH[s][1][4]));
    just_names2 = [i[1] for i in SEC_STG_COLS]
    
    y = Dict()
    for i in SEC_STG_COLS
        if i[2] == "int"
            y[i[1]] = @variable(sb, category=:Int)
        elseif i[2] == "bin"
            y[i[1]] = @variable(sb, category=:Bin)
        else
            y[i[1]] = @variable(sb)
        end
    end

    # @variable(sb, y[i = just_names2])

    cur_RHS = copy(SEC_STG_RHS)

    if perturb == "replace"
        for change in STOCH[s][2]
            # ["RHS1", "R0000021", 8.1918]
            if contains(change[1], "RHS")
                cur_RHS[change[2]] = change[3]
            elseif contains(change[2], "obj")
                # Faltan casos que no sean RHS u obj: por ahora no lo veo
                AUX2[change[1]] = get(change[3])
            else
                println("No estoy haciendo el cambio porque no es RHS u obj", change)
            end
            # aplicar el change, en una copia de las constraints, a la constraint correcta. Después del for, @constraint
        end
    elseif perturb == "add"
        println("ADD")
        for change in STOCH[s][2]
            # ["RHS1", "R0000021", 8.1918]
            if contains(change[1], "RHS")
                cur_RHS[change[2]] += change[3]
            elseif contains(change[2], "obj")
                # Faltan casos que no sean RHS u obj: por ahora no lo veo
                AUX2[change[1]] += get(change[3])
            else
                println("No estoy haciendo el cambio porque no es RHS u obj", change)
            end
            # aplicar el change, en una copia de las constraints, a la constraint correcta. Después del for, @constraint
        end
    else
        for change in STOCH[s][2]
            # ["RHS1", "R0000021", 8.1918]
            if contains(change[1], "RHS")
                cur_RHS[change[2]] *= change[3]
            elseif contains(change[2], "obj")
                # Faltan casos que no sean RHS u obj: por ahora no lo veo
                AUX2[change[1]] *= get(change[3])
            else
                println("No estoy haciendo el cambio porque no es RHS u obj", change)
            end
            # aplicar el change, en una copia de las constraints, a la constraint correcta. Después del for, @constraint
        end
    end
    @objective(sb, :Min, sum(AUX2[i[1]] * y[i[1]] for i=SEC_STG_COLS))
    # CONSTRAINTS
    for row in SEC_STG_ROWS
        comp, name = row
        if comp == "E"
            @constraint(sb, 
            sum(i[2] * var_getter(STG2_C, x, y, i[1], FIRST_STG_COLS) for i in SEC_STG_CONSTR[name]) == SEC_STG_RHS[name])
        elseif comp == "L"
            @constraint(sb,
            sum(i[2] * var_getter(STG2_C, x, y, i[1], FIRST_STG_COLS) for i in SEC_STG_CONSTR[name]) <= SEC_STG_RHS[name])
        else
            @constraint(sb,
            sum(i[2] * var_getter(STG2_C, x, y, i[1], FIRST_STG_COLS) for i in SEC_STG_CONSTR[name]) >= SEC_STG_RHS[name])
        end
    end

    for bound in CORE[4]
        if bound[3] in just_names2
            if bound[1] == "UI"
                @constraint(sb, x[bound[3]] <= get(bound[4]))
                for i in 1:length(SEC_STG_COLS)
                    if bound[3] == SEC_STG_COLS[i][1]
                        SEC_STG_COLS[i] = [SEC_STG_COLS[i][1], "int"]
                    end
                end
            elseif bound[1] == "UP"
                @constraint(sb, y[bound[3]] <= get(bound[4]))
            elseif bound[1] == "FX"
                @constraint(sb, y[bound[3]] == get(bound[4]))
            elseif bound[1] == "LO"
                @constraint(sb, y[bound[3]] >= get(bound[4]))
            end
        end
    end

end

println("m:\n")

println(m)

# status, objval, soln = DLP(m, GurobiSolver())



# DSPsolver.loadProblem(m);       # Load model m to DSP
# DSPsolver.solve(DSP_SOLVER_DD); # Solve problem using dual decomposition


# println("Upper Bound: ", DSPsolver.getPrimalBound());
# println("Lower Bound: ", DSPsolver.getDualBound());


