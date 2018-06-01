using StructJuMP, Gurobi


function var_getter(STG2_C, x, y, i)
    if i >= STG2_C
        return y[i]
    else
        return x[i]
    end
end


include("smps_parser.jl")

time = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\semi.time"
core = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\semi.core"
stoch = "C:\\Jose\\Universidad\\JULIA_MISTI\\SMPS_Parser\\semi2.stoch"

TIME = get_TIME(time)
CORE = get_CORE(core)
STOCH = get_STOCH(stoch)

numScens = length(STOCH)

m = StructuredModel(num_scenarios=numScens);

# CORE = [ROWS, COLUMNS, RHS, BOUNDS]

STG2_C, STG2_R = TIME[2][1], TIME[2][2]

#  w = (((z' invSigma * z) - .2 * (uno' * invSigma * z))/((uno' * invSigma * 1) * (z' * invSigma * z) - (z' * invSigma * uno)^2)) * (invSigma * uno)


FIRST_STG_COLS = []
for i in CORE[2]
    if i[1] == STG2_C
        break
    end
    push!(FIRST_STG_COLS, i[1])
end

FIRST_STG_COLS = unique(FIRST_STG_COLS)

# HERE I'M NOT YET CONSIDERING BOUNDS: PLEASE DO NOT FORGET
@variable(m, x[i = FIRST_STG_COLS])

# OBTENER FUNCIÓN OBJETIVO: PARA CADA LÍNEA DE COLUMNS VER SI SALE OBJECTRW EN i[2] o i[4] (revisar length antes) (1er caso meter i[3], en el 2do i[5])
# Las que no aparezcan deben ser un 0, y de ahí crear el producto punto entre x y estos valores, y minimizar

FIRST_STG_OBJECT = []
for i in CORE[2]
    if i[2] == "OBJECTRW" && i[1] < STG2_C
        push!(FIRST_STG_OBJECT, [i[1], i[3]])
    elseif length(i) > 3 && i[4] == "OBJECTRW" && i[1] < STG2_C
        push!(FIRST_STG_OBJECT, [i[1], i[3]])
    end
end

FIRST_STG_OBJECT = [[get(i[1]), get(i[2])] for i in FIRST_STG_OBJECT]
# println(FIRST_STG_OBJECT)



AUX = [i[1] for i in FIRST_STG_OBJECT]
for var in FIRST_STG_COLS
    if !(var in AUX)
        push!(FIRST_STG_OBJECT, [var, 0.0])
    end
end
 
AUX = Dict(i[1] => i[2] for i in FIRST_STG_OBJECT)

@objective(m, :Min, sum(AUX[i] * x[i] for i=FIRST_STG_COLS));


FIRST_STG_ROWS = []
for i in CORE[1]
    if i[2] != "OBJECTRW" && i[2] < STG2_R
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
        elseif length(col) > 3 && col[4] == name
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

for bound in CORE[4]
    if bound[3] < STG2_C
        if bound[1] == "UI" || bound[1] == "UP"
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
for i in CORE[2]
    if i[1] == STG2_C
        flag = true
    end
    if flag
        push!(SEC_STG_COLS, i[1])
    end
end
println(1)
SEC_STG_COLS = unique(SEC_STG_COLS)

SEC_STG_OBJECT = []
for i in CORE[2]
    if i[2] == "OBJECTRW" && i[1] >= STG2_C
        push!(SEC_STG_OBJECT, [i[1], i[3]])
    elseif length(i) > 3 && i[4] == "OBJECTRW" && i[1] >= STG2_C
        push!(SEC_STG_OBJECT, [i[1], i[3]])
    end
end
println(2)
SEC_STG_OBJECT = [[get(i[1]), get(i[2])] for i in SEC_STG_OBJECT]


AUX2 = [i[1] for i in SEC_STG_OBJECT]
for var in SEC_STG_COLS
    if !(var in AUX2)
        push!(SEC_STG_OBJECT, [var, 0.0])
    end
end
AUX2 = Dict(i[1] => i[2] for i in SEC_STG_OBJECT)
# println(sum(AUX2[i] for i=SEC_STG_COLS))
println(3)

SEC_STG_ROWS = [] #Van todas las restricciones en los escenarios
for i in CORE[1]
    if i[2] != "OBJECTRW"
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
        elseif length(col) > 3 && col[4] == name
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
# println(AUX2)
# println("\n\n")
# println(SEC_STG_COLS)
for s in 1:numScens
    sb = StructuredModel(parent=m, id = s, prob = get(STOCH[s][1][4]));

    # HERE I'M NOT YET CONSIDERING BOUNDS: PLEASE DO NOT FORGET
    @variable(sb, y[i = SEC_STG_COLS])

    # OBJECTIVE

    @objective(sb, :Min, sum(AUX2[i] * y[i] for i=SEC_STG_COLS))

    cur_RHS = copy(SEC_STG_RHS)

    for change in STOCH[s][2]

        # println(change)
        # ["RHS1", "R0000021", 8.1918]
        if contains(change[1], "RHS")
            cur_RHS[change[2]] = change[3]
            # Faltan casos que no sean RHS: por ahora no lo veo
        else
            println("No estoy haciendo el cambio porque no es RHS", change)
        end
        # aplicar el change, en una copia de las constraints, a la constraint correcta. Después del for, @constraint
    end

    # CONSTRAINTS
    for row in SEC_STG_ROWS
        comp, name = row
        if comp == "E"
            @constraint(sb, 
            sum(i[2] * var_getter(STG2_C, x, y, i[1]) for i in SEC_STG_CONSTR[name]) == SEC_STG_RHS[name])
        elseif comp == "L"
            @constraint(sb,
            sum(i[2] * var_getter(STG2_C, x, y, i[1]) for i in SEC_STG_CONSTR[name]) <= SEC_STG_RHS[name])
        else
            @constraint(sb,
            sum(i[2] * var_getter(STG2_C, x, y, i[1]) for i in SEC_STG_CONSTR[name]) >= SEC_STG_RHS[name])
        end
    end

    for bound in CORE[4]
        if bound[3] >= STG2_C
            if bound[1] == "UI" || bound[1] == "UP"
                @constraint(sb, y[bound[3]] <= get(bound[4]))
            elseif bound[1] == "FX"
                @constraint(sb, y[bound[3]] == get(bound[4]))
            elseif bound[1] == "LO"
                @constraint(sb, y[bound[3]] >= get(bound[4]))
            end
        end
    end


end


# status, objval, soln = DLP(m, GurobiSolver())