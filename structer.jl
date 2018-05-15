using StructJuMP, Gurobi

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
    if i[2] == "OBJECTRW"
        push!(FIRST_STG_OBJECT, [i[1], i[3]])
    elseif length(i) > 3 && i[4] == "OBJECTRW"
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
    if i[1] == STG2_R
        break
    end
    push!(FIRST_STG_ROWS, i)
end

for row in FIRST_STG_ROWS
    comp, name = row
    # add constraint indicated by this row: we need RHS and all columns which include this row first



end




