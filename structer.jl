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
# @variable(m, x[i = FIRST_STG_COLS] ¿>= 0?)


