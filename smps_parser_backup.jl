include("structs.jl")


# Vector{Union{Float64, String}}

function get_TIME(file_name::String)
    f = open(file_name) # En verdad no te van a dar full path: obtenerlo
    NAME = split(readline(f), " ")[end]

    PERIODS = []
    l = readline(f)
    l = readline(f)
    while !contains(l, "ENDATA")
        aux = split(l, " ")
        aux2 = String[]
        for elem in aux
            if length(elem) > 0
                push!(aux2, elem)
            end
        end
        push!(PERIODS, aux2)
        l = readline(f)
    end
    close(f)
    return PERIODS
end

function get_CORE(file_name::String)
    f = open(file_name)
    l = readline(f)
    while '*' == l[1]
        l = readline(f)
    end
    l = split(readline(f), " ")

    aux = String[]
    for elem in l
        if length(elem) > 0
            push!(aux, elem)
        end
    end
    NAME_CORE = aux[end]

    # if NAME_CORE != NAME
    #     # ERROR
    #     println("ERROR: Name mismatch in CORE, name is ", NAME_CORE, "\n\n")
    # end

    ROWS = []
    l = readline(f)
    while !contains(l, "COLUMNS")
        aux = split(l, " ")
        aux2 = []
        for elem in aux
            if length(elem) > 0
                if !isnull(tryparse(Float64, elem))
                    push!(aux2, tryparse(Float64, elem))
                else
                    push!(aux2, elem)
                end
            end
        end
        push!(ROWS, aux2)
        l = readline(f)
    end


    COLUMNS = []
    l = readline(f)
    while !contains(l, "RHS")
        # if (contains(l, "MARKER"))
        #     l = readline(f)
        #     continue
        # end
        aux = split(l, " ")
        aux2 = []
        for elem in aux
            if length(elem) > 0
                if !isnull(tryparse(Float64, elem))
                    push!(aux2, tryparse(Float64, elem))
                else
                    push!(aux2, elem)
                end
            end
        end
        push!(COLUMNS, aux2)
        l = readline(f)
    end

    RHS = []
    l = readline(f)
    while !contains(l, "BOUNDS")
        aux = split(l, " ")
        aux2 = []
        for elem in aux
            if length(elem) > 0
                if !isnull(tryparse(Float64, elem))
                    push!(aux2, tryparse(Float64, elem))
                else
                    push!(aux2, elem)
                end
            end
        end
        push!(RHS, aux2)
        l = readline(f)
    end

    BOUNDS = []
    l = readline(f)
    while !contains(l, "ENDATA")
        aux = split(l, " ")
        aux2 = []
        for elem in aux
            if length(elem) > 0
                if !isnull(tryparse(Float64, elem))
                    push!(aux2, tryparse(Float64, elem))
                else
                    push!(aux2, elem)
                end
            end
        end
        push!(BOUNDS, aux2)
        l = readline(f)
    end

    close(f)
    return [ROWS, COLUMNS, RHS, BOUNDS]
end

function get_STOCH(file_name::String)
    f = open(file_name)
    l = readline(f)
    while !contains(l, "STOCH")
        l = readline(f)
    end
    NAME_STOCH = split(l, " ")[end]

    # if NAME_STOCH != NAME
    #     # ERROR
    #     println("ERROR: Name mismatch in STOCH\n\n")
    # end

    SCENARIOS = []
    l = readline(f)
    l = split(l, " ")
    perturb = "replace"
    if length(l) == 3 && l[3] == "ADD"
        perturb = "add"
    elseif length(l) == 3 && l[3] == "MULTIPLY"
        perturb = "multiply"
    end
    l = readline(f)
    cont = 0
    while !contains(l, "ENDATA")
        l = split(l, r"\t| ")
        CUR_SCEN = []
        
        for elem in l
            if length(elem) > 0
                if !isnull(tryparse(Float64, elem))
                    push!(CUR_SCEN, tryparse(Float64, elem))
                else
                    push!(CUR_SCEN, elem)
                end
            end
        end
        l = readline(f)
        
        CUR_DATA = []
        c = 0
        while !contains(l, "SC") & !contains(l, "ENDATA")
            aux = []
            aux2 = split(l, " ")
            for elem in aux2
                if length(elem) > 0
                    if !isnull(tryparse(Float64, elem))
                        push!(aux, tryparse(Float64, elem))
                    else
                        push!(aux, elem)
                    end
                end
            end
            push!(CUR_DATA, aux)
            l = readline(f)
        end
        push!(SCENARIOS, [CUR_SCEN, CUR_DATA])
    end

    close(f)

    return SCENARIOS, perturb

end


function get_first_stg_cols(CORE, STG2_C::String)
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
end

function get_first_stg_object(CORE)
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
    return FIRST_STG_OBJECT, object_name
end

function get_first_stg_rows(CORE, STG2_R, object_name::SubString{String})
    FIRST_STG_ROWS = []
    for i in CORE[1]
        if i[2] == STG2_R
            break
        end
        if i[2] != object_name
            push!(FIRST_STG_ROWS, convert(Array{String,1}, i))
        end
    end
    return FIRST_STG_ROWS
end

function get_first_stg_rhs(CORE, FIRST_STG_ROWS, FIRST_STG_CONSTR)
    FIRST_STG_RHS = Dict{String,Float64}()

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
    return FIRST_STG_RHS, FIRST_STG_CONSTR
end

function get_sec_stg_cols(CORE, STG2_C)
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
end

function get_sec_stg_object(CORE, object_name)
    SEC_STG_OBJECT = []
    for i in CORE[2]
        if i[2] == object_name
            push!(SEC_STG_OBJECT, [i[1], i[3]])
        elseif length(i) > 3 && i[4] == object_name
            push!(SEC_STG_OBJECT, [i[1], i[5]])
        end
    end
    SEC_STG_OBJECT = [[get(i[1]), get(i[2])] for i in SEC_STG_OBJECT]
    return SEC_STG_OBJECT
end

function get_sec_stg_rows(CORE, STG2_R)
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
    return SEC_STG_ROWS
end

function get_sec_stg_rhs(CORE, SEC_STG_ROWS, SEC_STG_CONSTR)
    SEC_STG_RHS = Dict()
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
    return SEC_STG_RHS, SEC_STG_CONSTR
end

function create_scenarios(numScens::Int64, STOCH, SEC_STG_RHS, SEC_STG_CONSTR, SEC_STG_ROWS, SEC_STG_COLS::Array{Array{String,1},1}, FIRST_STG_COLS::Array{Array{String,1},1}, AUX2)
    SCENS = Stoch_Scenario[]
    for s in 1:numScens
        copy_RHS = copy(SEC_STG_RHS)
        copy_AUX2 = copy(AUX2)
        copy_SEC_STG_CONSTR = copy(SEC_STG_CONSTR)
        
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

        h = [copy_RHS[name[2]] for name in SEC_STG_ROWS]
        q = [copy_AUX2[(name[1], STOCH[s][1][2])] for name in SEC_STG_COLS]
        T = [[copy_SEC_STG_CONSTR[name[2]][col[1]] for col in FIRST_STG_COLS] for name in SEC_STG_ROWS]
        T = hcat(T...)'
        W = [[copy_SEC_STG_CONSTR[name[2]][col[1]] for col in SEC_STG_COLS] for name in SEC_STG_ROWS]
        W = hcat(W...)'
        comps = [row[1] for row in SEC_STG_ROWS]
        s = Stoch_Scenario(get(STOCH[s][1][4]), h, q, T, W, comps)
        push!(SCENS, s)
    end
    return SCENS
end



function linking_vars(SEC_STG_ROWS, FIRST_STG_COLS::Array{Array{String,1},1}, SEC_STG_CONSTR)
    LINKING_VARS = []
    for (comp, row) in SEC_STG_ROWS
        for (name, var_type) in FIRST_STG_COLS
            if haskey(SEC_STG_CONSTR[row], name) && SEC_STG_CONSTR[row][name] != 0 && !(name in LINKING_VARS)
                push!(LINKING_VARS, name)
            end
        end
    end
    return LINKING_VARS
end

function create_Lk(master_data::StageData, SEC_STG_COLS::Array{Array{String,1},1}, SCENARIO::Stoch_Scenario, BOUNDS, multi_thread::Bool)
    if multi_thread
        m = Model(solver=GurobiSolver(OutputFlag=0))
    else
        m = Model(solver=GurobiSolver(OutputFlag=0, Threads=1))
    end
    names1 = [var[1] for var in master_data.cols]
    # Faltan tipos (usar setcategory(x[i], :Int))
    @variable(m, z[i = names1])
    # for (name, var_type) in FIRST_STG_COLS
    #     if var_type == "int"
    #         setcategory(z[name], :Int)
    #     elseif var_type == "bin"
    #         setcategory(z[name], :Bin)
    #     end
    # end
    names2 = [var[1] for var in SEC_STG_COLS]
    @variable(m, y[j = names2])
    # for (name, var_type) in SEC_STG_COLS
    #     if var_type == "int"
    #         setcategory(y[name], :Int)
    #     elseif var_type == "bin"
    #         setcategory(y[name], :Bin)
    #     end
    # end

    @objective(m, :Min, sum(SCENARIO.q[j] * y[names2[j]] for j in 1:length(names2)))

    # agregar constraints Az ~ b
    for (constr_type, name) in master_data.rows
        if constr_type == "E"
            @constraint(m, sum(val * z[var] for (var, val) in master_data.constrs[name]) == master_data.rhs[name])
        elseif constr_type == "G"
            @constraint(m, sum(val * z[var] for (var, val) in master_data.constrs[name]) >= master_data.rhs[name])
        else
            @constraint(m, sum(val * z[var] for (var, val) in master_data.constrs[name]) <= master_data.rhs[name])
        end
    end
    
    # Agregar constraints Tz + Wy ~ h
    for r = 1:size(SCENARIO.T, 1)
        if SCENARIO.comps[r] == "E"
            @constraint(m, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) == SCENARIO.h[r])
        elseif SCENARIO.comps[r] == "G"
            @constraint(m, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) >= SCENARIO.h[r])
        else 
            @constraint(m, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) <= SCENARIO.h[r])
        end
    end

    add_bounds!(m, names1, z, BOUNDS, master_data.cols)
    add_bounds!(m, names2, y, BOUNDS, SEC_STG_COLS)
    

    return m
end

function init_master(master_data::StageData, solver, BOUNDS, multi_thread::Bool)
    if multi_thread
        master = Model(solver=solver())
    else
        master = Model(solver=solver(Threads=1))
    end

    names1 = [var[1] for var in master_data.cols]
    @variable(master, x[i = names1])
    # for (name, var_type) in FIRST_STG_COLS
    #     if var_type == "int"
    #         setcategory(x[name], :Int)
    #     elseif var_type == "bin"
    #         setcategory(x[name], :Bin)
    #     end
    # end
    @variable(master, θ >= master_data.L)
    master.colNames = [names1..., "θ"]
    @objective(master, :Min, sum(master_data.obj[name] * x[name] for name in names1) + θ)

    # agregar constraints Ax ~ b
    for (constr_type, name) in master_data.rows
        if constr_type == "E"
            @constraint(master, sum(val * x[var] for (var, val) in master_data.constrs[name]) == master_data.rhs[name])
        elseif constr_type == "G"
            @constraint(master, sum(val * x[var] for (var, val) in master_data.constrs[name]) >= master_data.rhs[name])
        else
            @constraint(master, sum(val * x[var] for (var, val) in master_data.constrs[name]) <= master_data.rhs[name])
        end
    end
    add_bounds!(master, names1, x, BOUNDS, master_data.cols)
    return master, x, θ
end

function add_bounds!(m::JuMP.Model, names::Array{String,1}, z, BOUNDS, COLS)
    non_neg = Dict(i[1] => true for i=COLS)
    for bound in BOUNDS
        if bound[3] in names
            if bound[1] == "UP" || bound[1] == "UI"
                @constraint(m, z[bound[3]] <= get(bound[4]))
            elseif bound[1] == "FX"
                @constraint(m, z[bound[3]] == get(bound[4]))
            elseif bound[1] == "LO"
                @constraint(m, z[bound[3]] >= get(bound[4]))
                non_neg[bound[3]] = false
            elseif bound[1] == "FR"
                non_neg[bound[3]] = false
            end
        end
    end
    for (key, val) in non_neg
        if val
            @constraint(m, z[key] >= 0)
        end
    end
end

function create_v_xs(x_hat::Array{Float64,1}, SCENS::Array{Stoch_Scenario,1}, BOUNDS, FIRST_STG_COLS::Array{Array{String,1},1}, SEC_STG_COLS::Array{Array{String,1},1}, multi_thread::Bool)
    v_xs = Subproblem[]
    for i in 1:length(SCENS)
        prob = create_v_x(x_hat, SCENS[i], BOUNDS, FIRST_STG_COLS, SEC_STG_COLS, multi_thread)
        push!(v_xs, prob)
    end
    v_xs
end

function create_v_x(x_hat::Array{Float64,1}, SCENARIO::Stoch_Scenario, BOUNDS, FIRST_STG_COLS::Array{Array{String,1},1}, SEC_STG_COLS::Array{Array{String,1},1}, multi_thread::Bool)
    if multi_thread
        v_x = Model(solver=GurobiSolver(OutputFlag=0, InfUnbdInfo=1))
    else
        v_x = Model(solver=GurobiSolver(OutputFlag=0, InfUnbdInfo=1, Threads=1))
    end

    names1 = [var[1] for var in FIRST_STG_COLS]
    @variable(v_x, z[i = names1])
    names2 = [var[1] for var in SEC_STG_COLS]
    @variable(v_x, y[j = names2])

    # for (name, var_type) in SEC_STG_COLS
    #     if var_type == "int"
    #         setcategory(y[name], :Int)
    #     elseif var_type == "bin"
    #         setcategory(y[name], :Bin)
    #     end
    # end

    @objective(v_x, :Min, sum(SCENARIO.q[j] * y[names2[j]] for j in 1:length(names2)))

    # Tz + Wy ~ h
    # Its duals are λ: not necessary to get reference
    constr_λ = JuMP.ConstraintRef[]
    # println("size es ", size(SCENARIO.W, 2))
    for r = 1:size(SCENARIO.T, 1)
        if SCENARIO.comps[r] == "E"
            push!(constr_λ, @constraint(v_x, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) == SCENARIO.h[r]))
        elseif SCENARIO.comps[r] == "G"
            push!(constr_λ, @constraint(v_x, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) >= SCENARIO.h[r]))
        else 
            push!(constr_λ, @constraint(v_x, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) <= SCENARIO.h[r]))
        end
    end

    # z = x_hat
    # Its duals are π: get reference to return v_x, constr_π
    @constraintref constr_π[1:length(names1)]
    
    for i in 1:length(names1)
        constr_π[i] = @constraint(v_x, z[names1[i]] == x_hat[i])
    end
    
    add_bounds!(v_x, names2, y, BOUNDS, SEC_STG_COLS)

    # println(constr_λ, " ", typeof(constr_λ), " ", typeof(constr_λ[1]))
    # exit(0)
    v_x = Subproblem(v_x, constr_π, constr_λ, y)
    # JuMP.JuMPArray{JuMP.Variable,1,Tuple{Array{String,1}}}
    return v_x

end

function create_v_xs_integer(x_hat::Array{Float64,1}, SCENS::Array{Stoch_Scenario,1}, BOUNDS, FIRST_STG_COLS::Array{Array{String,1},1}, SEC_STG_COLS::Array{Array{String,1},1}, multi_thread::Bool)
    v_xs = Subproblem[]
    for i in 1:length(SCENS)
        prob = create_v_x_integer(x_hat, SCENS[i], BOUNDS, FIRST_STG_COLS, SEC_STG_COLS, multi_thread)
        push!(v_xs, prob)
    end
    v_xs
end


function create_v_x_integer(x_hat::Array{Float64,1}, SCENARIO::Stoch_Scenario, BOUNDS, FIRST_STG_COLS::Array{Array{String,1},1}, SEC_STG_COLS::Array{Array{String,1},1}, multi_thread::Bool)
    if multi_thread
        v_x = Model(solver=GurobiSolver(OutputFlag=0, InfUnbdInfo=1))
    else
        v_x = Model(solver=GurobiSolver(OutputFlag=0, InfUnbdInfo=1, Threads=1))
    end

    names1 = [var[1] for var in FIRST_STG_COLS]
    @variable(v_x, z[i = names1])
    names2 = [var[1] for var in SEC_STG_COLS]
    @variable(v_x, y[j = names2])

    # for (name, var_type) in SEC_STG_COLS
    #     if var_type == "int"
    #         setcategory(y[name], :Int)
    #     elseif var_type == "bin"
    #         setcategory(y[name], :Bin)
    #     end
    # end

    @objective(v_x, :Min, sum(SCENARIO.q[j] * y[names2[j]] for j in 1:length(names2)))

    # Tz + Wy ~ h
    # Its duals are λ: not necessary to get reference
    constr_λ_integer = JuMP.ConstraintRef[]
    # println("size es ", size(SCENARIO.W, 2))
    for r = 1:size(SCENARIO.T, 1)
        if SCENARIO.comps[r] == "E"
            push!(constr_λ_integer, @constraint(v_x, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) == SCENARIO.h[r]))
        elseif SCENARIO.comps[r] == "G"
            push!(constr_λ_integer, @constraint(v_x, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) >= SCENARIO.h[r]))
        else 
            push!(constr_λ_integer, @constraint(v_x, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) <= SCENARIO.h[r]))
        end
    end

    # z = x_hat
    # Its duals are π: get reference to return v_x, constr_π
    @constraintref constr_π_integer[1:length(names1)]
    
    for i in 1:length(names1)
        constr_π_integer[i] = @constraint(v_x, z[names1[i]] == x_hat[i])
    end
    
    add_bounds!(v_x, names2, y, BOUNDS, SEC_STG_COLS)

    v_x = Subproblem(v_x, constr_π_integer, constr_λ_integer, y)
    v_x
end




function update_subproblems!(v_xs::Array{Subproblem,1}, x_hat::Array{Float64,1})
    for k in 1:length(v_xs)
        update_subproblem!(v_xs[k], x_hat)
    end
end

function update_subproblem!(v_x::Subproblem, x_hat::Array{Float64,1})
    for i in 1:length(x_hat)
        JuMP.setRHS(v_x.constr_π[i], x_hat[i])
    end
end


function binarize_first_stg_cols(FIRST_STG_COLS::Array{Array{String,1},1}, BOUNDS::Array{Any,1})
    for i in 1:length(FIRST_STG_COLS)
        for bound in BOUNDS
            if bound[3] == FIRST_STG_COLS[i][1]
                if (bound[1] == "UP" || bound[1] == "UI") && (get(bound[4]) == 1 || get(bound[4]) == 1.0)
                    FIRST_STG_COLS[i] = [FIRST_STG_COLS[i][1], "bin"]
                end
            end
        end
    end
    return FIRST_STG_COLS
end


function insert_absent_elements!(SEC_STG_CONSTR, SEC_STG_ROWS, FIRST_STG_COLS::Array{Array{String,1},1}, SEC_STG_COLS::Array{Array{String,1},1})
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
    nothing
end