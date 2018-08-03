
function get_TIME(file_name)
    f = open(file_name) # En verdad no te van a dar full path: obtenerlo
    NAME = split(readline(f), " ")[end]

    PERIODS = []
    l = readline(f)
    l = readline(f)
    while !contains(l, "ENDATA")
        aux = split(l, " ")
        aux2 = []
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

function get_CORE(file_name)
    f = open(file_name)
    l = readline(f)
    while '*' == l[1]
        l = readline(f)
    end
    l = split(readline(f), " ")

    aux = []
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

function get_STOCH(file_name)
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


function get_first_stg_cols(CORE, STG2_C)
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

function get_first_stg_rows(CORE, STG2_R, object_name)
    FIRST_STG_ROWS = []
    for i in CORE[1]
        if i[2] == STG2_R
            break
        end
        if i[2] != object_name
            push!(FIRST_STG_ROWS, i)
        end
    end
    return FIRST_STG_ROWS
end

function get_first_stg_rhs(CORE, FIRST_STG_ROWS)
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
    return FIRST_STG_RHS
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

function get_sec_stg_rhs(CORE, SEC_STG_ROWS)
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
    return SEC_STG_RHS
end

function create_scenarios(numScens, STOCH, SEC_STG_RHS, SEC_STG_CONSTR, SEC_STG_ROWS, SEC_STG_COLS, FIRST_STG_COLS)
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

        
        h = [copy_RHS[name[2]] for name in SEC_STG_ROWS]
        q = [copy_AUX2[(name[1], STOCH[s][1][2])] for name in SEC_STG_COLS]
        T = [[copy_SEC_STG_CONSTR[name[2]][col[1]] for col in FIRST_STG_COLS] for name in SEC_STG_ROWS]
        T = hcat(T...)'
        W = [[copy_SEC_STG_CONSTR[name[2]][col[1]] for col in SEC_STG_COLS] for name in SEC_STG_ROWS]
        W = hcat(W...)'
        comps = [row[1] for row in SEC_STG_ROWS]
        s = Stoch_Scenario(get(STOCH[s][1][4]), h, q, T, W, comps, STOCH[s][1][2])
        push!(SCENS, s)
    end
    return SCENS
end



function linking_vars(SEC_STG_ROWS, FIRST_STG_COLS, SEC_STG_CONSTR)
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

function create_Lk(FIRST_STG_COLS, FIRST_STG_ROWS, FIRST_STG_RHS, SEC_STG_COLS, SCENARIO, BOUNDS)
    m = Model(solver=GurobiSolver(OutputFlag=0))
    names1 = [var[1] for var in FIRST_STG_COLS]
    # Faltan tipos (usar setcategory(x[i], :Int))
    @variable(m, z[i = names1])
    for (name, var_type) in FIRST_STG_COLS
        if var_type == "int"
            setcategory(z[name], :Int)
        elseif var_type == "bin"
            setcategory(z[name], :Bin)
        end
    end
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
    for (constr_type, name) in FIRST_STG_ROWS
        if constr_type == "E"
            @constraint(m, sum(val * z[var] for (var, val) in FIRST_STG_CONSTR[name]) == FIRST_STG_RHS[name])
        elseif constr_type == "G"
            @constraint(m, sum(val * z[var] for (var, val) in FIRST_STG_CONSTR[name]) >= FIRST_STG_RHS[name])
        else
            @constraint(m, sum(val * z[var] for (var, val) in FIRST_STG_CONSTR[name]) <= FIRST_STG_RHS[name])
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

    add_bounds!(m, names1, z, BOUNDS, FIRST_STG_COLS)
    add_bounds!(m, names2, y, BOUNDS, SEC_STG_COLS)
    

    return m
end

function init_master(FIRST_STG_COLS, FIRST_STG_OBJECT, FIRST_STG_ROWS, FIRST_STG_RHS, solver, BOUNDS)
    master = Model(solver=solver(OutputFlag=0))

    names1 = [var[1] for var in FIRST_STG_COLS]
    # Por ahora sin bounds
    @variable(master, x[i = names1])
    for (name, var_type) in FIRST_STG_COLS
        if var_type == "int"
            setcategory(x[name], :Int)
        elseif var_type == "bin"
            setcategory(x[name], :Bin)
        end
    end
    @variable(master, θ >= L)
    @objective(master, :Min, sum(FIRST_STG_OBJECT[name] * x[name] for name in names1) + θ)

    # agregar constraints Ax ~ b
    for (constr_type, name) in FIRST_STG_ROWS
        if constr_type == "E"
            @constraint(master, sum(val * x[var] for (var, val) in FIRST_STG_CONSTR[name]) == FIRST_STG_RHS[name])
        elseif constr_type == "G"
            @constraint(master, sum(val * x[var] for (var, val) in FIRST_STG_CONSTR[name]) >= FIRST_STG_RHS[name])
        else
            @constraint(master, sum(val * x[var] for (var, val) in FIRST_STG_CONSTR[name]) <= FIRST_STG_RHS[name])
        end
    end

    add_bounds!(master, names1, x, BOUNDS, FIRST_STG_COLS)


    return master, x, θ
end


function add_bounds!(m, names, z, BOUNDS, COLS)
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



function create_v_xs(x_hat, SCENS, numScens, BOUNDS, FIRST_STG_COLS, SEC_STG_COLS)
    v_xs = [create_v_x(x_hat, SCENS[i], BOUNDS, FIRST_STG_COLS, SEC_STG_COLS) for i in 1:numScens]
end

function create_v_x(x_hat, SCENARIO, BOUNDS, FIRST_STG_COLS, SEC_STG_COLS)
    v_x = Model(solver=GurobiSolver(OutputFlag=0))
    # Faltan bounds y tipos (usar set_var_type! ?)
    names1 = [var[1] for var in FIRST_STG_COLS]
    # Faltan tipos (usar set_var_type! ?)
    @variable(v_x, z[i = names1])
    names2 = [var[1] for var in SEC_STG_COLS]
    @variable(v_x, y[j = names2])

    @objective(v_x, :Min, sum(SCENARIO.q[j] * y[names2[j]] for j in 1:length(names2)))

    # Tz + Wy ~ h
    # Its duals are λ: not necessary to get reference
    for r = 1:size(SCENARIO.T, 1)
        if SCENARIO.comps[r] == "E"
            @constraint(v_x, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) == SCENARIO.h[r])
        elseif SCENARIO.comps[r] == "G"
            @constraint(v_x, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) >= SCENARIO.h[r])
        else 
            @constraint(v_x, sum(SCENARIO.T[r, c] * z[names1[c]] for c in 1:size(SCENARIO.T, 2)) + 
            sum(SCENARIO.W[r, c] * y[names2[c]] for c in 1:size(SCENARIO.W, 2)) <= SCENARIO.h[r])
        end
    end

    # z = x_hat
    # Its duals are π: get reference to return v_x, constrs_π
    @constraint(v_x, constr_π[i=1:length(names1)], z[names1[i]] == x_hat[i])
    add_bounds!(v_x, names2, y, BOUNDS, SEC_STG_COLS)

    return v_x, constr_π

end


function update_subproblems!(v_xs, x_hat)
    for k in 1:length(v_xs)
        update_subproblem!(v_xs[k][1], x_hat, v_xs[k][2])
    end
end

function update_subproblem!(v_k, x_hat, constrs_π)
    for i in 1:length(x_hat)
        JuMP.setRHS(constrs_π[i], x_hat[i])
    end
end


