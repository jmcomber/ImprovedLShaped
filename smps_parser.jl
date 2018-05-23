

function get_TIME(file_name)
    # f = open("C:\\Jose\\Universidad\\WorkShopUAI\\semi.time")
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
    # f = open("C:\\Jose\\Universidad\\WorkShopUAI\\semi.core")
    f = open(file_name)
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
    # f = open("C:\\Jose\\Universidad\\WorkShopUAI\\semi2.stoch")
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
    l = readline(f)
    cont = 0
    while !contains(l, "ENDATA")
        l = split(l, " ")
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

    # for scen in SCENARIOS
    #     println(scen)
    # end
    close(f)

    return SCENARIOS

end


