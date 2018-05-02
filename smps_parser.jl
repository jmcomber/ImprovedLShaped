

# START TIME

f = open("C:\\Jose\\Universidad\\WorkShopUAI\\semi.time")
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

# END TIME

# START CORE

f = open("C:\\Jose\\Universidad\\WorkShopUAI\\semi.core")


l = split(readline(f), " ")

aux = []
for elem in l
    if length(elem) > 0
        push!(aux, elem)
    end
end
NAME_CORE = aux[end]

if NAME_CORE != NAME
    # ERROR
    println("ERROR: Name mismatch in CORE, name is ", NAME_CORE, "\n\n")
end

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

# END CORE

# START STOCH(S)

f = open("C:\\Jose\\Universidad\\WorkShopUAI\\semi2.stoch")

l = readline(f)
while !contains(l, "STOCH")
    l = readline(f)
end
NAME_STOCH = split(l, " ")[end]

if NAME_STOCH != NAME
    # ERROR
    println("ERROR: Name mismatch in STOCH\n\n")
end

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
    while !contains(l, "SC") & !contains(l, "ENDATA") & length(l) > 4
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



for scen in SCENARIOS
    println(scen)
end


# END STOCH(S)


#  Header line
# – columns 1 through 14: first word field
# – columns 15 through 24: second word field
# Data line
# – columns 2 and 3: code field
# – columns 5–12: first name field
# – columns 15–22: second name field
# – columns 25–36: first numeric field
# – columns 40–47: third name field
# – columns 50–61: second numeric field

#Pero todo separado por espacios, y al parecer así se usa ahora