module Permutations

 using GAP
 using Oscar
 using MacroTools

macro permNew3(ex)
    h = (ex).head
    arg = (ex).args
    res = []
    
    while h != :tuple
        arg1 = arg[1]
        arguments = deleteat!(arg,1)
        for i in 1:length(arguments)
            if !isexpr(arguments[i],Int64)
                arguments[i] = eval(arguments[i])
            end
        end
        insert!(res,1,Tuple(x for x in arguments))
        h = (arg1).head
        arg = (arg1).args
    end
    
    for i in 1:length(arg)
        if !isexpr(arg[i],Int64)
            arg[i] = eval(arg[i])
        end
    end
    insert!(res,1,Tuple(x for x in arg))
    
    n = length(res)
    if n == 1
        return res[1]
    end
    if n == 2
        return (res[1],res[2])
    end

    r = (res[1],res[2])
    for i in 3:n
        r = (r, res[i])
    end

    return r
end


macro permGens(n, gens)

    return [permNew(g) for g in gens]

end

       
end #module Permutations

using .Permutations
