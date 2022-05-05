module Permutations

 using GAP
 using Oscar
 using MacroTools

function HelloWorld()
     print("Hello World!")
 end
 
 macro perm(ex)
    MacroTools.postwalk(ex) do x
        @capture(x, f_(xs__)) || return x
        :($f, ($(xs...),))
    end
end


 macro perm(ex)
    MacroTools.postwalk(ex) do x
        @capture(x,  f_(args__) where f_ isa Expr ) || return x
        :($f, ($(xs...),))
    end
end
 
############
 
macro perm2(ex)
    MacroTools.postwalk(ex) do x
        @capture(x, f_(xs__)) || return x
        :($f, ($(xs...),))
    end
end
 
macro perm(ex)
    MacroTools.postwalk(ex) do x
        @capture(x, f_symbol(xs__) where isexpr(f_,String)) || return perm2(x)
        return x
    end
end

###############
# Without MacroTools
###############

# Working in some way
macro permNew(ex)

    h = (ex).head
    arg = (ex).args
    
    res = []

    while h != :tuple
        arg1 = arg[1]
        arguments = deleteat!(arg,1)
        insert!(res,1,arguments)
        h = (arg1).head
        arg = (arg1).args
    end
    
    insert!(res,1,arg)

    return res

end


macro permNew2(ex)

    h = (ex).head
    arg = (ex).args

    if h == :tuple
        return ex
    else
        @capture(ex, f_(xs__)) || print("hello world")
        res = [:($xs__)]
        
        return res
        
        h = f_.head

        while h != :tuple
            arg1 = arg[1]
            arguments = deleteat!(arg,1)
            insert!(res,1,arguments)
            h = (arg1).head
            arg = (arg1).args
        end

        insert!(res,1,arg)

        return res

    end
end


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


###############

macro perm2(ex)
    MacroTools.postwalk(ex) do x
        @capture(x, f_(xs__)) || return x
        :($f, ($(xs...),))
    end
end

macro perm3(ex)
    MacroTools.postwalk(ex) do x
        @capture(x, f_var(xs__)) || @capture(x, f_(xs__)) || return x
                                    :($f, ($(xs...),))
        :($f($(xs...)))
    end
end

macro test(ex)
    MacroTools.postwalk(ex) do x
        @capture(x, f_(f_(f_))) ||         @capture(x, f_(xs__)) || return x
                                         :($f, ($(xs...),))
        :($f($(xs...)))
    end
end

#function permutation_group(::Type{T}, n::Int) where T <: Oscar.GAPGroup
#  if n < 1
#    throw(ArgumentError("n must be a positive integer"))
#  end
#  return T(GAP.Globals.SymmetricGroup(_gap_filter(T), n))c
       
end #module Permutations

using .Permutations
