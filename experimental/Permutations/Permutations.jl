module Permutations

 using GAP
 using Oscar

macro perm(ex)
    h = (ex).head
    arg = (ex).args
    res = []
    
    while h != :tuple
        arg1 = arg[1]
        insert!(res,1,Expr(:vect,arg[2:length(arg)]...))
        h = (arg1).head
        arg = (arg1).args
    end
    
    insert!(res,1,Expr(:vect,arg...))

    return esc(:(Oscar.cperm($(res...))))
end


macro permGens(n, gens)

    return [permNew(g) for g in gens]

end

       
end #module Permutations

using .Permutations
