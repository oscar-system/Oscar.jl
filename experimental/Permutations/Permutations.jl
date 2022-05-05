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


macro perm(n,gens)

    s = symmetric_group(n)
    ores = Vector{Expr}(undef,length(gens.args))
    i = 1
    for g in gens.args
        h = g.head
        arg = g.args
        res = []
        
        while h != :tuple
            arg1 = arg[1]
            insert!(res,1,Expr(:vect,arg[2:length(arg)]...))
            h = (arg1).head
            arg = (arg1).args
        end
        
        insert!(res,1,Expr(:vect,arg...))

        ores[i] = esc(:(Oscar.cperm(symmetric_group($n),$(res...))))
        i = i + 1
    end

    return Expr(:vect,ores...)
end


function permgroup(n::Int64,gens::Vector{PermGroupElem})

    return PermGroup(GAP.Globals.Subgroup(GAP.Globals.SymmetricGroup(GAP.Obj(n)),GAP.Obj([GAP.Obj(x) for x in gens ])))
end

function size(G::PermGroup)
    return GAP.Globals.Size(G.X)
end

       
end #module Permutations

using .Permutations
