# This file implements a new way to input permutations in Julia. For example
# it is possible to create a permutation as follow
# pi = Oscar.Permutations.@perm (1,2,3)(4,5)(6,7,8)
# > (1,2,3)(4,5)(6,7,8)
# For this we use macros to modify the syntax tree of (1,2,3)(4,5)(6,7,8) such that
# Julia can deal with the expression.

module Permutations

 using GAP
 using Oscar

################################################################################
# Macro to input a permutation as
# julia> pi = Oscar.Permutations.@perm (1,2,3)(4,5)(6,7,8)
# > (1,2,3)(4,5)(6,7,8)
#
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

################################################################################
# Macro to input a list of permutation which are generated as elements of
# the symmetric_group(n)
# Example:
# julia> gens = Oscar.Permutations.@perm 14 [
#               (1,10)
#              (2,11)
#              (3,12)
#              (4,13)
#              (5,14)
#              (6,8)
#              (7,9)
#              (1,2,3,4,5,6,7)(8,9,10,11,12,13,14)
#              (1,2)(10,11)
#             ]
# > 9-element Vector{PermGroupElem}:
# (1,10)
# (2,11)
# (3,12)
# (4,13)
# (5,14)
# (6,8)
# (7,9)
# (1,2,3,4,5,6,7)(8,9,10,11,12,13,14)
# (1,2)(10,11)
#
# julia> gens[1].parent
# > Sym( [ 1 .. 14 ] )
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

################################################################################
# Generates a PermGroup with generators gens as a subgroup of symmetric_group(n)
# Example:
# julia> gens = Oscar.Permutations.@perm 14 [
#               (1,10)
#              (2,11)
#              (3,12)
#              (4,13)
#              (5,14)
#              (6,8)
#              (7,9)
#              (1,2,3,4,5,6,7)(8,9,10,11,12,13,14)
#              (1,2)(10,11)
#             ];
# julia> G = Oscar.Permutations.permgroup(14,gens)
# ><permutation group with 9 generators>
#
function permgroup(n::Int64,gens::Vector{PermGroupElem})

    return PermGroup(GAP.Globals.Subgroup(GAP.Globals.SymmetricGroup(GAP.Obj(n)),GAP.Obj([GAP.Obj(x) for x in gens ])))
end

################################################################################
# Computes the size of a PermGroup
# Example:
# julia> gens = Oscar.Permutations.@perm 14 [
#               (1,10)
#              (2,11)
#              (3,12)
#              (4,13)
#              (5,14)
#              (6,8)
#              (7,9)
#              (1,2,3,4,5,6,7)(8,9,10,11,12,13,14)
#              (1,2)(10,11)
#             ];
# julia> G = Oscar.Permutations.permgroup(14,gens);
# julia> Oscar.Permutations.size(G)
# > 645120
#
function size(G::PermGroup)
    return GAP.Globals.Size(G.X)
end

 export @perm

end #module Permutations

using .Permutations
