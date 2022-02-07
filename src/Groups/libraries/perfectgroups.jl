export
    number_perfect_groups,
    perfect_group,
    perfect_identification
    

###################################################################
# Perfect groups
###################################################################

"""
    perfect_group(::Type{T} = FPGroup, n::Int, i::Int)

Return the `i`-th group of order `n` and type `T` in the catalogue of
perfect groups in GAP's Perfect Groups Library.
The type `T` can be either `PermGroup` or `FPGroup`.
"""
function perfect_group(::Type{T}, n::Int, m::Int) where T <: GAPGroup
   N = number_perfect_groups(n)
   @assert m <= N "There are only $N perfect groups of order $n, up to isomorphism."
   if T==PermGroup
      G = T(GAP.Globals.PerfectGroup(GAP.Globals.IsPermGroup,n,m))
   elseif T==FPGroup
      G = T(GAP.Globals.PerfectGroup(GAP.Globals.IsSubgroupFpGroup,n,m))
   else
      throw(ArgumentError("Wrong type"))
   end

   return G
end

perfect_group(n::Int, m::Int) = perfect_group(PermGroup,n,m)

"""
    perfect_identification(G::GAPGroup)

Return `(n, m)` such that `G` is isomorphic with `perfect_group(n, m)`.
"""
perfect_identification(G::GAPGroup) = Tuple{Int,Int}(GAP.Globals.PerfectIdentification(G.X))

"""
    number_perfect_groups(n::Int)

Return the number of perfect groups of order `n`, up to isomorphism.
"""
number_perfect_groups(n::Int) = GAP.Globals.NumberPerfectGroups(n)
