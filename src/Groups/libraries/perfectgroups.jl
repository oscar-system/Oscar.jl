export
    number_perfect_groups,
    perfect_group,
    perfect_identification
    

###################################################################
# Perfect groups
###################################################################

"""
    perfect_group(n::Int, i::Int)
    perfect_group(::Type{T}, n::Int, i::Int)

Return the `i`-th group of order `n` and type `T` in the catalogue of perfect groups in the GAP Small Groups Library. The type `T` can be either ``PermGroup`` or ``FPGroup``. The group is given of type ``FPGroup`` if `T` is not specified.
"""
function perfect_group(::Type{T}, n::Int, m::Int) where T <: GAPGroup
   @assert m<= number_perfect_groups(n) "There are less than $m perfect groups of order $n."
   if T==PermGroup
      G = T(GAP.Globals.PerfectGroup(GAP.Globals.IsPermGroup,n,m))
   elseif T==FPGroup
      G = T(GAP.Globals.PerfectGroup(GAP.Globals.IsSubgroupFpGroup,n,m))
   else
      throw(ArgumentError("Wrong type"))
   end

   return G
end

perfect_group(n::Int, m::Int) = perfect_group(FPGroup,n,m)

"""
    perfect_identification(G::Group)

Return (`n`, `m`), where `G` = perfect_group(`n`,`m`).
"""
perfect_identification(G::GAPGroup) = GAP.gap_to_julia(GAP.Globals.PerfectIdentification(G.X))

"""
    number_perfect_groups(n::Int)

Return the number of perfect groups of order `n`.
"""
number_perfect_groups(n::Int) = GAP.Globals.NumberPerfectGroups(n)
