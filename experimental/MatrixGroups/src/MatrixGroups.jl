module MatrixGroups

using GAP
using Oscar
import Oscar:GAPWrap

export _wrap_for_gap

# The corresponding gap functions are defined in the file gap/OscarInterface/gap/ExperimentalMatrixGroups.g

################################################################################
"""
   _wrap_for_gap(m::MatrixElem)

Compute the JuliaMatrixRep of `m` in GAP.

# Examples
```jldoctest
julia> m = matrix(ZZ, [0 1 ; -1 0]);

julia> Oscar._wrap_for_gap(m)
GAP: <matrix object of dimensions 2x2 over Integer ring>
```
"""
_wrap_for_gap(m::MatrixElem) = GAP.Globals.MakeJuliaMatrixRep(m)


################################################################################
"""
    matrix_group(matrices::Vector{<:MatrixElem{T}}; check::Bool = true) where T <: Union{ZZRingElem, QQFieldElem, AbsSimpleNumFieldElem}

Construct a GAP group `G` where the generators on the GAP side are wrappers
of type `JuliaMatrixRep` around the given Oscar matrices `matrices`.

If `G` is not finite then an exception is thrown.

A nice monomorphism from `G` to a GAP matrix group `G2` over a finite field
is stored in `G`, such that calculations in `G` can be handled automatically
by transferring them to `G2`.

# Examples
```jldoctest
julia> m1 = matrix(QQ, [0 1 ; -1 0]);

julia> m2 = matrix(QQ, [ -1 0; 0 1]);

julia> Oscar.MatrixGroups.matrix_group([m1, m2])
GAP: <group with 2 generators>
```
"""
function matrix_group(matrices::Vector{<:MatrixElem{T}}; check::Bool = true) where T <: Union{ZZRingElem, QQFieldElem, AbsSimpleNumFieldElem}
     # Compute the reduction map to a matrix group over a finite field `F`.
     flag, res = Oscar._isomorphic_group_over_finite_field(matrices, check = check)

     if !flag
       error("Group is not finite")
     end
     G, _, F, OtoFq = res

     # Map the generating matrices over `F` to GAP matrices, create a GAP group.
     matrices_Fq = [matrix(x) for x in gens(G)]  # Oscar matrices over F
     iso = Oscar.iso_oscar_gap(base_ring(matrices_Fq[1]))
     gap_matrices_Fq = [map_entries(iso, m) for m in matrices_Fq]
     G2 = GAP.Globals.Group(GapObj(gap_matrices_Fq))

     # Create a GAP group of wrapped matrices in characteristic zero.
     gapMatrices = [Oscar.MatrixGroups._wrap_for_gap(m) for m in matrices]
     G = GAP.Globals.Group(GapObj(gapMatrices))

     # Create a nice monomorphism from `G` to `G2`.
     # (`GroupHomomorphismByFunction` admits computing images via the
     # reduction `OtoFq`,
     # computing preimages is possible via the nice monomorphism of `G2`
     # which is an action homomorphism.)
     JuliaGAPMap = GAP.Globals.GroupHomomorphismByFunction(G, G2,
       M -> map_entries(iso, Oscar._reduce(GAP.getbangproperty(M, :m), OtoFq)))
     GAP.Globals.SetIsBijective(JuliaGAPMap, true)
     GAP.Globals.SetFilterObj(JuliaGAPMap, GAP.Globals.IsPreImagesByAction)
     GAP.Globals.SetNiceMonomorphism(G, JuliaGAPMap)
     GAP.Globals.SetIsHandledByNiceMonomorphism(G, true);

     return G
end

function _lex_isless(a::T,b::T) where T<:MatElem{S} where S <: Union{ZZRingElem, QQFieldElem}
  @assert base_ring(a) === base_ring(b)
  @assert size(a) == size(b)
  for i in 1:nrows(a), j in 1:ncols(a)
    if a[i,j] != b[i,j]
      return a[i,j] < b[i,j]
    end
  end
  return false
end

function _lex_isEqual(a::T,b::T) where T<:MatElem
  @assert base_ring(a) === base_ring(b)
  @assert size(a) == size(b)
  for i in 1:nrows(a), j in 1:ncols(a)
    if a[i,j] != b[i,j]
      return false
    end
  end
  return true
end

end #module MatrixGroups


using .MatrixGroups
