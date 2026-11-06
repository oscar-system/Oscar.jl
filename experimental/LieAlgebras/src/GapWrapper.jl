################################################################################
#
#  Lie algebra modules
#
################################################################################

function lie_algebra_simple_module_struct_consts_gap(
  L::LieAlgebra, weight::WeightLatticeElem
)
  @req root_system(weight) == root_system(L) "Incompatible root systems"
  return lie_algebra_simple_module_struct_consts_gap(
    L, Int.(Oscar._vec(coefficients(weight)))
  )
end

function lie_algebra_simple_module_struct_consts_gap(
  L::LieAlgebra{C}, weight::Vector{Int}
) where {C<:FieldElem}
  R = coefficient_ring(L)
  isoR = Oscar.iso_oscar_gap(R)

  gapL = _gap_lie_algebra_with_root_system(L)
  dimL = GAPWrap.Dimension(gapL)
  @assert dimL == dim(L)
  basisL = GAPWrap.Basis(gapL)
  gapV = GAP.Globals.HighestWeightModule(gapL, GAP.Obj(weight))
  dimV = GAPWrap.Dimension(gapV)
  basisV = GAPWrap.Basis(gapV)

  struct_consts = Matrix{sparse_row_type(C)}(undef, dimL, dimV)
  for i in 1:dimL, j in 1:dimV
    entries = Tuple{Int,elem_type(R)}[]
    for (k, c) in enumerate(GAPWrap.Coefficients(basisV, basisL[i]^basisV[j]))
      # only convert the non-zeros, `iszero` is a GAP call unless `c` is a Julia integer
      !iszero(c) && push!(entries, (k, preimage(isoR, c)))
    end
    struct_consts[i, j] = sparse_row(R, entries)
  end

  return struct_consts
end

# The GAP counterpart of `L` with the simple roots of `root_system(L)`, so that GAP reads
# weights in the same order. GAP chooses its own if the cached isomorphism was created
# before the root system of `L`.
function _gap_lie_algebra_with_root_system(L::LieAlgebra)
  root_system(L) # `iso_oscar_gap` only passes a root system to GAP if `L` already has one
  isoL = Oscar.iso_oscar_gap(L)
  gapL = codomain(isoL)
  if !GAP.Globals.HasRootSystem(gapL)
    _iso_oscar_gap_set_root_system(L, gapL, isoL)
    return gapL
  end
  # GAP may have computed a root system of its own
  chev = GAP.Obj([GAP.Obj(isoL.(xs)) for xs in chevalley_basis(L)])
  GAP.Globals.ChevalleyBasis(gapL) == chev && return gapL
  # a fresh GAP copy of `L`, created with the root system of `L`
  return codomain(_iso_oscar_gap(L))
end
