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

  gapL = codomain(Oscar.iso_oscar_gap(L))
  dimL = GAPWrap.Dimension(gapL)
  @assert dimL == dim(L)
  basisL = GAPWrap.Basis(gapL)
  gapV = GAP.Globals.HighestWeightModule(gapL, GAP.Obj(weight))
  dimV = GAPWrap.Dimension(gapV)
  basisV = GAPWrap.Basis(gapV)

  struct_consts = Matrix{sparse_row_type(C)}(undef, dimL, dimV)
  for i in 1:dimL, j in 1:dimV
    struct_consts[i, j] = sparse_row(
      R,
      Tuple{Int,elem_type(R)}[
        (k, R(c)) for (k, c) in enumerate(
          map(
            c -> preimage(isoR, c),
            GAPWrap.Coefficients(GAPWrap.Basis(gapV), basisL[i]^basisV[j]),
          ),
        )
      ],
    )
  end

  return struct_consts
end
