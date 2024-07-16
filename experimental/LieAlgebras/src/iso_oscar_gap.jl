###############################################################################
#
#   Lie algebras
#
###############################################################################

function _iso_oscar_gap_lie_algebra_functions(
  LO::LieAlgebra{C}, LG::GapObj, coeffs_iso::MapFromFunc
) where {C<:FieldElem}
  basis_LG = GAPWrap.Basis(LG)

  f = function (x::LieAlgebraElem{C})
    cfs = GAP.Obj([coeffs_iso(c) for c in coefficients(x)])
    return GAPWrap.LinearCombination(basis_LG, cfs)
  end

  finv = function (x)
    cfs = Vector{GAP.Obj}(GAPWrap.Coefficients(basis_LG, x))
    return LO([preimage(coeffs_iso, c) for c in cfs])
  end

  return (f, finv)
end

function _iso_oscar_gap(LO::LinearLieAlgebra)
  coeffs_iso = Oscar.iso_oscar_gap(coefficient_ring(LO))
  LG = GAP.Globals.LieAlgebra(
    codomain(coeffs_iso),
    GAP.Obj([map_entries(coeffs_iso, xi) for xi in matrix_repr_basis(LO)]),
    GAP.Obj("basis"),
  )

  f, finv = _iso_oscar_gap_lie_algebra_functions(LO, LG, coeffs_iso)

  return MapFromFunc(LO, LG, f, finv)
end

function _iso_oscar_gap(LO::AbstractLieAlgebra; set_attributes::Bool=true)
  coeffs_iso = Oscar.iso_oscar_gap(coefficient_ring(LO))
  sc_table_G = [
    [
      [
        begin
          pairs = filter(pair -> !iszero(last(pair)), collect(enumerate(_matrix(xi * xj))))
          (map(first, pairs), GAP.Obj[coeffs_iso(c) for c in map(last, pairs)])
        end for xj in basis(LO)
      ] for xi in basis(LO)
    ]
    -1
    coeffs_iso(zero(coefficient_ring(LO)))
  ]

  LG = GAP.Globals.LieAlgebraByStructureConstants(
    codomain(coeffs_iso), GAP.Obj(sc_table_G; recursive=true)
  )

  f, finv = _iso_oscar_gap_lie_algebra_functions(LO, LG, coeffs_iso)

  if set_attributes && has_root_system(LO)
    # we need to construct the root system in GAP as otherwise it may detect a different order of simple roots
    RO = root_system(LO)
    RG = GAP.Globals.Objectify(
      GAP.Globals.NewType(
        GAP.Globals.NewFamily(GAP.Obj("RootSystemFam"), GAP.Globals.IsObject),
        GAP.evalstr("IsAttributeStoringRep and IsRootSystemFromLieAlgebra")),
      GAP.GapObj(Dict{Symbol,Any}()))
    GAP.Globals.SetUnderlyingLieAlgebra(RG, LG)

    cartan_trO = transpose(cartan_matrix(RO))
    transform_root(r::RootSpaceElem) = GAP.Obj(coefficients(r) * cartan_trO)[1]
    GAP.Globals.SetPositiveRoots(RG, GAP.Obj(transform_root.(positive_roots(RO))))
    GAP.Globals.SetNegativeRoots(RG, GAP.Obj(transform_root.(negative_roots(RO))))
    GAP.Globals.SetSimpleSystem(RG, GAP.Obj(transform_root.(simple_roots(RO))))
    can_basisG = GAP.Globals.CanonicalBasis(LG)
    pos_root_vectorsG = can_basisG[1:n_positive_roots(RO)]
    neg_root_vectorsG = can_basisG[(n_positive_roots(RO) + 1):(2 * n_positive_roots(RO))]
    csa_basisG = can_basisG[(2 * n_positive_roots(RO) + 1):end]
    GAP.Globals.SetPositiveRootVectors(RG, pos_root_vectorsG)
    GAP.Globals.SetNegativeRootVectors(RG, neg_root_vectorsG)
    GAP.Globals.SetCanonicalGenerators(
      RG, GAP.Obj([pos_root_vectorsG, neg_root_vectorsG, csa_basisG])
    )
    GAP.Globals.SetChevalleyBasis(
      LG, GAP.Obj([pos_root_vectorsG, neg_root_vectorsG, csa_basisG])
    )

    GAP.Globals.SetCartanMatrix(RG, GAP.Obj(cartan_trO))
    GAP.Globals.SetRootSystem(LG, RG)
  end

  return MapFromFunc(LO, LG, f, finv)
end
