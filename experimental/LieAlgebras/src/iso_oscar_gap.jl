###############################################################################
#
#   Lie algebras
#
###############################################################################

function _iso_oscar_gap_lie_algebra_functions(
  LO::LieAlgebra{C}, LG::GAP.GapObj, coeffs_iso::MapFromFunc
) where {C<:RingElement}
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
  coeffs_iso = Oscar.iso_oscar_gap(base_ring(LO))
  LG = GAP.Globals.LieAlgebra(
    codomain(coeffs_iso),
    GAP.Obj([map_entries(coeffs_iso, xi) for xi in matrix_repr_basis(LO)]),
    GAP.Obj("basis"),
  )

  f, finv = _iso_oscar_gap_lie_algebra_functions(LO, LG, coeffs_iso)

  return MapFromFunc(LO, LG, f, finv)
end

function _iso_oscar_gap(LO::AbstractLieAlgebra)
  coeffs_iso = Oscar.iso_oscar_gap(base_ring(LO))
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
    coeffs_iso(zero(base_ring(LO)))
  ]

  LG = GAP.Globals.LieAlgebraByStructureConstants(
    codomain(coeffs_iso), GAP.Obj(sc_table_G; recursive=true)
  )

  f, finv = _iso_oscar_gap_lie_algebra_functions(LO, LG, coeffs_iso)

  return MapFromFunc(LO, LG, f, finv)
end
