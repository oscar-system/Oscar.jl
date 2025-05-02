###############################################################################
#
#   Lie algebras
#
###############################################################################

function _iso_gap_oscar_lie_algebra(F::GapObj)
  if GAPWrap.IsFiniteDimensional(F)
    if GAPWrap.IsLieObjectCollection(F)
      return _iso_gap_oscar_linear_lie_algebra(F)
    else
      return _iso_gap_oscar_abstract_lie_algebra(F)
    end
  end

  error("no method found")
end

push!(Oscar._iso_gap_oscar_methods, "IsLieAlgebra" => _iso_gap_oscar_lie_algebra)

function _iso_gap_oscar_abstract_lie_algebra(
  LG::GapObj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(LG)];
  coeffs_iso::Map{GapObj}=Oscar.iso_gap_oscar(GAPWrap.LeftActingDomain(LG)),
)
  LO = _abstract_lie_algebra_from_GAP(LG, coeffs_iso, s)
  finv, f = _iso_oscar_gap_lie_algebra_functions(LO, LG, inv(coeffs_iso))

  return MapFromFunc(LG, LO, f, finv)
end

function _iso_gap_oscar_linear_lie_algebra(
  LG::GapObj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(LG)];
  coeffs_iso::Map{GapObj}=Oscar.iso_gap_oscar(GAPWrap.LeftActingDomain(LG)),
)
  LO = _linear_lie_algebra_from_GAP(LG, coeffs_iso, s)
  finv, f = _iso_oscar_gap_lie_algebra_functions(LO, LG, inv(coeffs_iso))

  return MapFromFunc(LG, LO, f, finv)
end

function _abstract_lie_algebra_from_GAP(
  LG::GapObj, coeffs_iso::Map{GapObj}, s::Vector{<:VarName}
)
  RO = codomain(coeffs_iso)
  dimL = GAPWrap.Dimension(LG)
  sc_table_G =
    (
      entry -> (entry[1], Vector{elem_type(RO)}(map(coeffs_iso, entry[2])))
    ).(
      Matrix{Tuple{Vector{Int},Vector{GAP.Obj}}}(
        (GAPWrap.StructureConstantsTable(GAPWrap.Basis(LG)))[1:dimL]
      )
    )

  struct_consts = Matrix{sparse_row_type(RO)}(undef, dimL, dimL)
  for i in 1:dimL, j in 1:dimL
    struct_consts[i, j] = sparse_row(
      RO, Tuple{Int,elem_type(RO)}[(k, RO(c)) for (k, c) in zip(sc_table_G[i, j]...)]
    )
  end

  LO = AbstractLieAlgebra{elem_type(RO)}(RO, struct_consts, Symbol.(s); check=false)
  return LO
end

function _linear_lie_algebra_from_GAP(
  LG::GapObj, coeffs_iso::Map{GapObj}, s::Vector{<:VarName}
)
  @req GAPWrap.IsLieObjectCollection(LG) "Input is not a linear Lie algebra."

  RO = codomain(coeffs_iso)
  basis = [
    map_entries(coeffs_iso, GAPWrap.UnderlyingRingElement(b)) for b in GAPWrap.Basis(LG)
  ]
  n = size(basis[1])[1]
  LO = LinearLieAlgebra{elem_type(RO)}(RO, n, basis, Symbol.(s); check=false)
  return LO
end
