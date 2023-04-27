###############################################################################
#
#   Lie algebras
#
###############################################################################

function _iso_gap_oscar(F::GAP.GapObj)
  if GAPWrap.IsLieAlgebra(F)
    if GAPWrap.IsLieObjectCollection(F)
      return _iso_gap_oscar_linear_lie_algebra(F)
    else
      return _iso_gap_oscar_abstract_lie_algebra(F)
    end
  end

  error("no method found")
end

function _iso_gap_oscar_abstract_lie_algebra(
  LG::GAP.GapObj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(LG)];
  cached::Bool=true,
)
  RG = GAPWrap.LeftActingDomain(LG)
  coeffs_iso = Oscar.iso_gap_oscar(RG)
  LO = _abstract_lie_algebra_from_GAP(LG, coeffs_iso, s; cached)
  finv, f = _iso_oscar_gap_lie_algebra_functions(LO, LG, inv(coeffs_iso))

  iso = MapFromFunc(f, finv, LG, LO)
  _set_iso_oscar_gap!(LO, inv(iso))
  return iso
end

function _iso_gap_oscar_linear_lie_algebra(
  LG::GAP.GapObj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(LG)];
  cached::Bool=true,
)
  RG = GAPWrap.LeftActingDomain(LG)
  coeffs_iso = Oscar.iso_gap_oscar(RG)
  LO = _linear_lie_algebra_from_GAP(LG, coeffs_iso, s; cached)
  finv, f = _iso_oscar_gap_lie_algebra_functions(LO, LG, inv(coeffs_iso))

  iso = MapFromFunc(f, finv, LG, LO)
  _set_iso_oscar_gap!(LO, inv(iso))
  return iso
end

function _abstract_lie_algebra_from_GAP(
  LG::GAP.GapObj, coeffs_iso::Map{GAP.GapObj}, s::Vector{<:VarName}; cached::Bool=true
)
  RO = codomain(coeffs_iso)
  dimL = GAPWrap.Dimension(LG)
  sc_table_G =
    (
      entry -> (entry[1], Vector{elem_type(RO)}(map(coeffs_iso, entry[2])))
    ).(
      Matrix{Tuple{Vector{Int},Vector{GAP.Obj}}}(
        (GAP.Globals.StructureConstantsTable(GAPWrap.Basis(LG)))[1:dimL]
      )
    )

  struct_consts = Matrix{SRow{elem_type(RO)}}(undef, dimL, dimL)
  for i in 1:dimL, j in 1:dimL
    struct_consts[i, j] = sparse_row(
      RO, Tuple{Int,elem_type(RO)}[(k, RO(c)) for (k, c) in zip(sc_table_G[i, j]...)]
    )
  end

  LO = AbstractLieAlgebra{elem_type(RO)}(RO, struct_consts, Symbol.(s); cached, check=false)
  return LO
end

function _linear_lie_algebra_from_GAP(
  LG::GAP.GapObj, coeffs_iso::Map{GAP.GapObj}, s::Vector{<:VarName}; cached::Bool=true
)
  @req GAPWrap.IsLieObjectCollection(LG) "Input is not a linear Lie algebra."

  RO = codomain(coeffs_iso)
  basis = [
    map_entries(coeffs_iso, GAPWrap.UnderlyingRingElement(b)) for b in GAPWrap.Basis(LG)
  ]
  n = size(basis[1])[1]
  LO = LinearLieAlgebra{elem_type(RO)}(RO, n, basis, Symbol.(s); cached)
  return LO
end

function lie_algebra(R::Ring, dynkin::Tuple{Char,Int}; cached::Bool=true)
  @req is_valid_dynkin(dynkin...) "Input not allowed by GAP."

  coeffs_iso = inv(Oscar.iso_oscar_gap(R))
  LG = GAP.Globals.SimpleLieAlgebra(
    GAP.Obj(string(dynkin[1])), dynkin[2], domain(coeffs_iso)
  )
  s = [Symbol("x_$i") for i in 1:GAPWrap.Dimension(LG)]
  LO = _abstract_lie_algebra_from_GAP(
    LG, coeffs_iso, s; cached
  )::AbstractLieAlgebra{elem_type(R)}

  finv, f = _iso_oscar_gap_lie_algebra_functions(LO, LG, inv(coeffs_iso))

  iso = MapFromFunc(f, finv, LG, LO)
  _set_iso_oscar_gap!(LO, inv(iso))
  return LO
end
