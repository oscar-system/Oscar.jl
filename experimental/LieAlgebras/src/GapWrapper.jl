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
    return GAP.Globals.LinearCombination(basis_LG, cfs)
  end

  finv = function (x)
    cfs = Vector{GAP.Obj}(GAPWrap.Coefficients(basis_LG, x))
    return LO([preimage(coeffs_iso, c) for c in cfs])
  end

  return (f, finv)
end

function _iso_oscar_gap(LO::LinearLieAlgebra{C}) where {C<:RingElement}
  _get_iso_oscar_gap!(LO) do
    coeffs_iso = Oscar.iso_oscar_gap(base_ring(LO))
    LG = GAP.Globals.LieAlgebra(
      codomain(coeffs_iso),
      GAP.Obj([map_entries(coeffs_iso, xi) for xi in matrix_repr_basis(LO)]),
      GAP.Obj("basis"),
    )

    f, finv = _iso_oscar_gap_lie_algebra_functions(LO, LG, coeffs_iso)

    MapFromFunc(f, finv, LO, LG)
  end
end

function _iso_oscar_gap(LO::AbstractLieAlgebra{C}) where {C<:RingElement}
  _get_iso_oscar_gap!(LO) do
    coeffs_iso = Oscar.iso_oscar_gap(base_ring(LO))
    sc_table_G = [
      [
        [
          begin
            pairs = filter(
              pair -> !iszero(last(pair)), collect(enumerate(Generic._matrix(xi * xj)))
            )
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

    MapFromFunc(f, finv, LO, LG)
  end
end

function _iso_gap_oscar(F::GAP.GapObj)
  if GAP.Globals.IsLieAlgebra(F)
    if GAP.Globals.IsLieObjectCollection(F)
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
  @req GAP.Globals.IsLieObjectCollection(LG) "Input is not a linear Lie algebra."

  RO = codomain(coeffs_iso)
  basis = [
    map_entries(coeffs_iso, GAP.Globals.UnderlyingRingElement(b)) for b in GAPWrap.Basis(LG)
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

################################################################################
#
#  Lie algebra modules
#
################################################################################

function lie_algebra_highest_weight_module_struct_consts_gap(
  L::LieAlgebra{C}, weight::Vector{Int}
) where {C<:RingElement}
  R = base_ring(L)
  isoR = Oscar.iso_oscar_gap(R)

  gapL = codomain(_iso_oscar_gap(L))
  dimL = GAPWrap.Dimension(gapL)
  @assert dimL == dim(L)
  basisL = GAPWrap.Basis(gapL)
  gapV = GAP.Globals.HighestWeightModule(gapL, GAP.Obj(weight))
  dimV = GAPWrap.Dimension(gapV)
  basisV = GAPWrap.Basis(gapV)

  struct_consts = Matrix{SRow{elem_type(R)}}(undef, dimL, dimV)
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
