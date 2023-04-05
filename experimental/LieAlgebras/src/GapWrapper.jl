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
    sc_table = [
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
      codomain(coeffs_iso), GAP.Obj(sc_table; recursive=true)
    )

    f, finv = _iso_oscar_gap_lie_algebra_functions(LO, LG, coeffs_iso)

    MapFromFunc(f, finv, LO, LG)
  end
end

function lie_algebra(
  gapL::GAP.GapObj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(gapL)];
  cached::Bool=true,
)
  lie_algebra(AbstractLieAlgebra, gapL, s; cached)
end

function lie_algebra(
  ::Type{AbstractLieAlgebra},
  gapL::GAP.GapObj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(gapL)];
  cached::Bool=true,
)
  gapR = GAPWrap.LeftActingDomain(gapL)
  isoR = Oscar.iso_gap_oscar(gapR)
  R = codomain(isoR)
  dimL = GAPWrap.Dimension(gapL)
  comm_table_L =
    (
      entry -> (entry[1], Vector{elem_type(R)}(map(isoR, entry[2])))
    ).(
      Matrix{Tuple{Vector{Int},Vector{GAP.Obj}}}(
        (GAP.Globals.StructureConstantsTable(GAPWrap.Basis(gapL)))[1:dimL]
      )
    )

  struct_consts = Matrix{SRow{elem_type(R)}}(undef, dimL, dimL)
  for i in 1:dimL, j in 1:dimL
    struct_consts[i, j] = sparse_row(
      R, Tuple{Int,elem_type(R)}[(k, R(c)) for (k, c) in zip(comm_table_L[i, j]...)]
    )
  end

  L = AbstractLieAlgebra{elem_type(R)}(R, struct_consts, Symbol.(s); cached, check=false)
  # _set_gap_object!(L, gapL)
  return L
end

function lie_algebra(
  ::Type{LinearLieAlgebra},
  gapL::GAP.GapObj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(gapL)];
  cached::Bool=true,
)
  @req GAP.Globals.IsLieObjectCollection(gapL) "Input is not a linear Lie algebra."

  gapR = GAPWrap.LeftActingDomain(gapL)
  isoR = Oscar.iso_gap_oscar(gapR)
  R = codomain(isoR)

  basis = [
    map_entries(isoR, GAP.Globals.UnderlyingRingElement(b)) for b in GAPWrap.Basis(gapL)
  ]
  n = size(basis[1])[1]

  L = LinearLieAlgebra{elem_type(R)}(R, n, basis, Symbol.(s); cached)
  # _set_gap_object!(L, gapL)
  return L
end

function lie_algebra(R::Ring, dynkin::Tuple{Char,Int}; cached::Bool=true)
  @req is_valid_dynkin(dynkin...) "Input not allowed by GAP."

  isoR = Oscar.iso_oscar_gap(R)
  gapL = GAP.Globals.SimpleLieAlgebra(GAP.Obj(string(dynkin[1])), dynkin[2], codomain(isoR))

  L = lie_algebra(AbstractLieAlgebra, gapL; cached)::AbstractLieAlgebra{elem_type(R)}
  set_attribute!(L, :dynkin => dynkin)
  return L
end

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
