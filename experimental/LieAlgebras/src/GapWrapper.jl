function lie_algebra(
  gapL::GAP.Obj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(gapL)];
  cached::Bool=true,
)
  lie_algebra(AbstractLieAlgebra, gapL, s; cached)
end

function lie_algebra(
  ::Type{AbstractLieAlgebra},
  gapL::GAP.Obj,
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
  _set_gap_object!(L, gapL)
  return L
end

function lie_algebra(
  ::Type{LinearLieAlgebra},
  gapL::GAP.Obj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(gapL)];
  cached::Bool=true,
)
  @req GAP.Globals.IsLieObjectCollection(gapL) "Input is not a linear Lie algebra."

  gapR = GAPWrap.LeftActingDomain(gapL)
  isoR = Oscar.iso_gap_oscar(gapR)
  R = codomain(isoR)
  dimL = GAPWrap.Dimension(gapL)

  basis = [
    map_entries(isoR, GAP.Globals.UnderlyingRingElement(b)) for b in GAPWrap.Basis(gapL)
  ]
  n = size(basis[1])[1]

  L = LinearLieAlgebra{elem_type(R)}(R, n, basis, Symbol.(s); cached)
  _set_gap_object!(L, gapL)
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

function gap_lie_algebra_by_struct_consts(L::LieAlgebra{C}) where {C<:RingElement}
  R = base_ring(L)
  isoR = Oscar.iso_oscar_gap(R)

  gap_sc_table = [
    [
      [
        begin
          pairs = filter(
            pair -> !iszero(last(pair)), collect(enumerate(coefficients(xi * xj)))
          )
          (map(first, pairs), GAP.Obj[isoR(c) for c in map(last, pairs)])
        end for xj in basis(L)
      ] for xi in basis(L)
    ]
    -1
    isoR(zero(R))
  ]

  gapL = GAP.Globals.LieAlgebraByStructureConstants(
    codomain(isoR), GAP.Obj(gap_sc_table; recursive=true)
  )
  return gapL
end

function gap_lie_algebra_by_matrices(L::LinearLieAlgebra{C}) where {C<:RingElement}
  R = base_ring(L)
  isoR = Oscar.iso_oscar_gap(R)

  gapL = GAP.Globals.LieAlgebra(
    codomain(isoR),
    GAP.Obj([map_entries(isoR, xi) for xi in matrix_repr_basis(L)]),
    GAP.Obj("basis"),
  )

  return gapL
end

function lie_algebra_highest_weight_module_struct_consts_gap(
  L::LieAlgebra{C}, weight::Vector{Int}
) where {C<:RingElement}
  R = base_ring(L)
  isoR = Oscar.iso_oscar_gap(R)

  gapL = _gap_object(L)
  dimL = GAPWrap.Dimension(gapL)
  @assert dimL == dim(L)
  basisL = GAP.Globals.BasisVectors(GAPWrap.Basis(gapL))
  gapV = GAP.Globals.HighestWeightModule(gapL, GAP.Obj(weight))
  dimV = GAPWrap.Dimension(gapV)
  basisV = GAP.Globals.BasisVectors(GAPWrap.Basis(gapV))

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
