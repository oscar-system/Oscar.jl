function lie_algebra_struct_consts_gap(R::Ring, dynkin::Tuple{Char,Int})
  @req is_valid_dynkin(dynkin...) "Input not allowed by GAP."

  isoR = Oscar.iso_oscar_gap(R)
  gapL = GAP.Globals.SimpleLieAlgebra(GAP.Obj(string(dynkin[1])), dynkin[2], codomain(isoR))
  dimL = GAPWrap.Dimension(gapL)
  comm_table_L =
    (
      entry -> (entry[1], Vector{elem_type(R)}(map(c -> preimage(isoR, c), entry[2])))
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

  return struct_consts, gapL
end

function gap_lie_algebra_by_struct_consts(L::LieAlgebra{C}) where {C<:RingElement}
  R = base_ring(L)
  isoR = Oscar.iso_oscar_gap(R)

  gap_sc_table = [
    [
      [
        begin
          pairs = filter(
            pair -> !iszero(last(pair)), collect(enumerate(Generic._matrix(xi * xj)))
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
