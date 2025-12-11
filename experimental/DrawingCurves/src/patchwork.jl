function _patchworks2polynomial(
  hs::Polymake.BigObject; t=QQ(1, 2), points::Union{AbstractMatrix,String,Nothing}=nothing
)
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  return _patchworks2polynomial(hs, R; t, points)
end
function _patchworks2polynomial(
  hs::Polymake.BigObject,
  R::MPolyRing;
  t=QQ(1, 2),
  points::Union{AbstractMatrix,String,Nothing}=nothing,
)
  pws = Polymake._lookup_multi(hs, "PATCHWORK")

  if hs.MONOMIALS !== nothing
    @req ncols(hs.MONOMIALS) == 3 "Wrong ambient dimension of hypersurface"
    pts = matrix(ZZ, hs.MONOMIALS)
    pts = pts[:, 2:end]
  else
    @req points !== nothing "Unable to determine points"
    if points isa AbstractMatrix
      pts = matrix(ZZ, points)
    elseif points isa String
      pts = matrix(ZZ, Polymake.load(points))
      pts = pts[:, 2:end]
    end
    @req ncols(pts) == 2 "Wrong embedding dimension for points"
  end

  delta = convert(Int, maximum([sum(pts[i, :]) for i in 1:nrows(pts)]))

  ds = subdivision_of_points(Polymake.tropical.dual_subdivision(hs))
  mc = maximal_cells(IncidenceMatrix, ds)
  # Since we still get faulty data...
  mc = incidence_matrix([collect(row(mc, i)) for i in 1:nrows(mc)])
  sop = subdivision_of_points(pts, mc)
  # @assert is_regular(sop) "Subdivision is not regular"
  w = min_weights(sop)

  result = elem_type(R)[]
  for pw in pws
    signs = convert(Vector{Int}, Polymake.@convert_to Vector{Int} pw.SIGNS)
    push!(result, _patchwork2polynomial(R, signs, pts, w, delta, t))
  end
  return result
end

function _patchwork2polynomial(
  R::QQMPolyRing,
  sgns::AbstractVector,
  exps::ZZMatrix,
  w::Vector{Int},
  delta::Int,
  t::QQFieldElem,
)
  @req length(sgns) == length(w) "Dimensions do not agree"
  @req length(sgns) == nrows(exps) "Dimensions do not agree"

  ctx = MPolyBuildCtx(R)
  for i in 1:nrows(exps)
    coeff = t^w[i]
    if sgns[i] == 1
      coeff = -coeff
    end
    pt = exps[i, :]
    zc = delta - sum(pt)
    exp = convert(Vector{Int}, vcat(pt, [zc]))
    push_term!(ctx, coeff, exp)
  end
  return finish(ctx)
end
