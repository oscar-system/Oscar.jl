struct MonomialBasis
  ZZx::ZZMPolyRing
  set_mon::Set{ZZMPolyRingElem}
  dimension::Int
  no_minkowski::Set{Vector{Int}}
  polytope::Oscar.Polymake.BigObjectAllocated
end

function MonomialBasis(
  ZZx::ZZMPolyRing, set_mon::Set{ZZMPolyRingElem}, no_minkowski::Set{Vector{ZZRingElem}}
)
  vertices = degrees.(collect(set_mon))
  vertices_hom = transpose(reduce(hcat, [prepend!(vec, 1) for vec in vertices])) # homogenoues coordinate system
  poly = Oscar.Polymake.polytope.Polytope(; POINTS=vertices_hom)
  return MonomialBasis(ZZx, set_mon, length(set_mon), no_minkowski, poly)
end

function Base.show(io::IO, monomial_basis::MonomialBasis)
  println(io, "MonomialBasis")
  println(io, "Dimension: ", monomial_basis.dimension)
  println(io, "Generators within semi-group: ", monomial_basis.no_minkowski)
  println(
    io,
    "First 10 Monomials in degrevlex: ",
    sort(
      collect(monomial_basis.set_mon);
      lt=get_monomial_order_lt("degrevlex", monomial_basis.ZZx),
    )[1:min(end, 10)],
  )
  # print(io, "Volume polytope: ", monomial_basis.polytope.VOLUME)
end
