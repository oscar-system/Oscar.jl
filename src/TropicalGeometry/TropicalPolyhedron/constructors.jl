const MinOrMax = Union{typeof(min), typeof(max)}

struct TropicalPolyhedron{M<:MinOrMax} <: PolyhedralObject{QQFieldElem}
  pm_tpolytope::Polymake.BigObject

  TropicalPolyhedron{M}(p::Polymake.BigObject) where {M<:MinOrMax} = new{M}(p)
end

pm_object(P::TropicalPolyhedron) = P.pm_tpolytope

convention(P::TropicalPolyhedron{typeof(min)}) = min
convention(P::TropicalPolyhedron{typeof(max)}) = max

tropical_convex_hull(V::AbstractVector{T}) where {T<:AbstractVector{<:TropicalSemiringElem}} = tropical_convex_hull(matrix(V))

function tropical_convex_hull(V::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax}
  p = Polymake.tropical.Polytope{M.instance}(POINTS=V)

  return TropicalPolyhedron{M}(p)
end

function tropical_convex_hull(convention::MinOrMax, V::QQMatrix)
  p = Polymake.tropical.Polytope{convention}(POINTS=V)

  return TropicalPolyhedron{typeof(convention)}(p)
end

