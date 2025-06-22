const MinOrMax = Union{typeof(min), typeof(max)}

struct TropicalPolyhedron{M<:MinOrMax} <: PolyhedralObject{QQFieldElem}
  pm_tpolytope::Polymake.BigObject

  TropicalPolyhedron{M}(p::Polymake.BigObject) where {M<:MinOrMax} = new{M}(p)
end
struct TropicalPointConfiguration{M<:MinOrMax} <: PolyhedralObject{QQFieldElem}
  pm_tpolytope::Polymake.BigObject

  TropicalPointConfiguration{M}(p::Polymake.BigObject) where {M<:MinOrMax} = new{M}(p)
end

pm_object(P::TropicalPolyhedron) = P.pm_tpolytope
pm_object(P::TropicalPointConfiguration) = P.pm_tpolytope

convention(P::TropicalPolyhedron{typeof(min)}) = min
convention(P::TropicalPolyhedron{typeof(max)}) = max
convention(P::TropicalPointConfiguration{typeof(min)}) = min
convention(P::TropicalPointConfiguration{typeof(max)}) = max

tropical_convex_hull(P::TropicalPointConfiguration{M}) where {M<:MinOrMax} = TropicalPolyhedron{M}(pm_object(P))
tropical_point_configuration(P::TropicalPolyhedron{M}) where {M<:MinOrMax} = TropicalPointConfiguration{M}(pm_object(P))

tropical_convex_hull(V::AbstractVector{T}) where {T<:AbstractVector{<:TropicalSemiringElem}} = tropical_convex_hull(matrix(V))
tropical_point_configuration(V::AbstractVector{T}) where {T<:AbstractVector{<:TropicalSemiringElem}} = tropical_point_configuration(matrix(V))

tropical_convex_hull(V::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax} = _tropical_convex_hull(TropicalPolyhedron, V)
tropical_point_configuration(V::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax} = _tropical_convex_hull(TropicalPointConfiguration, V)

tropical_convex_hull(convention::MinOrMax, V::QQMatrix) = _tropical_convex_hull(TropicalPolyhedron, convention, V)
tropical_point_configuration(convention::MinOrMax, V::QQMatrix) = _tropical_convex_hull(TropicalPointConfiguration, convention, V)

function _tropical_convex_hull(as::Union{Type{TropicalPolyhedron},Type{TropicalPointConfiguration}}, V::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax}
  p = Polymake.tropical.Polytope{M.instance}(POINTS=V)

  return as{M}(p)
end

function _tropical_convex_hull(as::Union{Type{TropicalPolyhedron},Type{TropicalPointConfiguration}}, convention::MinOrMax, V::QQMatrix)
  p = Polymake.tropical.Polytope{convention}(POINTS=V)

  return as{typeof(convention)}(p)
end

function tropical_polyhedron(A::MatElem{TropicalSemiringElem{M}}, B::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax}
  p = Polymake.tropical.Polytope{M.instance}(INEQUALITIES=(A,B))

  return TropicalPolyhedron{M}(p)
end
function tropical_polyhedron(convention::MinOrMax, A::QQMatrix, B::QQMatrix)
  p = Polymake.tropical.Polytope{convention}(INEQUALITIES=(A,B))

  return TropicalPolyhedron{typeof(convention)}(p)
end

