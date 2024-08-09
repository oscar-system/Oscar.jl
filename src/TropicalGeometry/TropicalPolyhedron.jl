const MinOrMax = Union{typeof(min), typeof(max)}

struct TropicalPolyhedron{M<:MinOrMax} <: PolyhedralObject{QQFieldElem}
  pm_tpolytope::Polymake.BigObject

  TropicalPolyhedron{M}(p::Polymake.BigObject) where {M<:MinOrMax} = new{M}(p)
end

pm_object(P::TropicalPolyhedron) = P.pm_tpolytope

convention(P::TropicalPolyhedron{typeof(min)}) = min
convention(P::TropicalPolyhedron{typeof(max)}) = max

function tropical_convex_hull(V::AbstractVector{T}) where {T<:AbstractVector{<:TropicalSemiringElem}}
  M = convention(V)
  
  p = Polymake.tropical.Polytope{M}(POINTS=matrix(V))
  return TropicalPolyhedron{typeof{M}}(p)
end

function tropical_convex_hull(V::MatElem{<:TropicalSemiringElem})
  M = convention(V)
  p = Polymake.tropical.Polytope{M}(POINTS=V)

  return TropicalPolyhedron{typeof(M)}(p)
end

function tropical_convex_hull(convention::MinOrMax, V::QQMatrix)
  p = Polymake.tropical.Polytope{convention}(POINTS=V)

  return TropicalPolyhedron{typeof(convention)}(p)
end

function vertices(as::Type{PointVector{T}}, P::TropicalPolyhedron) where {T<:TropicalSemiringElem}
  n = length(P.pm_tpolytope.VERTICES_IN_POINTS)

  return SubObjectIterator{as}(P, _vertex_of_polyhedron, n)
end

vertices(P::TropicalPolyhedron{M}) where {M<:MinOrMax} = vertices(PointVector{TropicalSemiringElem{M}}, P)

function _vertex_of_polyhedron(::Type{PointVector{TropicalSemiringElem{M}}}, P::TropicalPolyhedron, i::Int) where {M<:MinOrMax}
  T = tropical_semiring(convention(P))

  return point_vector(
    T,
    T.(P.pm_tpolytope.POINTS[i,:])
  )
end

n_vertices(P::TropicalPolyhedron) = size(P.VERTICES, 1)-1

function Base.show(io::IO, P::TropicalPolyhedron{T}) where {T<:MinOrMax}
  if T == typeof(min)
    print(io, "Min tropical polyhedron")
  else
    print(io, "Max tropical polyhedron")
  end
end

