const MinOrMax = Union{typeof(min), typeof(max)}

struct TropicalPolyhedron{M<:MinOrMax} <: PolyhedralObject{QQFieldElem}
  pm_tpolytope::Polymake.BigObject

  @doc raw"""
      TropicalPolyhedron{M}(p::Polymake.BigObject) where M<:MinOrMax

  Constructs a `TropicalPolyhedron` corresponding to a `Polymake.BigObject` with tropical addition `M`.
  """
  TropicalPolyhedron{M}(p::Polymake.BigObject) where {M<:MinOrMax} = new{M}(p)
end
struct TropicalPointConfiguration{M<:MinOrMax} <: PolyhedralObject{QQFieldElem}
  pm_tpolytope::Polymake.BigObject

  @doc raw"""
      TropicalPointConfiguration{M}(p::Polymake.BigObject) where M<:MinOrMax

  Constructs a `TropicalPointConfiguration` corresponding to a `Polymake.BigObject` with tropical addition `M`.
  """
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

@doc raw"""
    tropical_convex_hull(V::MatElem{TropicalSemiringElem{M}}) where M<:MinOrMax

Constructs the tropical convex hull of the vertices `V` in the tropical projective torus
with `M` as tropical addition.
"""
tropical_convex_hull(V::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax} = _tropical_convex_hull(TropicalPolyhedron, V)
@doc raw"""
    tropical_point_configuration(V::MatElem{TropicalSemiringElem{M}}) where M<:MinOrMax

Constructs the tropical point configuration of the vertices `V` in the tropical projective torus
with `M` as tropical addition.
"""
tropical_point_configuration(V::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax} = _tropical_convex_hull(TropicalPointConfiguration, V)

@doc raw"""
    tropical_convex_hull(convention::MinOrMax, V::QQMatrix)

Constructs the tropical convex hull of the vertices `V` in the tropical projective torus
with `convention` as tropical addition.
"""
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

@doc raw"""
    tropical_polyhedron(A::MatElem{TropicalSemiringElem{M}}, B::MatElem{TropicalSemiringElem{M}}) where M<:MinOrMax

Constructs a `TropicalPolyhedron` as the set of points in the tropical projective torus satisfying the inequalities $A\odot x\leq B\odot x$
using `M` as the tropical addition.
"""
function tropical_polyhedron(A::MatElem{TropicalSemiringElem{M}}, B::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax}
  p = Polymake.tropical.Polytope{M.instance}(INEQUALITIES=(A,B))

  return TropicalPolyhedron{M}(p)
end
function tropical_polyhedron(convention::MinOrMax, A::QQMatrix, B::QQMatrix)
  p = Polymake.tropical.Polytope{convention}(INEQUALITIES=(A,B))

  return TropicalPolyhedron{typeof(convention)}(p)
end

struct ProjectiveTropicalPolyhedron{M<:MinOrMax} <: PolyhedralObject{QQFieldElem}
  
end

