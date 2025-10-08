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

convention(::TropicalPolyhedron{typeof(min)}) = min
convention(::TropicalPolyhedron{typeof(max)}) = max
convention(::TropicalPointConfiguration{typeof(min)}) = min
convention(::TropicalPointConfiguration{typeof(max)}) = max

tropical_convex_hull(P::TropicalPointConfiguration{M}) where {M<:MinOrMax} = TropicalPolyhedron{M}(pm_object(P))
tropical_point_configuration(P::TropicalPolyhedron{M}) where {M<:MinOrMax} = TropicalPointConfiguration{M}(pm_object(P))

tropical_convex_hull(V::AbstractVector{T}) where {T<:AbstractVector{<:TropicalSemiringElem}} = tropical_convex_hull(matrix(V))
tropical_point_configuration(V::AbstractVector{T}) where {T<:AbstractVector{<:TropicalSemiringElem}} = tropical_point_configuration(matrix(V))

@doc raw"""
    tropical_convex_hull(V::MatElem{TropicalSemiringElem{M}}) where M<:MinOrMax

Construct the tropical convex hull of the vertices `V` in the tropical projective torus
with `M` as tropical addition.

# Examples
Defines a tropical unit ball in the 2-dimensional tropical projective torus with `max` as addition:
```jldoctest tropical_convex_hull
julia> T = tropical_semiring(max);

julia> A = T.(identity_matrix(QQ, 3));

julia> P = tropical_convex_hull(A)
Max tropical polyhedron in tropical projective torus of dimension 2
```

The vertices of a tropical convex hull may also be given by the rows of a tropical matrix:
```jldoctest tropical_convex_hull
julia> A = T[0 -4 -1; 0 -1 0];

julia> P = tropical_convex_hull(A)
Max tropical polyhedron in tropical projective torus of dimension 2
```
"""
tropical_convex_hull(V::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax} = _tropical_convex_hull(TropicalPolyhedron, V)
@doc raw"""
    tropical_point_configuration(V::MatElem{TropicalSemiringElem{M}}) where M<:MinOrMax

Construct the tropical point configuration of the vertices `V` in the tropical projective torus
with `M` as tropical addition.

# Examples
Defines a point configuration containing a tropical unit ball in the 2-dimensional tropical projective torus with `max` as addition:
```jldoctest tropical_point_configuration
julia> T = tropical_semiring(max);

julia> A = T.(identity_matrix(QQ, 3));

julia> P = tropical_point_configuration(A)
Max tropical point configuration in tropical projective torus of dimension 2
```
The tropical point configuration may also be given by the rows of a tropical matrix:
```jldoctest tropical_point_configuration
julia> A = T[0 -4 -1; 0 -1 0];

julia> P = tropical_point_configuration(A)
Max tropical point configuration in tropical projective torus of dimension 2
```
"""
tropical_point_configuration(V::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax} = _tropical_convex_hull(TropicalPointConfiguration, V)

@doc raw"""
    tropical_convex_hull(V::QQMatrix[,convention::MinOrMax=min])

Construct the tropical convex hull of the vertices `V` in the tropical projective torus
with `convention` as tropical addition.

# Examples
Define a point configuration containing a tropical unit ball in the 2-dimensional tropical projective torus with `max` as addition:
```jldoctest
julia> A = identity_matrix(QQ, 3);

julia> P = tropical_convex_hull(A,max)
Max tropical polyhedron in tropical projective torus of dimension 2
```
"""
tropical_convex_hull(V::QQMatrix, convention::MinOrMax=min) = _tropical_convex_hull(TropicalPolyhedron, convention, V)
@doc raw"""
    tropical_point_configuration(V::QQMatrix[, convention::MinOrMax=min])

Construct the tropical point configuration of the vertices `V` in the tropical projective torus
with `convention` as tropical addition.

# Examples
Define a point configuration containing a tropical unit ball in the 2-dimensional tropical projective torus with `max` as addition:
```jldoctest
julia> A = identity_matrix(QQ, 3);

julia> P = tropical_point_configuration(A,max)
Max tropical point configuration in tropical projective torus of dimension 2
```
"""
tropical_point_configuration(V::QQMatrix, convention::MinOrMax=min) = _tropical_convex_hull(TropicalPointConfiguration, convention, V)

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

Construct a `TropicalPolyhedron` as the set of points in the tropical projective torus satisfying the inequalities $A\odot x\leq B\odot x$
using `M` as the tropical addition.
"""
function tropical_polyhedron(A::MatElem{TropicalSemiringElem{M}}, B::MatElem{TropicalSemiringElem{M}}) where {M<:MinOrMax}
  p = Polymake.tropical.Polytope{M.instance}(INEQUALITIES=(A,B))

  return TropicalPolyhedron{M}(p)
end

function tropical_polyhedron(A::QQMatrix, B::QQMatrix, convention::MinOrMax=min)
  p = Polymake.tropical.Polytope{convention}(INEQUALITIES=(A,B))

  return TropicalPolyhedron{typeof(convention)}(p)
end

