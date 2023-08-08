@doc raw"""
    AffineRationalPoint{CoeffType<:RingElem, ParentType<:AbsSpec}

A rational point represented in terms of a vector of coordinates.

Two rational points are considered equal if their parents have the same ambient
space and their coordinates agree.

# Examples
```jldoctest
julia> A2 = affine_space(GF(2), [:x, :y]);

julia> (x, y) = coordinates(A2);

julia> X = algebraic_set(x*y);

julia> X([1, 0])
Rational point
  of V(x*y)
  with coordinates (1, 0)

```
"""
struct AffineRationalPoint{S<:RingElem, T<:AbsSpec} <: AbsAffineRationalPoint{S, T}
  # not mutable since this should be efficient
  coordinates::Vector{S}
  parent::T

  function AffineRationalPoint(parent::T, coordinates::Vector{S}; check::Bool=true) where {S <: RingElem, T<:AbsSpec}
    r = new{S, T}(coordinates, parent)
    @check begin
      n = dim(ambient_space(parent))
      k = base_ring(parent)

      all(Oscar.parent(i)===k for i in coordinates) || error("coordinates do not lie in the base ring")
      n == length(coordinates) || error("$(coordinates) must have length $(n)")
      _in(r, parent)|| error("$(coordinates) does not lie in $(parent)")
    end
    return r
  end
end

@doc raw"""
    parent(p::AffineRationalPoint) -> AbsSpec

Return the parent ``X`` of the rational point ``p \in X``.
"""
parent(p::AffineRationalPoint) = p.parent

@doc raw"""
    coordinates(p::AffineRationalPoint{S,T}) -> Vector{S}

Return the coordinates of the rational point `p`.

The coordinates are with respect to the ambient space of its parent.
"""
coordinates(p::AffineRationalPoint) = p.coordinates

function (X::AbsSpec)(coordinates::Vector; check::Bool=true)
  k = base_ring(X)
  coordinates = k.(coordinates)
  return AffineRationalPoint(X, coordinates, check=check)
end

function Base.show(io::IO, ::MIME"text/plain", P::AbsAffineRationalPoint)
  io = pretty(io)
  println(io, "Rational point")
  print(io, Indent())
  println(io, "of ", parent(P))
  print(io, "with coordinates (")
  join(io, coordinates(P), ", ")
  print(io, ")", Dedent())
end

function Base.show(io::IO, P::AbsAffineRationalPoint)
  print(io,"(")
  join(io, coordinates(P), ", ")
  print(io,")")
end

function ==(P::AbsAffineRationalPoint{S,T}, Q::AbsAffineRationalPoint{S,U}) where {S,T,U}
  ambient_space(parent(P)) == ambient_space(parent(Q)) || return false
  return coordinates(P) == coordinates(Q)
end


function Base.hash(P::AbsAffineRationalPoint, h::UInt)
  return xor(hash(coordinates(P), h), hash(parent(P), h))
end

base_ring(P::AbsAffineRationalPoint) = base_ring(parent(P))
coefficient_ring(P::AbsAffineRationalPoint) = base_ring(parent(P))
coordinate(P::AbsAffineRationalPoint, i::Int) = coordinates(P)[i]
getindex(P::AbsAffineRationalPoint, i::Int) = coordinate(P, i)

@doc raw"""
    ideal(P::AbsAffineRationalPoint)

Return the maximal ideal associated to `P` in the  ambient coordinate ring
of its parent.
"""
function ideal(P::AbsAffineRationalPoint)
  R = ambient_coordinate_ring(parent(P))
  V = ambient_coordinates(parent(P))
  I = ideal(R, [V[i] - coordinate(P, i) for i in 1:length(V)])
  set_attribute!(I, :is_prime=>true)
  set_attribute!(I, :is_maximal=>true)
  set_attribute!(I, :is_absolutely_prime=>true)
  return I
end

function _in(P::AbsAffineRationalPoint, X::AbsSpec{<:Any,<:MPolyRing})
  # X is affine space
  return ambient_space(parent(P)) == X
end

function evaluate(f::MPolyRingElem, P::AbsAffineRationalPoint)
  return evaluate(f, coordinates(P))
end

function _in(P::AbsAffineRationalPoint, X::AbsSpec{<:Any,<:MPolyQuoRing})
  ambient_space(parent(P)) == ambient_space(X) || return false
  c = coordinates(P)
  for f in gens(ambient_closure_ideal(X))
    iszero(evaluate(f, c)) || return false
  end
  return true
end

function _in(P::AbsAffineRationalPoint, X::AbsAffineAlgebraicSet{<:Any,<:MPolyQuoRing})
  # we can do this without computing the vanishing ideal of X
  ambient_space(parent(P)) == ambient_space(X) || return false
  c = coordinates(P)
  for f in gens(fat_ideal(X))
    iszero(evaluate(f, c)) || return false
  end
  return true
end

function _in(P::AbsAffineRationalPoint, X::AbsSpec)
  # slow fallback for affine opens
  # this should be improved
  return issubset(scheme(P), X)
end

function Base.in(P::AbsAffineRationalPoint, X::AbsSpec)
  parent(P) === X && return true
  return _in(P, X)
end

function Base.in(P::AbsAffineRationalPoint, X::AbsCoveredScheme)
  return any(any(P in U for U in C) for C in coverings(X))
end

@doc raw"""
    scheme(P::AbsAffineRationalPoint) -> AbsSpec

Return the rational point ``P`` viewed as a reduced, affine subscheme
of its ambient affine space.
"""
function scheme(P::AbsAffineRationalPoint)
  I = ideal(P)
  R = ambient_coordinate_ring(parent(P))
  return Spec(R, I)
end

@doc raw"""
    closed_embedding(P::AbsAffineRationalPoint) -> ClosedEmbedding

Return the closed embedding of `P` into its parent scheme `X`.
"""
function closed_embedding(P::AbsAffineRationalPoint)
  I = ideal(P)
  X = parent(P)
  return ClosedEmbedding(X, I)
end

function (f::AbsSpecMor)(P::AbsAffineRationalPoint)
  @req domain(f) == parent(P) "$(P) not in domain"
  x = coordinates(codomain(f))
  g = pullback(f)
  p = coordinates(P)
  imgs = [evaluate(lift(g(y)),p) for y in x]
  return AffineRationalPoint(codomain(f), imgs, check=false)
end

function MPolyComplementOfKPointIdeal(P::AbsAffineRationalPoint)
  R = ambient_coordinate_ring(parent(P))
  return MPolyComplementOfKPointIdeal(R, coordinates(P))
end

function is_smooth_at(X::AbsSpec{<:Field}, P::AbsAffineRationalPoint)
  @req P in X "not a point on X"
  U = MPolyComplementOfKPointIdeal(P)
  R = OO(parent(P))
  XU = Spec(localization(R, U)[1])
  return is_smooth(XU)
end

@doc raw"""
    is_smooth(P::AbsAffineRationalPoint)

Return whether ``P`` is a smooth point of its parent ``X``.
"""
is_smooth(P::AbsAffineRationalPoint) = is_smooth_at(parent(P),P)


@doc raw"""
    tangent_space(X::AbsSpec{<:Field}, P::AbsAffineRationalPoint) -> AlgebraicSet

Return the Zariski tangent space of `X` at its rational point `P`.

See also [`tangent_space(P::AbsAffineRationalPoint{<:Field})`](@ref)
"""
function tangent_space(X::AbsSpec{<:Field}, P::AbsAffineRationalPoint)
  @req P in X "the point needs to lie on the algebraic set"
  J = jacobi_matrix(gens(ambient_closure_ideal(X)))
  v = coordinates(P)
  JP = map_entries(x->evaluate(x, v), J)
  V = ambient_coordinates(X)
  R = ambient_coordinate_ring(X)
  T = elem_type(R)[ sum([(V[i]-v[i])*JP[i,j] for i in 1:nrows(JP)], init=zero(R)) for j in 1:ncols(JP)]
  return algebraic_set(ideal(R,T), is_radical=true, check=false)
end

@doc raw"""
    tangent_space(P::AbsAffineRationalPoint{<:FieldElem}) -> AlgebraicSet

Return the Zariski tangent space of the parent of `P` at its point `P`.

See also [`tangent_space(X::AbsSpec{<:Field}, P::AbsAffineRationalPoint)`](@ref)
"""
tangent_space(P::AbsAffineRationalPoint{<:FieldElem}) = tangent_space(parent(P), P)

@doc raw"""
    is_du_val_singularity(P::AbsAffineRationalPoint{<:Field})

Return whether the parent of `P` has hat most a Du Val singularity at `P`.

Note that this includes the case that ``P`` is a smooth point.
"""
is_du_val_singularity(P::AbsAffineRationalPoint{<:FieldElem}) = is_du_val_singularity(parent(P), ideal(P))

@doc raw"""
    decide_du_val_singularity(P::AbsAffineRationalPoint{<:Field})

Return whether the parent of `P` has a Du Val singularity at `P`.

# Examples
```jldoctest
julia> A3 = affine_space(QQ, [:x, :y, :z]);

julia> (x, y, z) = ambient_coordinates(A3);

julia> X = subscheme(A3, ideal([x^2+y^2-z^2]));

julia> Oscar.decide_du_val_singularity(X([0,0,0]))
(true, (:A, 1))

```
"""
function decide_du_val_singularity(P::AbsAffineRationalPoint{<:FieldElem,<:Any})
  d = decide_du_val_singularity(parent(P), ideal(P))
  return d[1][1],d[1][3]
end
