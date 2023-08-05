@doc raw"""
    AffineRationalPoint{CoeffType<:RingElem, ParentType<:AbsSpec}

A $k$-rational point ``P`` of an affine scheme ``X``.

We refer to ``X`` as the parent of ``P``.

# Examples
```jldoctest
julia> A2 = affine_space(GF(2), [:x, :y]);
Rational point
  of Affine 2-space over GF(2) with coordinates [x, y]
  with coordinates (1, 0)

julia> A2([1, 0])

julia> (x, y) = coordinates(A2);

julia> X = algebraic_set(x*y);

julia> X([1,0])
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

Return the maximal ideal associated to `P` in the coordinate ring of its parent.
"""
function ideal(P::AbsAffineRationalPoint)
  R = ambient_coordinate_ring(parent(P))
  V = ambient_coordinates(parent(P))
  return ideal(R, [V[i] - coordinate(P, i) for i in 1:length(V)])
end

function _in(P::AbsAffineRationalPoint, X::AbsSpec{<:Any,<:MPolyRing})
  # X is affine space
  return ambient_space(parent(P)) == X
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
  return issubset(affine_scheme(P), X)
end

function Base.in(P::AbsAffineRationalPoint, X::AbsAffineAlgebraicSet)
  parent(P) === X && return true
  return _in(P, X)
end

@doc raw"""
    scheme(P::AbsAffineRationalPoint) -> AbsSpec

Return the rational point ``P`` viewed as a reduced, affine subscheme
of its ambient affine space.
"""
function scheme(P::AbsAffineRationalPoint)
  I = saturated_ideal(ideal(P))
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
