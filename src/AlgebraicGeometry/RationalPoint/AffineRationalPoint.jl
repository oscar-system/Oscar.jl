function (X::AbsAffineScheme{BRT})(p::AbsProjectiveRationalPoint) where BRT
  # certainly not the fastest for the standard affine covering
  # ... but it will work in general.
  P = codomain(p)
  h = homogenization_map(P, X)
  q = elem_type(BRT)[]
  for g in gens(OO(X))
    n,d = lift.(h(g))
    c = evaluate(n,coordinates(p))*inv(evaluate(d,coordinates(p)))
    push!(q, c)
  end
  return X(q)
end

# implement the AbsRationalPoint interface
parent(p::AffineRationalPoint) = p.parent

coordinate(P::AbsAffineRationalPoint, i::Int) = coordinates(P)[i]
Base.getindex(P::AbsAffineRationalPoint, i::Int) = coordinate(P, i)

@doc raw"""
    coordinates(p::AffineRationalPoint{S,T}) -> Vector{S}

Return the coordinates of the rational point `p`.

The coordinates are with respect to the ambient space of its ambient scheme.
"""
coordinates(p::AffineRationalPoint) = p.coordinates

function (X::AbsAffineScheme)(coordinates::Vector; check::Bool=true)
  k = base_ring(X)
  coordinates = k.(coordinates)
  return AffineRationalPoint(rational_point_set(X), coordinates, check=check)
end

function (X::AbsAffineScheme)(p::AffineRationalPoint; check::Bool=true)
  if codomain(p) === X
    return X
  end
  return X(coordinates(p);check=check)
end

function Base.show(io::IO, ::MIME"text/plain", P::AbsAffineRationalPoint)
  io = pretty(io)
  println(io, "Rational point")
  print(io, Indent())
  println(io, "of ", Lowercase(), ambient_scheme(P), Dedent())
  print(io, "with coordinates (")
  join(io, coordinates(P), ", ")
  print(io, ")")
end

function Base.show(io::IO, P::AbsAffineRationalPoint)
  print(io,"(")
  join(io, coordinates(P), ", ")
  print(io,")")
end

function ==(P::AbsAffineRationalPoint{S,T}, Q::AbsAffineRationalPoint{S,U}) where {S,T,U}
  ambient_space(P) == ambient_space(Q) || return false
  return coordinates(P) == coordinates(Q)
end


function Base.hash(P::AbsAffineRationalPoint, h::UInt)
  return xor(hash(coordinates(P), h), hash(ambient_scheme(P), h))
end


@doc raw"""
    ideal(P::AbsAffineRationalPoint)

Return the maximal ideal associated to `P` in the coordinate ring
of its ambient space.
"""
function ideal(P::AbsAffineRationalPoint)
  R = OO(ambient_space(P))
  V = coordinates(ambient_space(P))
  I = ideal(R, [V[i] - coordinate(P, i) for i in 1:length(V)])
  set_attribute!(I, :is_prime=>true)
  set_attribute!(I, :is_maximal=>true)
  set_attribute!(I, :is_absolutely_prime=>true)
  return I
end

function _in(P::AbsAffineRationalPoint, X::AbsAffineScheme{<:Any,<:MPolyRing})
  # X is affine space
  return ambient_space(P) == X
end

function evaluate(f::MPolyRingElem, P::AbsAffineRationalPoint)
  return evaluate(f, coordinates(P))
end

function _in(P::AbsAffineRationalPoint, X::AbsAffineScheme{<:Any,<:MPolyQuoRing})
  ambient_space(P) == ambient_space(X) || return false
  c = coordinates(P)
  for f in gens(saturated_ideal(defining_ideal(X)))
    iszero(evaluate(f, c)) || return false
  end
  return true
end

function _in(P::AbsAffineRationalPoint, X::AbsAffineAlgebraicSet{<:Any,<:MPolyQuoRing})
  # we can do this without computing the vanishing ideal of X
  ambient_space(P) == ambient_space(X) || return false
  c = coordinates(P)
  for f in gens(fat_ideal(X))
    iszero(evaluate(f, c)) || return false
  end
  return true
end

function _in(P::AbsAffineRationalPoint, X::AbsAffineScheme)
  # slow fallback for affine opens
  # this should be improved
  return is_subscheme(scheme(P), X)
end

function Base.in(P::AbsAffineRationalPoint, X::AbsAffineScheme)
  codomain(P) === X && return true
  return _in(P, X)
end

function Base.in(P::AbsAffineRationalPoint, X::AbsCoveredScheme)
  return any(any(P in U for U in C) for C in coverings(X))
end

@doc raw"""
    scheme(P::AbsAffineRationalPoint) -> AbsAffineScheme

Return the rational point ``P`` viewed as a reduced, affine subscheme
of its ambient affine space.
"""
function scheme(P::AbsAffineRationalPoint)
  I = ideal(P)
  R = OO(ambient_space(P))
  return spec(R, I)
end

@doc raw"""
    closed_embedding(P::AbsAffineRationalPoint) -> ClosedEmbedding

Return the closed embedding of `P` into its ambient scheme `X`.
"""
function closed_embedding(P::AbsAffineRationalPoint)
  I = ideal(P)
  X = ambient_scheme(P)
  return ClosedEmbedding(X, ideal(OO(X),OO(X).(gens(I))))
end


function (f::AbsAffineSchemeMor{<:AbsAffineScheme{S},<:AbsAffineScheme{S},<:Any,<:Any,Nothing})(P::AbsAffineRationalPoint) where {S}
  # The Nothing type parameter assures that the base morphism is trivial.
  @req domain(f) == codomain(P) "$(P) not in domain"
  @req base_ring(domain(f)) == base_ring(codomain(f)) "schemes must be defined over the same base ring. Try to map the point as an ideal instead"
  x = coordinates(codomain(f))
  g = pullback(f)
  p = coordinates(P)
  imgs = [evaluate(lift(g(y)),p) for y in x]
  return codomain(f)(imgs; check=false)
end

function MPolyComplementOfKPointIdeal(P::AbsAffineRationalPoint)
  R = OO(ambient_space(P))
  return MPolyComplementOfKPointIdeal(R, coordinates(P))
end

function is_smooth(X::AbsAffineScheme{<:Field}, P::AbsAffineRationalPoint)
  @req P in X "not a point on X"
  U = MPolyComplementOfKPointIdeal(P)
  R = OO(codomain(P))
  XU = spec(localization(R, U)[1])
  return is_smooth(XU)
end

@doc raw"""
    is_smooth(P::AbsAffineRationalPoint)

Return whether ``P`` is a smooth point of its ambient scheme ``X``.
"""
is_smooth(P::AbsAffineRationalPoint) = is_smooth(codomain(P),P)


@doc raw"""
    tangent_space(X::AbsAffineScheme{<:Field}, P::AbsAffineRationalPoint) -> AlgebraicSet

Return the Zariski tangent space of `X` at its rational point `P`.

See also [`tangent_space(P::AbsAffineRationalPoint{<:Field})`](@ref)
"""
function tangent_space(X::AbsAffineScheme{<:Field}, P::AbsAffineRationalPoint)
  @req P in X "the point needs to lie on the algebraic set"
  J = jacobian_matrix(gens(saturated_ideal(defining_ideal(X))))
  v = coordinates(P)
  JP = map_entries(x->evaluate(x, v), J)
  V = ambient_coordinates(X)
  R = ambient_coordinate_ring(X)
  T = elem_type(R)[ sum([(V[i]-v[i])*JP[i,j] for i in 1:nrows(JP)], init=zero(R)) for j in 1:ncols(JP)]
  return algebraic_set(ideal(R,T), is_radical=true, check=false)
end

@doc raw"""
    tangent_space(P::AbsAffineRationalPoint{<:FieldElem}) -> AlgebraicSet

Return the Zariski tangent space of the ambient scheme of `P` at its point `P`.

See also [`tangent_space(X::AbsAffineScheme{<:Field}, P::AbsAffineRationalPoint)`](@ref)
"""
tangent_space(P::AbsAffineRationalPoint{<:FieldElem}) = tangent_space(codomain(P), P)

@doc raw"""
    is_du_val_singularity(P::AbsAffineRationalPoint{<:Field})

Return whether the ambient scheme of `P` has at most a Du Val singularity at `P`.

Note that this includes the case that ``P`` is a smooth point.
"""
is_du_val_singularity(P::AbsAffineRationalPoint{<:FieldElem}) = is_du_val_singularity(codomain(P), ideal(P))

@doc raw"""
    decide_du_val_singularity(P::AbsAffineRationalPoint{<:Field})

Return whether the ambient scheme of `P` has a Du Val singularity at `P`.

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
  d = decide_du_val_singularity(codomain(P), ideal(P))
  return d[1][1],d[1][3]
end

function stalk(X::AbsAffineScheme, P::AbsAffineRationalPoint)
  P in X || error("not a point of X")
  S = MPolyComplementOfKPointIdeal(P)
  return localization(OO(X),S)
end
