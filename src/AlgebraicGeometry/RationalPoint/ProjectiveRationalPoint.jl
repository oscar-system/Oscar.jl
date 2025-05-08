
@doc raw"""
    coordinates(p::AbsProjectiveRationalPoint{S,T}) -> Vector{S}

Return the homogeneous coordinates of the rational point `p`.
"""
coordinates(p::AbsProjectiveRationalPoint)

################################################################################
#
# AbsProjectivePoint Interface
#
################################################################################

parent(p::ProjectiveRationalPoint) = p.parent
coordinates(p::ProjectiveRationalPoint) = p.coordinates

################################################################################

homogeneous_coordinates(p::AbsProjectiveRationalPoint) = coordinates(p)

function dehomogenization(p::AbsProjectiveRationalPoint, i::Int)
  c = coordinates(p)
  s = inv(c[i])
  w = weights(codomain(p))
  any(!is_one, w) && error("cannot dehomogenize with weights")
  return [s*c[j] for j in 1:length(c) if j!=i]
end

function (X::AbsProjectiveScheme)(coordinates::Vector; check::Bool=true)
  k = base_ring(X)
  coordinates = k.(coordinates)
  return ProjectiveRationalPoint(rational_point_set(X), coordinates; check=check)
end

function (X::AbsProjectiveScheme)(p::AbsProjectiveRationalPoint; check::Bool=true)
  codomain(p) === X && return p
  return X(coordinates(p); check=check)
end

function Base.show(io::IO, ::MIME"text/plain", P::AbsProjectiveRationalPoint)
  io = pretty(io)
  println(io, "Projective rational point")
  print(io, Indent())
  println(io, "of ", codomain(P), Dedent())
  print(io, "with coordinates (")
  join(io, coordinates(P), " : ")
  print(io, ")")
end

function Base.show(io::IO, P::AbsProjectiveRationalPoint)
  print(io,"(")
  join(io, coordinates(P), " : ")
  print(io,")")
end

base_ring(P::AbsProjectiveRationalPoint) = base_ring(codomain(P))
coefficient_ring(P::AbsProjectiveRationalPoint) = base_ring(codomain(P))
coordinate(P::AbsProjectiveRationalPoint, i::Int) = coordinates(P)[i]
Base.getindex(P::AbsProjectiveRationalPoint, i::Int) = coordinate(P, i)
Base.setindex!(a::AbsProjectiveRationalPoint, v, i::Int) = coordinates(a)[i] = v

@doc raw"""
    ideal(P::AbsProjectiveRationalPoint)

Return the homogeneous ideal associated to `P`
in the homogeneous coordinate ring of its ambient space.
"""
function ideal(P::AbsProjectiveRationalPoint)
  R = ambient_coordinate_ring(codomain(P))
  V = gens(R)
  w = weights(codomain(P))
  n = length(V)
  m = elem_type(R)[]
  # first non-zero coordinate
  i = 1
  while iszero(P[i])
    i = i+1
  end
  for j in 1:n
    if i==j
      continue
    end
    push!(m, P[i]^w[j]*V[j]^w[i] - P[j]^w[i]*V[i]^w[j])
  end
  return ideal(R, m)
end

function _in(P::AbsProjectiveRationalPoint, X::AbsProjectiveScheme{<:Any,<:MPolyDecRing{<:Any,<:MPolyRing}})
  # X is projective space
  return ambient_space(P) == X
end

function _in(P::AbsProjectiveRationalPoint, X::AbsProjectiveScheme{<:Any,<:MPolyQuoRing{<:MPolyDecRingElem}})
  ambient_space(P) == ambient_space(X) || return false
  c = coordinates(P)
  for f in gens(defining_ideal(X))
    iszero(evaluate(f, c)) || return false
  end
  return true
end

# We can do this without computing the vanishing ideal of X
function _in(P::AbsProjectiveRationalPoint, X::AbsProjectiveAlgebraicSet{<:Any,<:MPolyQuoRing{<:MPolyDecRingElem}})
  ambient_space(P) == ambient_space(X) || return false
  c = coordinates(P)
  for f in gens(fat_ideal(X))
    iszero(evaluate(f, c)) || return false
  end
  return true
end

function Base.in(P::AbsProjectiveRationalPoint, X::AbsProjectiveScheme)
  codomain(P) === X && return true
  return _in(P, X)
end

@doc raw"""
    scheme(P::AbsProjectiveRationalPoint) -> AbsProjectiveScheme

Return the rational point ``P`` viewed as a reduced, projective subscheme
of its ambient projective space.
"""
function scheme(P::AbsProjectiveRationalPoint)
  I = ideal(P)
  R = ambient_coordinate_ring(codomain(P))
  return proj(R, I)
end

function (f::ProjectiveSchemeMor{<:Any,<:Any,<:Any,Nothing})(P::AbsProjectiveRationalPoint)
  @req domain(f) == codomain(P) "$(P) not in domain"
  @req base_ring(domain(f)) == base_ring(codomain(f)) "schemes must be defined over the same base ring"
  x = homogeneous_coordinates(codomain(f))
  g = pullback(f)
  p = coordinates(P)
  imgs = [evaluate(lift(g(y)),p) for y in x]
  return ProjectiveRationalPoint(codomain(f), imgs, check=false)
end

function is_smooth(X::AbsProjectiveScheme{<:Field}, P::AbsProjectiveRationalPoint)
  @req P in X "not a point on X"
  X = codomain(P)
  S = standard_covering(X)
  i = findfirst(!iszero,coordinates(P))
  U = S[i]
  return is_smooth(U,U(dehomogenization(P,i)))
end

is_smooth(P::AbsProjectiveRationalPoint) = is_smooth(codomain(P),P)

function _is_projective(a::Vector{T}) where {T <: AbstractAlgebra.FieldElem}
  return !all(iszero, a)
end

function _is_projective(a::Vector{T}) where {T}
  error("not implemented")
end

function _is_projective(a::Vector{ZZRingElem})
  all(iszero, a) && return false
  return isone(gcd(a))
end

Nemo.parent_type(::Type{AbsProjectiveRationalPoint{S,T}}) where {S,T} = T

function ==(a::AbsProjectiveRationalPoint{S, T}, b::AbsProjectiveRationalPoint{S, U}) where {S<:Union{FieldElem,ZZRingElem},T, U}
  ambient_space(a) == ambient_space(b) || return false
  n = length(coordinates(a))
  i = 1
  w = weights(codomain(a))
  while iszero(a[i])
    iszero(b[i]) || return false
    i += 1
  end
  iszero(b[i]) && return false
  ca = a[i]
  cb = b[i]
  i += 1
  while i <= n
    ca^w[i] * b[i] == cb^w[i] * a[i] || return false
    i += 1
  end
  return true
end

@doc raw"""
    normalize!(a::AbsProjectiveRationalPoint{<:FieldElem})

Normalize `a` such that its first non-zero coordinate is one.
"""
function normalize!(a::AbsProjectiveRationalPoint{<:FieldElem})
  i = 1
  while iszero(a[i])
    i += 1
  end
  ca = inv(a[i])
  a[i] = one(coefficient_ring(a))
  w = weights(codomain(a))
  any(!is_one, w) && error("cannot normalize with weights")
  i += 1
  l = length(w)
  while i <= l
    a[i] = ca^w[i]*a[i]
    i += 1
  end
  return a
end

@doc raw"""
    normalize!(a::AbsProjectiveRationalPoint{ZZRingElem})

Normalize `a` such that its first non-zero coordinate is positive.
"""
function normalize!(a::AbsProjectiveRationalPoint{ZZRingElem})
  #projective elements are primitive, ie. the gcd is 1
  #they can be changed only by units, so we can normalize...
  i = 1
  while iszero(a[i])
    i += 1
  end
  w = weights(codomain(a))
  any(!is_one, w) && error("cannot normalize with weights")
  ca = sign(a[i])
  while i <= dim(codomain(a))
    a[i] = ca*a[i]  #TODO: weights
    i += 1
  end
  return a
end

function Base.hash(a::AbsProjectiveRationalPoint, u::UInt=UInt(123432))
  w = weights(codomain(a))
  if any(!isone, w)
    return u
  end
  normalize!(a)
  return hash(coordinates(a), u)
end
