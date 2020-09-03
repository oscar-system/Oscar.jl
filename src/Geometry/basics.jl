

abstract type ProjSpc{T} end
abstract type ProjSpcElem{T} end

module Geometry
using Oscar
import AbstractAlgebra, Nemo
import Base: ==, show, hash

export projective_space, normalize!, weights, coordinate_ring

struct ProjSpc{T}  <: Oscar.ProjSpc{T}
  R::AbstractAlgebra.Ring
  n::Int
  Rx::MPolyRing{T}
end

function Base.show(io::IO, P::ProjSpc)
  w = weights(P)
  if all(isone, w)
    print(io, "Projective space of dim $(P.n) over $(P.R)\n")
  else
    print(io, "Weigthed projective space of dim $(P.n) over $(P.R) and weights $(w)\n")
  end
end

function projective_space(R::AbstractAlgebra.Ring, n::Int, name::Symbol=:x)
  Sx = PolynomialRing(R, name => 0:n)[1]
  Rx = grade(Sx, [1 for i=0:n])
  return ProjSpc(R, n, Rx), gens(Rx)
end

function projective_space(R::AbstractAlgebra.Ring, n::Array{<:Integer, 1}, name::Symbol = :x)
  Sx = PolynomialRing(R, name => 0:length(n)-1)[1]
  Rx = grade(Sx, n)
  return ProjSpc(R, length(n)-1, Rx), gens(Rx)
end

function coordinate_ring(P::ProjSpc)
  return P.Rx
end

function weights(P::ProjSpc)
  return Int[x[1] for x = P.Rx.d]
end

function isweighted(P::ProjSpc)
  return !all(x->isone(x[1]), P.Rx.d)
end

function Oscar.dim(P::ProjSpc)
  return P.n
end

struct ProjSpcElem{T} <: Oscar.ProjSpcElem{T}
  v::Array{T, 1}
  parent::ProjSpc{T}

  function ProjSpcElem(P::ProjSpc{T}, a::Array{T, 1}) where {T}
    @assert length(a) == dim(P)+1
    @assert _isprojective(a)
    r = new{T}(a, P)
    return r
  end
end

function _isprojective(a::Array{T, 1}) where {T <: AbstractAlgebra.FieldElem}
  return !all(iszero, a) 
end

function _isprojective(a::Array{T, 1}) where {T}
  return !all(iszero, a)
end

function _isprojective(a::Array{fmpz, 1})
  all(iszero, a) && return false
  return isone(gcd(a))
end

parent(a::ProjSpcElem) = a.parent
Nemo.elem_type(::ProjSpc{T}) where {T} = ProjSpcElem{T}
Nemo.parent_type(::ProjSpcElem{T}) where {T} = ProjSpc{T}
Nemo.base_ring(P::ProjSpc) = P.R

Base.getindex(a::ProjSpcElem, i::Int) = a.v[i+1]
Base.setindex!(a::ProjSpcElem, v, i::Int) = a.v[i+1] = v

function (P::ProjSpc)(a::Array{<:Any, 1})
  R = base_ring(P)
  @assert length(a) == dim(P)+1
  return ProjSpcElem(P, elem_type(R)[R(x) for x = a])
end

function Base.show(io::IO, a::ProjSpcElem)
  print(io, "(", a[0])
  for i=1:dim(parent(a))
    print(io, " : ", a[i])
  end
  print(io, ")")
end

function ==(a::ProjSpcElem{T}, b::ProjSpcElem{T}) where {T <: AbstractAlgebra.FieldElem}
  i = 0
  w = weights(parent(a))
  while iszero(a[i])
    iszero(b[i]) || return false
    i += 1
  end
  iszero(b[i]) && return false
  ca = a[i]
  cb = b[i]
  ci = i
  i += 1
  while i<= dim(parent(a)) #TODO: weights!!!!
    ca^w[i+1]*b[i]^w[ci+1] == cb^w[i+1]*a[i]^w[ci+1] || return false
    i += 1
  end
  return true
end

function ==(a::ProjSpcElem{fmpz}, b::ProjSpcElem{fmpz})
  i = 0
  while iszero(a[i])
    iszero(b[i]) || return false
    i += 1
  end
  iszero(b[i]) && return false
  ca = a[i]
  cb = b[i]
  (ca == cb || ca == -cb) || return false
  i += 1
  while i<= dim(parent(a)) #TODO: weights!!!!
    ca^w[i+1]*b[i]^w[ci+1] == cb^w[i+1]*a[i]^w[ci+1] || return false
    i += 1
  end
  return true
end


function normalize!(a::ProjSpcElem{T}) where {T <: AbstractAlgebra.FieldElem}
  i = 0
  while iszero(a[i])
    i += 1
  end
  ca = inv(a[i])
  a[i] = one(base_ring(parent(a)))
  w = weights(parent(a))
  any(x->!isone(x), w) && error("cannot normalize with weights")
  i += 1
  while i <= dim(parent(a))
    a[i] = ca*a[i]  #TODO: weights
    i += 1
  end
  return a
end

function normalize!(a::ProjSpcElem{fmpz})
  #projective elements are primitive, ie. the gcd is 1
  #they can be changed only by units, so we can normalize...
  i = 0
  while iszero(a[i])
    i += 1
  end
  w = weights(parent(a))
  any(x->!isone(x), w) && error("cannot normalize with weights")
  ca = sign(a[i])
  while i <= dim(parent(a))
    a[i] = ca*a[i]  #TODO: weights
    i += 1
  end
  return a
end

function Base.hash(a::ProjSpcElem, u::UInt=UInt(123432))
  if isweighted(parent(a))
    return u
  end
  normalize!(a)
  return hash(a.v, u)
end

function Oscar.evaluate(f::MPolyElem, a::ProjSpcElem)
  return evaluate(f, a.v)
end

end

using .Geometry
export projective_space




