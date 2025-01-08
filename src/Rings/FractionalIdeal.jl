mutable struct FractionalIdeal{S <: Ideal, T <: RingElem}
  num::S
  den::T
end

#TODO get expressify working on all of the ideals
function Base.show(io::IOContext, a::FractionalIdeal)
  print(io, "(")
  show(io, a.num)
  print(io, ")/(")
  show(io, a.den)
  print(io, ")")
end

function fractional_ideal(n::Ideal, d::RingElem)
  base_ring(n) === parent(d) || error("wrong parents")
  iszero(d) && error("zero divisor detected in denominator")
  return FractionalIdeal{typeof(n), typeof(d)}(n, d)
end

function numerator(a::FractionalIdeal)
  return a.num
end

function denominator(b::FractionalIdeal)
  return b.den
end

function base_ring(a::FractionalIdeal)
  return base_ring(numerator(a))
end

base_ring_type(::Type{FractionalIdeal{S, T}}) where {S, T} = base_ring_type(S)

function iszero(a::FractionalIdeal)
  return iszero(numerator(a))
end

function +(a::FractionalIdeal, b::FractionalIdeal)
  return fractional_ideal(a.num*b.den + b.num*a.den, a.den*b.den)
end

function *(a::FractionalIdeal, b::FractionalIdeal)
  return fractional_ideal(a.num*b.num, a.den*b.den)
end

function *(a::FractionalIdeal{S, T}, b::T) where {S <: Ideal, T <: RingElem}
  return fractional_ideal(a.num*b, a.den)
end

function *(b::T, a::FractionalIdeal{S, T}) where {S <: Ideal, T <: RingElem}
  return fractional_ideal(a.num*b, a.den)
end

function ==(a::FractionalIdeal, b::FractionalIdeal)
  return a.num*b.den == b.num*a.den
end

function Base.hash(a::FractionalIdeal, h::UInt)
  b = 0x4680fb583e498597 % UInt
  # there is nothing better to include in the hash
  h = hash(base_ring(a), h)
  return xor(b, h)
end

function ^(a::FractionalIdeal, b::Int)
  return fractional_ideal(a.num^b, a.den^b)
end
