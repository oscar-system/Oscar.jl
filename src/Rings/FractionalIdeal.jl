export fractional_ideal

import Nemo: numerator, denominator, iszero

mutable struct FractionalIdeal{S <: Ideal, T <: RingElem}
  num::S
  den::T
end

function fractional_ideal(n::Ideal, d::RingElem)
  iszero(d) && error("zero divisor detected in denominator")
  return FractionalIdeal{typeof(n), typeof(d)}(n, d)
end

function numerator(a::FractionalIdeal)
  return a.num
end

function denominator(b::FractionalIdeal)
  return b.den
end

function iszero(a::FractionalIdeal)
  return iszero(numerator(a))
end

function +(a::FractionalIdeal, b::FractionalIdeal)
  return fractional_ideal(a.num*b.den + b.num*a.den, a.den*b.den)
end

function *(a::FractionalIdeal, b::FractionalIdeal)
  return fractional_ideal(a.num*b.num, a.den*b.den)
end

function ==(a::FractionalIdeal, b::FractionalIdeal)
  return a.num*b.den == b.num*a.den
end

function ^(a::FractionalIdeal, b::Int)
  return fractional_ideal(a.num^b, a.den^b)
end
