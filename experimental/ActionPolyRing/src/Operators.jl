#######################################
#
#  Unary Operators 
#
#######################################

Base.:-(apre::ActionPolyRingElem) = parent(apre)(-data(apre))

#######################################
#
#  Binary Operators 
#
#######################################

function Base.:+(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(data(apre1) + data(apre2))
end

function Base.:-(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(data(apre1) - data(apre2))
end

function Base.:*(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(data(apre1) * data(apre2))
end

function Base.:/(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(data(apre1) / data(apre2))
end

function divexact(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(divexact(data(apre1), data(apre2)))
end

function mod(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(mod(data(apre1), data(apre2)))
end

function gcd(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2) 
  return parent(apre1)(gcd(data(apre1), data(apre2)))
end

function lcm(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
   check_parent(apre1, apre2)
   return parent(apre1)(lcm(data(apre1), data(apre2)))
end

#######################################
#
#  Comparison 
#
#######################################

Base.:(==)(dpre1::DifferencePolyRingElem{T}, dpre2::DifferencePolyRingElem{T}) where {T<:RingElement} = data(dpre1) == data(dpre2) && parent(dpre1) === parent(dpre2)

function Base.hash(dpre::DifferencePolyRingElem, h::UInt)
  b = 0x475b3fa701aa3148 % UInt
  h = hash(parent(dpre), h)
  h = hash(data(dpre), h)
  return xor(h, b)
end

function Base.isless(dpre1::DifferencePolyRingElem{T}, dpre2::DifferencePolyRingElem{T}) where {T<:RingElement}
  check_parent(dpre1, dpre2)
  dpr = parent(dpre1)
  vtj = __vtj(dpr)
  @req haskey(vtj, dpre1) && haskey(vtj, dpre2) "Not jet variables in comparison"
  m = nelementary_symbols(dpr)
  ind1, ind2 = vtj[dpre1], vtj[dpre2]
  v1, v2 = vcat([ind1[1] == j ? one(ZZ) : zero(ZZ) for j in 1:m], ind1[2]), vcat([ind2[1] == j ? one(ZZ) : zero(ZZ) for j in 1:m], ind2[2])
  M = riquier_matrix(ranking(dpr))
  return isless(M * v1, M * v2)
end

#######################################
#
#  Unsafe functions 
#
#######################################

function zero!(a::DifferencePolyRingElem{T}) where {T <: RingElement}
  a.p = zero!(a.p)
  return a
end

function one!(a::DifferencePolyRingElem{T}) where {T <: RingElement}
  a.p = one!(a.p)
  return a
end

function neg!(z::DifferencePolyRingElem{T}, a::DifferencePolyRingElem{T}) where {T <: RingElement}
  if parent(data(z)) == parent(data(a))
    z.p = neg!(z.p, a.p)
  else
    z.p = -a.p
  end
  return z
end

function fit!(a::DifferencePolyRingElem, n::Int)
  fit!(data(a), n)
end

function add!(a::DifferencePolyRingElem{T}, b::DifferencePolyRingElem{T}, c::DifferencePolyRingElem{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = add!(data(a), data(b), data(c))
  else
    a.p = data(b + c)
  end
  return a
end

function add!(a::DifferencePolyRingElem{T}, b::DifferencePolyRingElem{T}, c::RingElement) where {T <: RingElement}
  if parent(data(a)) == parent(data(b))
    a.p = add!(data(a), data(b), c)
  else
    a.p = data(b + c)
  end
  return a
end

add!(a::DifferencePolyRingElem{T}, b::RingElement, c::DifferencePolyRingElem{T}) where {T <: RingElement} = add!(a, c, b)

function sub!(a::DifferencePolyRingElem{T}, b::DifferencePolyRingElem{T}, c::DifferencePolyRingElem{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = sub!(data(a), data(b), data(c))
  else
    a.p = data(b - c)
  end
  return a
end

function sub!(a::DifferencePolyRingElem{T}, b::DifferencePolyRingElem{T}, c::RingElement) where {T <: RingElement}
  if parent(data(a)) == parent(data(b))
    a.p = sub!(data(a), data(b), c)
  else
    a.p = data(b - c)
  end
  return a
end

function sub!(a::DifferencePolyRingElem{T}, b::RingElement, c::DifferencePolyRingElem{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(c))
    a.p = sub!(data(a), b, data(c))
  else
    a.p = data(b - c)
  end
  return a
end

function mul!(a::DifferencePolyRingElem{T}, b::DifferencePolyRingElem{T}, c::DifferencePolyRingElem{T}) where {T <: RingElement}
  if parent(data(a)) == parent(data(b)) == parent(data(c))
    a.p = mul!(data(a), data(b), data(c))
  else
    a.p = data(b * c)
  end
  return a
end

function mul!(a::DifferencePolyRingElem{T}, b::DifferencePolyRingElem{T}, c::RingElement) where {T <: RingElement}
  if parent(data(a)) == parent(data(b))
    a.p = mul!(data(a), data(b), c)
  else
    a.p = data(b * c)
  end
  return a
end

mul!(a::DifferencePolyRingElem{T}, b::RingElement, c::DifferencePolyRingElem{T}) where {T <: RingElement} = mul!(a, c, b)

