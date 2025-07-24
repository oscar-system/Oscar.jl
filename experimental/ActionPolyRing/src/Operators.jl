#######################################
#
#  Unary Operators 
#
#######################################

Base.:-(apre::ActionPolyRingElem) = parent(apre)(-__poly(apre))

#######################################
#
#  Binary Operators 
#
#######################################

function Base.:+(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(__poly(apre1) + __poly(apre2))
end

function Base.:-(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(__poly(apre1) - __poly(apre2))
end

function Base.:*(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(__poly(apre1) * __poly(apre2))
end

function Base.:/(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(__poly(apre1) / __poly(apre2))
end

function divexact(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(divexact(__poly(apre1), __poly(apre2)))
end

#######################################
#
#  Comparison 
#
#######################################

Base.:(==)(dpre1::DifferencePolyRingElem{T}, dpre2::DifferencePolyRingElem{T}) where {T<:RingElem} = __poly(dpre1) == __poly(dpre2) && parent(dpre1) === parent(dpre2)

function Base.hash(dpre::DifferencePolyRingElem, h::UInt)
  b = 0x475b3fa701aa3148 % UInt
  h = hash(parent(dpre), h)
  h = hash(__poly(dpre), h)
  return xor(h, b)
end

