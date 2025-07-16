#######################################
#
#  Unary Operators 
#
#######################################

Base.:-(apre::ActionPolyRingElem) = parent(apre)(-apre.p)

#######################################
#
#  Binary Operators 
#
#######################################

function Base.:+(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(apre1.p + apre2.p)
end

function Base.:-(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(apre1.p - apre2.p)
end

function Base.:*(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(apre1.p * apre2.p)
end

function Base.:/(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(apre1.p / apre2.p)
end

function divexact(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(divexact(apre1.p, apre2.p))
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

function Base.isless(dpre1::DifferencePolyRingElem{T}, dpre2::DifferencePolyRingElem{T}) where {T<:RingElem}
  return 0
end
