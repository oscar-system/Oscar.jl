#######################################
#
#  Unary Operators 
#
#######################################

Base.:-(apre::ActionPolyRingElem) = parent(apre)(-apre.upoly_ring_elem)

#######################################
#
#  Binary Operators 
#
#######################################

function Base.:+(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(apre1.upoly_ring_elem + apre2.upoly_ring_elem)
end

function Base.:*(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  check_parent(apre1, apre2)
  return parent(apre1)(apre1.upoly_ring_elem * apre2.upoly_ring_elem)
end
