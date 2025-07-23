# For a ranking of jet variables, e.g. u1[123], u3[041], ... one needs:
# - A monomial ordering to compare the elementary symbols u_1, ..., u_n
#   (this is basically a monomial ordering of R[u_1, ..., u_n]
# - An ordering of the multiindices 
# - A decision on which of the two has priority, e.g. position-over-term (pot)
#   or term-over-position (top). Here, the position of the jet variable corresponding
#   to (i, [a_1, ..., a_n]) is the first coordinate, i.e. i.

function Base.isless(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElem}
  @req is_monomial(apre1) && is_monomial(apre2) "Not monomials in comparison"
  
  return 0
end
