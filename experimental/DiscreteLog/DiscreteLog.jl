module DiscreteLog

  using GAP
  using Oscar

################################################################################
# Input: b::T and x::T where {T <: FinFieldElem} and b and x are elements from
#        the same finite field
# Output: Integer s such that b^s = x
#         If no such x exists, an error is thrown
# Example:
# julia> F = GF(3,4);
# julia> a = gen(F)^21;
# julia> Oscar.DiscreteLog.disc_log(gen(F), a)
# > 21
#
function disc_log(b::T, x::T) where {T <: FinFieldElem}
  @assert parent(b) === parent(x)
  return disc_log_bs_gs(b, x, order(parent(b)))
end

end #module DiscreteLog

using .DiscreteLog
