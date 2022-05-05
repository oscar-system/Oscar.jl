module DiscreteLog

  using GAP
  using Oscar

 ################################################################################
 # Computes the discrete_log of TODO
 # Example:
 #
function disc_log(b::T, x::T) where {T <: FinFieldElem}
  @assert parent(b) === parent(x)
  return disc_log_bs_gs(b, x, order(parent(b)))
end


 end #module DiscreteLog

 using .DiscreteLog
