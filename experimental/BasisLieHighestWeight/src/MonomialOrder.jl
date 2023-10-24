function get_monomial_ordering_lt(
  ordering_input::Union{Symbol,Function},
  ZZx::ZZMPolyRing,
  weights_alpha::Vector{Vector{QQFieldElem}},
)::Function
  """
  Returns the desired monomial_ordering function less than, i.e. return true <=> mon1 < mon2
  """
  if isa(ordering_input, Function)
    choosen_monomial_order = ordering_input
  elseif isweighted(ordering_input)
    choosen_monomial_order = monomial_ordering(
      ZZx, ordering_input, Int[Int(sum(w)) for w in weights_alpha]
    )
  else
    choosen_monomial_order = monomial_ordering(ZZx, ordering_input)
  end
  return (mon1::ZZMPolyRingElem, mon2::ZZMPolyRingElem) ->
    (cmp(choosen_monomial_order, mon1, mon2) < 0)
end
