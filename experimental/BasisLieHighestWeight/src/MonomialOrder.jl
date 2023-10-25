function get_monomial_ordering_lt(
  ordering_input::Union{Symbol,Function}, ZZx::ZZMPolyRing
)::Function
  """
  Returns the desired monomial_ordering function less than, i.e. return true <=> mon1 < mon2
  """
  if isa(ordering_input, Function)
    choosen_monomial_order = ordering_input
  else
    choosen_monomial_order = monomial_ordering(ZZx, ordering_input)
  end
  return (mon1::ZZMPolyRingElem, mon2::ZZMPolyRingElem) ->
    (cmp(choosen_monomial_order, mon1, mon2) < 0)
end
