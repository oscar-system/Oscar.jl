function get_monomial_ordering(
  ordering_input::Union{Symbol,Function},
  ZZx::ZZMPolyRing,
  weights_alpha::Vector{Vector{QQFieldElem}},
)
  if _is_weighted(ordering_input)
    choosen_monomial_order = monomial_ordering(
      ZZx, ordering_input, Int[Int(sum(w)) for w in weights_alpha]
    )
  else
    choosen_monomial_order = monomial_ordering(ZZx, ordering_input)
  end
  return choosen_monomial_order
end
