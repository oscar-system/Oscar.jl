function get_monomial_ordering(
  ordering_input::Union{Symbol,Function},
  ZZx::ZZMPolyRing,
  operators::Vector{RootSpaceElem},
)
  if _is_weighted(ordering_input)
    choosen_monomial_order = monomial_ordering(
      ZZx, ordering_input, [Int(height(alpha)) for alpha in operators]
    )
  else
    choosen_monomial_order = monomial_ordering(ZZx, ordering_input)
  end
  return choosen_monomial_order
end
