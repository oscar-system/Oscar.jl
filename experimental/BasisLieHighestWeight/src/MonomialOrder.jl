function construct_abs_gen_ordering(
  symb_monomial_ordering::Symbol,
  operators::Vector{RootSpaceElem},
)
  if _is_weighted(symb_monomial_ordering)
    return WSymbOrdering(symb_monomial_ordering, 1:length(operators), [Int(height(alpha)) for alpha in operators])
  else
    return SymbOrdering(symb_monomial_ordering, 1:length(operators))
  end
end
