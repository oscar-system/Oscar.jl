# The scheme on which the toric divisor is defined
scheme(td::ToricDivisor) = toric_variety(td)

# The scheme-theoretic weil divisor corresponding to the toric divisor
function underlying_divisor(td::ToricDivisor; check::Bool=false)
  if has_attribute(td, :underlying_divisor)
    return get_attribute(td, :underlying_divisor)::WeilDivisor
  end
  X = scheme(td)
  generating_divisors = _torusinvariant_weil_divisors(X; check)
  result = sum(a*D for (a, D) in zip(coefficients(td), generating_divisors))
  set_attribute!(td, :underlying_divisor=>result)
  return result
end

# User facing version to get the underlying scheme-theoretic Weil divisor
function forget_toric_structure(td::ToricDivisor)
  return underlying_divisor(td)::WeilDivisor
end
