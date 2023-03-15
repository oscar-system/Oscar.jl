export complement_ideal
export complement_scheme
export ideal_type
export image_ideal
export underlying_morphism

function is_constant(a::MPolyLocRingElem) 
  reduce_fraction(a)
  return is_constant(numerator(a)) && is_constant(denominator(a))
end

### Already implemented in AA -- but probably buggy?
#function is_zero_divisor(f::MPolyRingElem)
#  iszero(f) && return true
#  if is_constant(f)
#    c = AbstractAlgebra.coefficients(f)[1]
#    return is_zero_divisor(c)
#  end
#  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
#end

function is_zero_divisor(f::MPolyLocRingElem)
  iszero(f) && return true
  if is_constant(f)
    c = first(AbstractAlgebra.coefficients(numerator(f)))
    return is_zero_divisor(c)
  end
  return is_zero_divisor(numerator(f))
end

### The following method is only implemented when the coefficient ring is a field.
# The code should be valid generically, but the Singular backend needed for the 
# ideal quotient is probably buggy for non-fields.
function is_zero_divisor(f::MPolyQuoRingElem{<:MPolyRingElem{<:FieldElem}})
  iszero(f) && return true
  b = simplify(f)
  # The next block is basically useless when the coefficient ring is 
  # a field, because it is merely another `is_zero`-check. However, 
  # once more functionality is working, it will actually do stuff and 
  # the above signature can be widened.
  if is_constant(lift(b))
    c = first(AbstractAlgebra.coefficients(lift(b)))
    return is_zero_divisor(c)
  end
  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
end

function is_zero_divisor(f::MPolyQuoLocRingElem{<:Field})
  iszero(f) && return true
  # The next block is basically useless when the coefficient ring is 
  # a field, because it is merely another `is_zero`-check. However, 
  # once more functionality is working, it will actually do stuff and 
  # the above signature can be widened.
  if is_constant(lifted_numerator(f)) && is_constant(lifted_denominator(f))
    c = first(coefficients(lift(numerator(f))))
    return is_zero_divisor(c)
  end
  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
end



