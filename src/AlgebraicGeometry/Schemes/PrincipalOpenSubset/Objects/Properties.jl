########################################################################
# Properties of PrincipalOpenSubsets                                   #
########################################################################
@attr Any function is_dense(U::PrincipalOpenSubset)
  return !is_zero_divisor(complement_equation(U))
end

