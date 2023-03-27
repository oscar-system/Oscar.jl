########################################################################
# Properties of PrincipalOpenSubsets                                   #
########################################################################
@attr function is_dense(U::PrincipalOpenSubset)
  return !is_zero_divisor(complement_equation(U))
end

