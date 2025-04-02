@attr Bool is_trivial(tdc::ToricDivisorClass) = iszero(divisor_class(tdc))


@doc raw"""
    is_effective(tdc::ToricDivisorClass)

Determine whether the toric divisor class `tdc` is effective, that is if a toric divisor
in this divisor class is linearly equivalent to an effective toric divisor.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety,2)
Normal toric variety

julia> tdc = toric_divisor_class(P2, [1])
Divisor class on a normal toric variety

julia> is_effective(tdc)
true

julia> tdc2 = toric_divisor_class(P2, [-1])
Divisor class on a normal toric variety

julia> is_effective(tdc2)
false
```
"""
@attr Bool function is_effective(tdc::ToricDivisorClass)
  amb = toric_variety(tdc)
  pi = matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(amb))
  coeffs = coefficients(toric_divisor(tdc))
  P = polyhedron((-identity_matrix(QQ, nrows(pi)), zeros(QQ, nrows(pi))), (transpose(pi), transpose(pi)*coeffs))
  # If the polyhedron is empty, there cannot be a effective representative.
  is_feasible(P) || return false
  # We check whether there is an integral point using a MILP
  milp = mixed_integer_linear_program(P, ones(QQFieldElem, nrows(pi)))
  return !isnothing(optimal_solution(milp))
end
