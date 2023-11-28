@doc raw"""
    scheme(td::ToricDivisor)

Every toric divisor has an underlying scheme-theoretic Weil
divisor. This method returns the scheme, on which said
scheme-theoretic divisor is defined. In the case at hand,
this is by design the toric variety at hand, which knows
its underlying scheme. The latter can be accessed with
the function `underlying_toric_structure`.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> td = toric_divisor(P3, [0, 1, 0, 0])
Torus-invariant, prime divisor on a normal toric variety

julia> scheme(td)
Normal toric variety

julia> forget_toric_structure(scheme(td))
(Scheme over QQ covered with 4 patches, Hom: scheme over QQ covered with 4 patches -> normal toric variety)
```
"""
scheme(td::ToricDivisor) = toric_variety(td)


@doc raw"""
    forget_toric_structure(td::ToricDivisor)

Every toric divisor has an underlying scheme-theoretic Weil
divisor. This method returns said divisor.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> td = toric_divisor(P3, [0, 1, 0, 0])
Torus-invariant, prime divisor on a normal toric variety

julia> forget_toric_structure(td)
Effective weil divisor
  on normal, 3-dimensional toric variety
with coefficients in integer ring
given as the formal sum of
  1 * sheaf of ideals
```
"""
function forget_toric_structure(td::ToricDivisor)
  return underlying_divisor(td)::WeilDivisor
end

# For method delegation.
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
