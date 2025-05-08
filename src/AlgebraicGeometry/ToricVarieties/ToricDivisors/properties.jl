@doc raw"""
    is_cartier(td::ToricDivisor)

Check if the divisor `td` is Cartier.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> td = toric_divisor(F4, [1,0,0,0])
Torus-invariant, prime divisor on a normal toric variety

julia> is_cartier(td)
true
```
"""
@attr Bool function is_cartier(td::ToricDivisor)
  try
    pm_object(td).CARTIER
  catch
    rethrow(ArgumentError("polymake could not determine whether divisor is Cartier"))
  end
  @req !isnothing(pm_object(td).CARTIER) "polymake could not determine whether divisor is Cartier"
  return pm_object(td).CARTIER::Bool
end


@doc raw"""
    is_principal(td::ToricDivisor)

Determine whether the toric divisor `td` is principal.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> td = toric_divisor(F4, [1,0,0,0])
Torus-invariant, prime divisor on a normal toric variety

julia> is_principal(td)
false
```
"""
@attr Bool is_principal(td::ToricDivisor) = pm_object(td).PRINCIPAL

@attr Bool is_trivial(td::ToricDivisor) = all([c == 0 for c in coefficients(td)])


@doc raw"""
    is_basepoint_free(td::ToricDivisor)

Determine whether the toric divisor `td` is basepoint free.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> td = toric_divisor(F4, [1,0,0,0])
Torus-invariant, prime divisor on a normal toric variety

julia> is_basepoint_free(td)
true
```
"""
@attr Bool is_basepoint_free(td::ToricDivisor) = pm_object(td).BASEPOINT_FREE


@doc raw"""
    is_effective(td::ToricDivisor)

Determine whether the toric divisor `td` is effective,
i.e. if all of its coefficients are non-negative.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety,2)
Normal toric variety

julia> td = toric_divisor(P2, [1,-1,0])
Torus-invariant, non-prime divisor on a normal toric variety

julia> is_effective(td)
false

julia> td2 = toric_divisor(P2, [1,2,3])
Torus-invariant, non-prime divisor on a normal toric variety

julia> is_effective(td2)
true
```
"""
@attr Bool is_effective(td::ToricDivisor) = all(>=(0), coefficients(td))


@doc raw"""
    is_integral(td::ToricDivisor)

Determine whether the toric divisor `td` is integral.
# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> td = toric_divisor(F4, [1,0,0,0])
Torus-invariant, prime divisor on a normal toric variety

julia> is_integral(td)
true
```
"""
@attr Bool is_integral(td::ToricDivisor) = pm_object(td).INTEGRAL


@doc raw"""
    is_ample(td::ToricDivisor)

Determine whether the toric divisor `td` is ample.
# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> td = toric_divisor(F4, [1,0,0,0])
Torus-invariant, prime divisor on a normal toric variety

julia> is_ample(td)
false
```
"""
@attr Bool function is_ample(td::ToricDivisor)
  @req is_complete(toric_variety(td)) "Ampleness test is only implemented for complete toric varieties. ([Def 6.1.9 CLS11])"
  @req is_cartier(td) "Definition of ample divisor requires divisor to be Cartier. ([Def 6.1.9 CLS11])"
  return pm_object(td).AMPLE::Bool
end


@doc raw"""
    is_very_ample(td::ToricDivisor)

Determine whether the toric divisor `td` is very ample.
# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> td = toric_divisor(F4, [1,0,0,0])
Torus-invariant, prime divisor on a normal toric variety

julia> is_very_ample(td)
false
```
"""
@attr Bool function is_very_ample(td::ToricDivisor)
  @req is_complete(toric_variety(td)) "Very ampleness test is only implemented for complete toric varieties. ([Def 6.1.9 CLS11])"
  return pm_object(td).VERY_AMPLE::Bool
end


@doc raw"""
    is_nef(td::ToricDivisor)

Determine whether the toric divisor `td` is nef.
# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> td = toric_divisor(F4, [1,0,0,0])
Torus-invariant, prime divisor on a normal toric variety

julia> is_nef(td)
true
```
"""
@attr Bool is_nef(td::ToricDivisor) = pm_object(td).NEF


@doc raw"""
    is_q_cartier(td::ToricDivisor)

Determine whether the toric divisor `td` is Q-Cartier.
# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> td = toric_divisor(F4, [1,0,0,0])
Torus-invariant, prime divisor on a normal toric variety

julia> is_q_cartier(td)
true
```
"""
@attr Bool is_q_cartier(td::ToricDivisor) = pm_object(td).Q_CARTIER


@doc raw"""
    is_prime(td::ToricDivisor)

Determine whether the toric divisor `td` is a prime divisor.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> td = toric_divisor(F4, [1,0,0,0])
Torus-invariant, prime divisor on a normal toric variety

julia> is_prime(td)
true
```
"""
@attr Bool function is_prime(td::ToricDivisor)
    if sum(coefficients(td)) != 1
        return false
    else
        return all(y -> is_zero(y) || is_one(y), coefficients(td))
    end
end
