@doc Markdown.doc"""
    is_cartier(td::ToricDivisor)

Checks if the divisor `td` is Cartier.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(F4, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_cartier(td)
true
```
"""
@attr Bool is_cartier(td::ToricDivisor) = pm_tdivisor(td).CARTIER
export is_cartier


@doc Markdown.doc"""
    is_principal(td::ToricDivisor)

Determine whether the toric divisor `td` is principal.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(F4, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_principal(td)
false
```
"""
@attr Bool is_principal(td::ToricDivisor) = pm_tdivisor(td).PRINCIPAL
export is_principal

@attr Bool is_trivial(td::ToricDivisor) = all([c == 0 for c in coefficients(td)])
export is_trivial


@doc Markdown.doc"""
    is_basepoint_free(td::ToricDivisor)

Determine whether the toric divisor `td` is basepoint free.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(F4, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_basepoint_free(td)
true
```
"""
@attr Bool is_basepoint_free(td::ToricDivisor) = pm_tdivisor(td).BASEPOINT_FREE
export is_basepoint_free


@doc Markdown.doc"""
    is_effective(td::ToricDivisor)

Determine whether the toric divisor `td` is effective,
i.e. if all of its coefficients are non-negative.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety,2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(P2, [1,-1,0])
A torus-invariant, non-prime divisor on a normal toric variety

julia> is_effective(td)
false

julia> td2 = ToricDivisor(P2, [1,2,3])
A torus-invariant, non-prime divisor on a normal toric variety

julia> is_effective(td2)
true
```
"""
@attr Bool is_effective(td::ToricDivisor) = all(c -> (c >= 0), coefficients(td))
export is_effective


@doc Markdown.doc"""
    is_linearly_equivalent_to_effective_toric_divisor(td::ToricDivisor)

Determine whether the toric divisor `td` is linearly equivalent to an effective toric divisor.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety,2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(P2, [1,-1,0])
A torus-invariant, non-prime divisor on a normal toric variety

julia> is_effective(td)
false

julia> is_linearly_equivalent_to_effective_toric_divisor(td)
true
```
"""
@attr Bool is_linearly_equivalent_to_effective_toric_divisor(td::ToricDivisor) = pm_tdivisor(td).EFFECTIVE
export is_linearly_equivalent_to_effective_toric_divisor


@doc Markdown.doc"""
    is_integral(td::ToricDivisor)

Determine whether the toric divisor `td` is integral.
# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(F4, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_integral(td)
true
```
"""
@attr Bool is_integral(td::ToricDivisor) = pm_tdivisor(td).INTEGRAL
export is_integral


@doc Markdown.doc"""
    is_ample(td::ToricDivisor)

Determine whether the toric divisor `td` is ample.
# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(F4, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_ample(td)
false
```
"""
@attr Bool is_ample(td::ToricDivisor) = pm_tdivisor(td).AMPLE
export is_ample


@doc Markdown.doc"""
    is_very_ample(td::ToricDivisor)

Determine whether the toric divisor `td` is very ample.
# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(F4, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_very_ample(td)
false
```
"""
@attr Bool is_very_ample(td::ToricDivisor) = pm_tdivisor(td).VERY_AMPLE
export is_very_ample


@doc Markdown.doc"""
    is_nef(td::ToricDivisor)

Determine whether the toric divisor `td` is nef.
# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(F4, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_nef(td)
true
```
"""
@attr Bool is_nef(td::ToricDivisor) = pm_tdivisor(td).NEF
export is_nef


@doc Markdown.doc"""
    is_q_cartier(td::ToricDivisor)

Determine whether the toric divisor `td` is Q-Cartier.
# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(F4, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_q_cartier(td)
true
```
"""
@attr Bool is_q_cartier(td::ToricDivisor) = pm_tdivisor(td).Q_CARTIER
export is_q_cartier


@doc Markdown.doc"""
    is_prime(td::ToricDivisor)

Determine whether the toric divisor `td` is a prime divisor.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(F4, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_prime(td)
true
```
"""
@attr Bool function is_prime(td::ToricDivisor)
    if sum(coefficients(td)) != 1
        return false
    else
        return all(y -> (y == 1 || y == 0), coefficients(td))
    end
end
export is_prime
