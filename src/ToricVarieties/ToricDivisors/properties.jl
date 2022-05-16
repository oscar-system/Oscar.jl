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
istrivial(td::ToricDivisor) = is_principal(td)
export is_principal, istrivial


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
    iseffective(td::ToricDivisor)

Determine whether the toric divisor `td` is effective.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(F4, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> iseffective(td)
true
```
"""
@attr Bool iseffective(td::ToricDivisor) = pm_tdivisor(td).EFFECTIVE
export iseffective


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
