@doc Markdown.doc"""
    iscartier(td::ToricDivisor)

Checks if the divisor `td` is Cartier.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> iscartier(td)
true
```
"""
@attr Bool iscartier(td::ToricDivisor) = pm_tdivisor(td).CARTIER
export iscartier


@doc Markdown.doc"""
    isprincipal(td::ToricDivisor) 

Determine whether the toric divisor `td` is principal.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> isprincipal(td)
false
```
"""
@attr Bool isprincipal(td::ToricDivisor) = pm_tdivisor(td).PRINCIPAL
export isprincipal


@doc Markdown.doc"""
    is_basepoint_free(td::ToricDivisor) 

Determine whether the toric divisor `td` is basepoint free.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
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
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> iseffective(td)
true
```
"""
@attr Bool iseffective(td::ToricDivisor) = pm_tdivisor(td).EFFECTIVE
export iseffective


@doc Markdown.doc"""
    isintegral(td::ToricDivisor) 

Determine whether the toric divisor `td` is integral.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> isintegral(td)
true
```
"""
@attr Bool isintegral(td::ToricDivisor) = pm_tdivisor(td).INTEGRAL
export isintegral


@doc Markdown.doc"""
    isample(td::ToricDivisor) 

Determine whether the toric divisor `td` is ample.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> isample(td)
false
```
"""
@attr Bool isample(td::ToricDivisor) = pm_tdivisor(td).AMPLE
export isample


@doc Markdown.doc"""
    is_very_ample(td::ToricDivisor) 

Determine whether the toric divisor `td` is very ample.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_very_ample(td)
false
```
"""
@attr Bool is_very_ample(td::ToricDivisor) = pm_tdivisor(td).VERY_AMPLE
export is_very_ample


@doc Markdown.doc"""
    isnef(td::ToricDivisor) 

Determine whether the toric divisor `td` is nef.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> isnef(td)
true
```
"""
@attr Bool isnef(td::ToricDivisor) = pm_tdivisor(td).NEF
export isnef


@doc Markdown.doc"""
    is_q_cartier(td::ToricDivisor) 

Determine whether the toric divisor `td` is Q-Cartier.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_q_cartier(td)
true
```
"""
@attr Bool is_q_cartier(td::ToricDivisor) = pm_tdivisor(td).Q_CARTIER
export is_q_cartier


@doc Markdown.doc"""
    isprime(td::ToricDivisor) 

Determine whether the toric divisor `td` is a prime divisor.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> isprime(td)
true
```
"""
@attr Bool function isprime(td::ToricDivisor)
    if sum(coefficients(td)) != 1
        return false
    else
        return all(y -> (y == 1 || y == 0), coefficients(td))
    end
end
export isprime
