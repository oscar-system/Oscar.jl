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
function iscartier(td::ToricDivisor)
    return get_attribute!(td, :iscartier) do
        return pm_tdivisor(td).CARTIER
    end::Bool
end
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
function isprincipal(td::ToricDivisor)
    return get_attribute!(td, :isprincipal) do
        return pm_tdivisor(td).PRINCIPAL
    end::Bool
end
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
function is_basepoint_free(td::ToricDivisor)
    return get_attribute!(td, :is_basepoint_free) do
        return pm_tdivisor(td).BASEPOINT_FREE::Bool
    end
end
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
function iseffective(td::ToricDivisor)
    return get_attribute!(td, :iseffective) do
        return pm_tdivisor(td).EFFECTIVE
    end::Bool
end
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
function isintegral(td::ToricDivisor)
    return get_attribute!(td, :isintegral) do
        return pm_tdivisor(td).INTEGRAL
    end::Bool
end
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
function isample(td::ToricDivisor)
    return get_attribute!(td, :isample) do
        return pm_tdivisor(td).AMPLE
    end::Bool
end
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
function is_very_ample(td::ToricDivisor)
    return get_attribute!(td, :is_very_ample) do
        return pm_tdivisor(td).VERY_AMPLE::Bool
    end
end
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
function isnef(td::ToricDivisor)
    return get_attribute!(td, :isnef) do
        return pm_tdivisor(td).NEF
    end::Bool
end
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
function is_q_cartier(td::ToricDivisor)
    return get_attribute!(td, :is_q_cartier) do
        return pm_tdivisor(td).Q_CARTIER::Bool
    end
end
export is_q_cartier


@doc Markdown.doc"""
    is_prime_divisor(td::ToricDivisor) 

Determine whether the toric divisor `td` is a prime divisor.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus-invariant, prime divisor on a normal toric variety

julia> is_prime_divisor(td)
true
```
"""
function is_prime_divisor(td::ToricDivisor)
    return get_attribute!(td, :is_prime_divisor) do    
        if sum(coefficients(td)) != 1
            return false
        else
            return all(y -> (y == 1 || y == 0), coefficients(td))
        end
    end::Bool
end
export is_prime_divisor
