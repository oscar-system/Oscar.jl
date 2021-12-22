@doc Markdown.doc"""
    iscartier(td::ToricDivisor)

Checks if the divisor `td` is Cartier.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> iscartier(td)
true
```
"""
function iscartier(td::ToricDivisor)
    if !has_attribute(td, :iscartier)
        set_attribute!(td, :iscartier, pm_tdivisor(td).CARTIER::Bool)
    end
    return get_attribute(td, :iscartier)
end
export iscartier


@doc Markdown.doc"""
    isprincipal(td::ToricDivisor) 

Determine whether the toric divisor `td` is principal.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isprincipal(td)
false
```
"""
function isprincipal(td::ToricDivisor)
    if !has_attribute(td, :isprincipal)
        set_attribute!(td, :isprincipal, pm_tdivisor(td).PRINCIPAL::Bool)
    end
    return get_attribute(td, :isprincipal)
end
export isprincipal


@doc Markdown.doc"""
    isbasepoint_free(td::ToricDivisor) 

Determine whether the toric divisor `td` is basepoint free.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isbasepoint_free(td)
true
```
"""
function isbasepoint_free(td::ToricDivisor)
    if !has_attribute(td, :isbasepoint_free)
        set_attribute!(td, :isbasepoint_free, pm_tdivisor(td).BASEPOINT_FREE::Bool)
    end
    return get_attribute(td, :isbasepoint_free)
end
export isbasepoint_free


@doc Markdown.doc"""
    iseffective(td::ToricDivisor) 

Determine whether the toric divisor `td` is effective.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> iseffective(td)
true
```
"""
function iseffective(td::ToricDivisor)
    if !has_attribute(td, :iseffective)
        set_attribute!(td, :iseffective, pm_tdivisor(td).EFFECTIVE::Bool)
    end
    return get_attribute(td, :iseffective)
end
export iseffective


@doc Markdown.doc"""
    isintegral(td::ToricDivisor) 

Determine whether the toric divisor `td` is integral.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isintegral(td)
true
```
"""
function isintegral(td::ToricDivisor)
    if !has_attribute(td, :isintegral)
        set_attribute!(td, :isintegral, pm_tdivisor(td).INTEGRAL::Bool)
    end
    return get_attribute(td, :isintegral)
end
export isintegral


@doc Markdown.doc"""
    isample(td::ToricDivisor) 

Determine whether the toric divisor `td` is ample.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isample(td)
false
```
"""
function isample(td::ToricDivisor)
    if !has_attribute(td, :isample)
        set_attribute!(td, :isample, pm_tdivisor(td).AMPLE::Bool)
    end
    return get_attribute(td, :isample)
end
export isample


@doc Markdown.doc"""
    isvery_ample(td::ToricDivisor) 

Determine whether the toric divisor `td` is very ample.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isvery_ample(td)
false
```
"""
function isvery_ample(td::ToricDivisor)
    if !has_attribute(td, :isvery_ample)
        set_attribute!(td, :isvery_ample, pm_tdivisor(td).VERY_AMPLE::Bool)
    end
    return get_attribute(td, :isvery_ample)
end
export isvery_ample


@doc Markdown.doc"""
    isnef(td::ToricDivisor) 

Determine whether the toric divisor `td` is nef.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isnef(td)
true
```
"""
function isnef(td::ToricDivisor)
    if !has_attribute(td, :isnef)
        set_attribute!(td, :isnef, pm_tdivisor(td).NEF::Bool)
    end
    return get_attribute(td, :isnef)
end
export isnef


@doc Markdown.doc"""
    isq_cartier(td::ToricDivisor) 

Determine whether the toric divisor `td` is Q-Cartier.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isq_cartier(td)
true
```
"""
function isq_cartier(td::ToricDivisor)
    if !has_attribute(td, :isq_cartier)
        set_attribute!(td, :isq_cartier, pm_tdivisor(td).Q_CARTIER::Bool)
    end
    return get_attribute(td, :isq_cartier)
end
export isq_cartier


@doc Markdown.doc"""
    isprime_divisor(td::ToricDivisor) 

Determine whether the toric divisor `td` is a prime divisor.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isprime_divisor(td)
true
```
"""
function isprime_divisor(td::ToricDivisor)
    if !has_attribute(td, :isprime_divisor)
        if sum(coefficients(td)) != 1
            set_attribute!(td, :isprime_divisor, false)
        else
            set_attribute!(td, :isprime_divisor, all(y -> (y == 1 || y == 0), coefficients(td)))
        end
    end
    return get_attribute(td, :isprime_divisor)
end
export isprime_divisor
