@doc Markdown.doc"""
    iscartier(td::ToricDivisor)

Checks if the divisor `td` is Cartier.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> iscartier(td)
true
```
"""
iscartier(td::ToricDivisor) = pm_tdivisor(td).CARTIER::Bool
export iscartier


@doc Markdown.doc"""
    isprincipal(td::ToricDivisor) 

Determine whether the toric divisor `td` is principal.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isprincipal(td)
false
```
"""
isprincipal(td::ToricDivisor) = pm_tdivisor(td).PRINCIPAL::Bool
export isprincipal


@doc Markdown.doc"""
    isbasepoint_free(td::ToricDivisor) 

Determine whether the toric divisor `td` is basepoint free.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isbasepoint_free(td)
true
```
"""
isbasepoint_free(td::ToricDivisor) = pm_tdivisor(td).BASEPOINT_FREE::Bool
export isbasepoint_free


@doc Markdown.doc"""
    iseffective(td::ToricDivisor) 

Determine whether the toric divisor `td` is effective.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> iseffective(td)
true
```
"""
iseffective(td::ToricDivisor) = pm_tdivisor(td).EFFECTIVE::Bool
export iseffective


@doc Markdown.doc"""
    isintegral(td::ToricDivisor) 

Determine whether the toric divisor `td` is integral.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isintegral(td)
true
```
"""
isintegral(td::ToricDivisor) = pm_tdivisor(td).INTEGRAL::Bool
export isintegral


@doc Markdown.doc"""
    isample(td::ToricDivisor) 

Determine whether the toric divisor `td` is ample.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isample(td)
false
```
"""
isample(td::ToricDivisor) = pm_tdivisor(td).AMPLE::Bool
export isample


@doc Markdown.doc"""
    isvery_ample(td::ToricDivisor) 

Determine whether the toric divisor `td` is very ample.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isvery_ample(td)
false
```
"""
isvery_ample(td::ToricDivisor) = pm_tdivisor(td).VERY_AMPLE::Bool
export isvery_ample


@doc Markdown.doc"""
    isnef(td::ToricDivisor) 

Determine whether the toric divisor `td` is nef.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isnef(td)
true
```
"""
isnef(td::ToricDivisor) = pm_tdivisor(td).NEF::Bool
export isnef


@doc Markdown.doc"""
    isq_cartier(td::ToricDivisor) 

Determine whether the toric divisor `td` is Q-Cartier.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isq_cartier(td)
true
```
"""
isq_cartier(td::ToricDivisor) = pm_tdivisor(td).Q_CARTIER::Bool
export isq_cartier


@doc Markdown.doc"""
    isprime_divisor(td::ToricDivisor) 

Determine whether the toric divisor `td` is a prime divisor.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isprime_divisor(td)
true
```
"""
function isprime_divisor(td::ToricDivisor)
    if sum(coefficients(td)) != 1
        return false
    end
    return all(y -> (y == 1 || y == 0), coefficients(td))
end
export isprime_divisor
