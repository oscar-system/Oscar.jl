struct ToricDivisor
    ambient_variety::NormalToricVarietyType
    polymake_divisor::Polymake.BigObject
end

function pm_tdivisor(td::ToricDivisor)
    return td.polymake_divisor
end


@doc Markdown.doc"""
    ToricDivisor(coeffs::AbstractVector, X::NormalToricVarietyType)

Construct a toric divisor as a sum of the prime divisors on the normal toric
variety `X`.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor([1,0,0,0], H)
A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function ToricDivisor(coeffs::AbstractVector, X::NormalToricVarietyType)
    ptd = Polymake.fulton.TDivisor(COEFFICIENTS=coeffs)
    if length(coeffs) != nprime_divisors(X)
        throw(ArgumentError("Number of coefficients needs to match number of prime divisors!"))
    end
    pmntv = pm_ntv(X)
    Polymake.add(pmntv, "DIVISOR", ptd)
    return ToricDivisor(X, ptd)
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, td::ToricDivisor)
    # fan = get_polyhedral_fan(ntv)
    pmntv = pm_ntv(td.ambient_variety)
    ambdim = pmntv.FAN_AMBIENT_DIM
    print(io, "A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension $(ambdim)")
end



@doc Markdown.doc"""
    isample(td::ToricDivisor) 

Determine whether the toric divisor `td` is ample.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor([1,0,0,0], H)
A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> isample(td)
false
```
"""
isample(td::ToricDivisor) = pm_tdivisor(td).AMPLE::Bool


@doc Markdown.doc"""
    isbasepoint_free(td::ToricDivisor) 

Determine whether the toric divisor `td` is basepoint free.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor([1,0,0,0], H)
A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> isbasepoint_free(td)
true
```
"""
isbasepoint_free(td::ToricDivisor) = pm_tdivisor(td).BASEPOINT_FREE::Bool


@doc Markdown.doc"""
    iscartier(td::ToricDivisor) 

Determine whether the toric divisor `td` is cartier.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor([1,0,0,0], H)
A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> iscartier(td)
true
```
"""
iscartier(td::ToricDivisor) = pm_tdivisor(td).CARTIER::Bool


@doc Markdown.doc"""
    iseffective(td::ToricDivisor) 

Determine whether the toric divisor `td` is effective.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor([1,0,0,0], H)
A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> iseffective(td)
true
```
"""
iseffective(td::ToricDivisor) = pm_tdivisor(td).EFFECTIVE::Bool


@doc Markdown.doc"""
    isintegral(td::ToricDivisor) 

Determine whether the toric divisor `td` is integral.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor([1,0,0,0], H)
A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> isintegral(td)
true
```
"""
isintegral(td::ToricDivisor) = pm_tdivisor(td).INTEGRAL::Bool


@doc Markdown.doc"""
    isnef(td::ToricDivisor) 

Determine whether the toric divisor `td` is nef.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor([1,0,0,0], H)
A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> isnef(td)
true
```
"""
isnef(td::ToricDivisor) = pm_tdivisor(td).NEF::Bool


@doc Markdown.doc"""
    isprincipal(td::ToricDivisor) 

Determine whether the toric divisor `td` is principal.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor([1,0,0,0], H)
A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> isprincipal(td)
false
```
"""
isprincipal(td::ToricDivisor) = pm_tdivisor(td).PRINCIPAL::Bool


@doc Markdown.doc"""
    isq_cartier(td::ToricDivisor) 

Determine whether the toric divisor `td` is Q-Cartier.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor([1,0,0,0], H)
A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> isq_cartier(td)
true
```
"""
isq_cartier(td::ToricDivisor) = pm_tdivisor(td).Q_CARTIER::Bool


@doc Markdown.doc"""
    issemiample(td::ToricDivisor) 

Determine whether the toric divisor `td` is semiample.
"""
issemiample(td::ToricDivisor) = pm_tdivisor(td).SEMIAMPLE::Bool


@doc Markdown.doc"""
    isvery_ample(td::ToricDivisor) 

Determine whether the toric divisor `td` is very ample.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor([1,0,0,0], H)
A torus invariant divisor on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> isvery_ample(td)
false
```
"""
isvery_ample(td::ToricDivisor) = pm_tdivisor(td).VERY_AMPLE::Bool
    
    
