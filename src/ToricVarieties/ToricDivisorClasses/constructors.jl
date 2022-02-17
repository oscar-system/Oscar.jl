#########################
# 1: The Julia type for ToricDivisorClasses
#########################

@attributes mutable struct ToricDivisorClass
    toric_variety::AbstractNormalToricVariety
    class::GrpAbFinGenElem
    function ToricDivisorClass(toric_variety::AbstractNormalToricVariety, class::GrpAbFinGenElem)
        if parent(class) !== class_group(toric_variety)
            throw(ArgumentError("The class must belong to the class group of the toric variety."))
        end
        return new(toric_variety, class)
    end
end
export ToricDivisorClass


######################
# 2: Generic constructors
######################

@doc Markdown.doc"""
    ToricDivisorClass(v::AbstractNormalToricVariety, coeffs::Vector{fmpz})

Construct the toric divisor class associated to a list of integers which specify an element of the class group of the normal toric variety `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, class_group(P2)([fmpz(1)]))
A divisor class on a normal toric variety
```
"""
ToricDivisorClass(v::AbstractNormalToricVariety, coeffs::Vector{fmpz}) = ToricDivisorClass(v, class_group(v)(coeffs))

@doc Markdown.doc"""
    ToricDivisorClass(v::AbstractNormalToricVariety, coeffs::Vector{Int})

Construct the toric divisor class associated to a list of integers which specify an element of the class group of the normal toric variety `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, class_group(P2)([1]))
A divisor class on a normal toric variety
```
"""
ToricDivisorClass(v::AbstractNormalToricVariety, coeffs::Vector{Int}) = ToricDivisorClass(v, class_group(v)([fmpz(c) for c in coeffs]))


######################
# 3: Special constructors
######################

@doc Markdown.doc"""
    ToricDivisorClass(td::ToricDivisor)

Construct the toric divisor class associated to the element ... of the class group of the normal toric variety `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> td = ToricDivisor(P2, [1,2,3])
A torus-invariant, non-prime divisor on a normal toric variety

julia> tdc = ToricDivisorClass(td)
A divisor class on a normal toric variety
```
"""
function ToricDivisorClass(td::ToricDivisor)
    f = map_from_torusinvariant_weil_divisor_group_to_class_group(toric_variety(td))
    class = f(sum(coefficients(td) .* gens(domain(f))))
    return ToricDivisorClass(toric_variety(td), class)
end


########################
# 4: Addition and scalar multiplication
########################

@doc Markdown.doc"""
    Base.:+(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)

Return the sum of the toric divisor classes `tdc1` and `tdc2`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc1 = ToricDivisorClass(P2, [1])
A divisor class on a normal toric variety

julia> tdc2 = ToricDivisorClass(P2, [2])
A divisor class on a normal toric variety

julia> tdc1+tdc2
A divisor class on a normal toric variety
```
"""
function Base.:+(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
    # check input
    if toric_variety(tdc1) !== toric_variety(tdc2)
        throw(ArgumentError("The toric classes must be defined on identically the same toric variety."))
    end
    
    # return the new divisor class
    return ToricDivisorClass(toric_variety(tdc1), divisor_class(tdc1) + divisor_class(tdc2))
end


@doc Markdown.doc"""
    Base.:-(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)

Return the difference of the toric divisor classes `tdc1` and `tdc2`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc1 = ToricDivisorClass(P2, [1])
A divisor class on a normal toric variety

julia> tdc2 = ToricDivisorClass(P2, [2])
A divisor class on a normal toric variety

julia> tdc1-tdc2
A divisor class on a normal toric variety
```
"""
function Base.:-(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
    # check input
    if toric_variety(tdc1) !== toric_variety(tdc2)
        throw(ArgumentError("The toric classes must be defined on identically the same toric variety."))
    end
    
    # return the new divisor class
    return ToricDivisorClass(toric_variety(tdc1), divisor_class(tdc1) - divisor_class(tdc2))
end


@doc Markdown.doc"""
    Base.:*(c::fmpz, td::ToricDivisorClass)

Return `c`-times the toric divisor `tdc`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, [1])
A divisor class on a normal toric variety

julia> fmpz(2)*tdc
A divisor class on a normal toric variety
```
"""
Base.:*(c::fmpz, tdc::ToricDivisorClass) = ToricDivisorClass(toric_variety(tdc), c * divisor_class(tdc))
Base.:*(c::Int, tdc::ToricDivisorClass) = ToricDivisorClass(toric_variety(tdc), fmpz(c) * divisor_class(tdc))


########################
# 5: Equality
########################

@doc Markdown.doc"""
    Base.:(==)(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)

Returns true if the toric divisor classes `tdc1` and `tdc2` are equal and false otherwise.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc1 = ToricDivisorClass(P2, [1])
A divisor class on a normal toric variety

julia> tdc2 = ToricDivisorClass(P2, [2])
A divisor class on a normal toric variety

julia> tdc1 == tdc2
false
```
"""
function Base.:(==)(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
    if toric_variety(tdc1) !== toric_variety(tdc2)
        return false
    end
    return iszero(divisor_class(tdc1) - divisor_class(tdc2))
end


######################
# 6: Display
######################s

function Base.show(io::IO, td::ToricDivisorClass)
    join(io, "A divisor class on a normal toric variety")
end
