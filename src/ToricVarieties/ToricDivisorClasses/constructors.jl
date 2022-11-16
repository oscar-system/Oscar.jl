#########################
# 1: The Julia type for ToricDivisorClasses
#########################

@attributes mutable struct ToricDivisorClass
    toric_variety::AbstractNormalToricVariety
    class::GrpAbFinGenElem
    function ToricDivisorClass(toric_variety::AbstractNormalToricVariety, class::GrpAbFinGenElem)
        if parent(class) !== class_group(toric_variety)
            throw(ArgumentError("The class must belong to the class group of the toric variety"))
        end
        return new(toric_variety, class)
    end
end
export ToricDivisorClass


######################
# 2: Generic constructors
######################

@doc Markdown.doc"""
    ToricDivisorClass(v::AbstractNormalToricVariety, coeffs::Vector{T}) where {T <: IntegerUnion}

Construct the toric divisor class associated to a list of integers which specify an element of the class group of the normal toric variety `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, class_group(P2)([fmpz(1)]))
A divisor class on a normal toric variety
```
"""
ToricDivisorClass(v::AbstractNormalToricVariety, coeffs::Vector{T}) where {T <: IntegerUnion} = ToricDivisorClass(v, class_group(v)([fmpz(c) for c in coeffs]))


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

julia> td = ToricDivisor(P2, [1, 2, 3])
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

function Base.:+(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
    if toric_variety(tdc1) !== toric_variety(tdc2)
        throw(ArgumentError("The divisor classes must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return ToricDivisorClass(toric_variety(tdc1), divisor_class(tdc1) + divisor_class(tdc2))
end


function Base.:-(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
    if toric_variety(tdc1) !== toric_variety(tdc2)
        throw(ArgumentError("The divisor classes must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return ToricDivisorClass(toric_variety(tdc1), divisor_class(tdc1) - divisor_class(tdc2))
end


Base.:*(c::T, tdc::ToricDivisorClass) where {T <: IntegerUnion} = ToricDivisorClass(toric_variety(tdc), fmpz(c) * divisor_class(tdc))


########################
# 5: Equality
########################

function Base.:(==)(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
    return toric_variety(tdc1) === toric_variety(tdc2) && iszero(divisor_class(tdc1) - divisor_class(tdc2))
end


######################
# 6: Display
######################s

function Base.show(io::IO, td::ToricDivisorClass)
    join(io, "A divisor class on a normal toric variety")
end
