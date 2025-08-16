#########################
# 1: The Julia type for ToricDivisorClasses
#########################

@attributes mutable struct ToricDivisorClass
    toric_variety::NormalToricVarietyType
    class::FinGenAbGroupElem
    function ToricDivisorClass(toric_variety::NormalToricVarietyType, class::FinGenAbGroupElem)
        @req parent(class) === class_group_with_map(toric_variety)[1] "The class must belong to the class group of the toric variety"
        return new(toric_variety, class)
    end
end


######################
# 2: Generic constructors
######################

@doc raw"""
    toric_divisor_class(v::NormalToricVarietyType, class::FinGenAbGroupElem)

Construct the toric divisor class associated to a group
element of the class group of the normal toric variety `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> cl = class_group_with_map(P2)[1];

julia> tdc = toric_divisor_class(P2, cl([1]))
Divisor class on a normal toric variety
```
"""
toric_divisor_class(v::NormalToricVarietyType, class::FinGenAbGroupElem) = ToricDivisorClass(v, class)


@doc raw"""
    toric_divisor_class(v::NormalToricVarietyType, coeffs::Vector{T}) where {T <: IntegerUnion}

Construct the toric divisor class associated to a list of integers which
specify an element of the class group of the normal toric variety `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> cl = class_group_with_map(P2)[1];

julia> tdc = toric_divisor_class(P2, cl([ZZRingElem(1)]))
Divisor class on a normal toric variety
```
"""
toric_divisor_class(v::NormalToricVarietyType, coeffs::Vector{T}) where {T <: IntegerUnion} = ToricDivisorClass(v, class_group_with_map(v)[1]([ZZRingElem(c) for c in coeffs]))


######################
# 3: Special constructors
######################

@doc raw"""
    toric_divisor_class(td::ToricDivisor)

Construct the toric divisor class associated to the element ... of the class group of the normal toric variety `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> td = toric_divisor(P2, [1, 2, 3])
Torus-invariant, non-prime divisor on a normal toric variety

julia> tdc = toric_divisor_class(td)
Divisor class on a normal toric variety
```
"""
function toric_divisor_class(td::ToricDivisor)
    f = map_from_torusinvariant_weil_divisor_group_to_class_group(toric_variety(td))
    class = f(sum(coefficients(td) .* gens(domain(f))))
    return toric_divisor_class(toric_variety(td), class)
end


########################
# 4: Addition and scalar multiplication
########################

function Base.:+(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
    @req toric_variety(tdc1) === toric_variety(tdc1) "The divisor classes must be defined on the same toric variety"
    return toric_divisor_class(toric_variety(tdc1), divisor_class(tdc1) + divisor_class(tdc2))
end


function Base.:-(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
    @req toric_variety(tdc1) === toric_variety(tdc1) "The divisor classes must be defined on the same toric variety"
    return toric_divisor_class(toric_variety(tdc1), divisor_class(tdc1) - divisor_class(tdc2))
end

Base.:-(td :: ToricDivisorClass) = toric_divisor_class(toric_variety(td), -divisor_class(td))

Base.:*(c::T, tdc::ToricDivisorClass) where {T <: IntegerUnion} = toric_divisor_class(toric_variety(tdc), ZZRingElem(c) * divisor_class(tdc))


########################
# 5: Equality and hash
########################

function Base.:(==)(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
    return toric_variety(tdc1) === toric_variety(tdc2) && divisor_class(tdc1) == divisor_class(tdc2)
end

function Base.hash(tdc::ToricDivisorClass, h::UInt)
    b = 0x118eb1fba136490c % UInt
    h = hash(toric_variety(tdc), h)
    h = hash(divisor_class(tdc), h)
    return xor(h, b)
end


######################
# 6: Display
######################s

function Base.show(io::IO, td::ToricDivisorClass)
    join(io, "Divisor class on a normal toric variety")
end
