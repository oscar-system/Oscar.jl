########################
# 1: The Julia type for toric line bundles
########################

abstract type ToricCoherentSheaf end

@attributes mutable struct ToricLineBundle <: ToricCoherentSheaf
    toric_variety::AbstractNormalToricVariety
    divisor_class::GrpAbFinGenElem
    function ToricLineBundle(toric_variety::AbstractNormalToricVariety, class::GrpAbFinGenElem)
        if parent(class) !== picard_group(toric_variety)
            throw(ArgumentError("The class must belong to the Picard group of the toric variety"))
        end
        return new(toric_variety, class)
    end
end


########################
# 2: Generic constructors
########################

@doc raw"""
    toric_line_bundle(v::AbstractNormalToricVariety, class::GrpAbFinGenElem)

Construct the line bundle on the abstract normal toric variety `v` with class `c`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = toric_line_bundle(P2, picard_group(P2)([1]))
Toric line bundle on a normal toric variety
```
"""
toric_line_bundle(v::AbstractNormalToricVariety, class::GrpAbFinGenElem) = ToricLineBundle(v, class)


@doc raw"""
    toric_line_bundle(v::AbstractNormalToricVariety, c::Vector{T}) where {T <: IntegerUnion}

Construct the line bundle on the abstract normal toric variety `v` with class `c`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = toric_line_bundle(v, [ZZRingElem(2)])
Toric line bundle on a normal toric variety
```
"""
function toric_line_bundle(v::AbstractNormalToricVariety, input_class::Vector{T}) where {T <: IntegerUnion}
    class = picard_group(v)(input_class)
    return ToricLineBundle(v, class)
end


########################
# 3: Special constructor
########################

@doc raw"""
    toric_line_bundle(v::AbstractNormalToricVariety, d::ToricDivisor)

Construct the toric variety associated to a (Cartier) torus-invariant divisor `d` on the normal toric variety `v`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = toric_line_bundle(v, toric_divisor(v, [1, 2, 3]))
Toric line bundle on a normal toric variety
```
"""
function toric_line_bundle(v::AbstractNormalToricVariety, d::ToricDivisor)
    if !is_cartier(d)
        throw(ArgumentError("The toric divisor must be Cartier to define a toric line bundle"))
    end
    f = map_from_torusinvariant_cartier_divisor_group_to_picard_group(v)
    class = f(sum(coefficients(d)[i] * gens(domain(f))[i] for i in 1:length(gens(domain(f)))))
    l = ToricLineBundle(v, class)
    set_attribute!(l, :toric_divisor, d)
    return l
end

@doc raw"""
    toric_line_bundle(d::ToricDivisor)

Construct the toric variety associated to a (Cartier) torus-invariant divisor `d`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = toric_divisor(v, [1, 2, 3]);

julia> l = toric_line_bundle(d)
Toric line bundle on a normal toric variety
```
"""
toric_line_bundle(d::ToricDivisor) = toric_line_bundle(toric_variety(d), d)


########################
# 4: Tensor products
########################

function Base.:*(l1::ToricLineBundle, l2::ToricLineBundle)
    if toric_variety(l1) !== toric_variety(l2)
        throw(ArgumentError("The line bundles must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return toric_line_bundle(toric_variety(l1), divisor_class(l1) + divisor_class(l2))
end
Base.:inv(l::ToricLineBundle) = toric_line_bundle(toric_variety(l), (-1)*divisor_class(l))
Base.:^(l::ToricLineBundle, p::ZZRingElem) = toric_line_bundle(toric_variety(l), p * divisor_class(l))
Base.:^(l::ToricLineBundle, p::Int) = l^ZZRingElem(p)


########################
# 5: Equality
########################

function Base.:(==)(l1::ToricLineBundle, l2::ToricLineBundle)
    return toric_variety(l1) === toric_variety(l2) && iszero(divisor_class(l1) - divisor_class(l2))
end


########################
# 6: Display
########################

function Base.show(io::IO, line_bundle::ToricLineBundle)

    # initiate properties string
    properties_string = ["Toric"]

    # collect known properties
    if has_attribute(line_bundle, :toric_divisor)
        td = toric_divisor(line_bundle)
        push_attribute_if_exists!(properties_string, td, :is_principal, "trivial")
        push_attribute_if_exists!(properties_string, td, :is_basepoint_free, "basepoint-free")
        ample_cb!(a, b) = push_attribute_if_exists!(a, b, :is_ample, "ample")
        push_attribute_if_exists!(properties_string, td, :is_very_ample, "very-ample"; callback=ample_cb!)
    end

    # print
    push!(properties_string, "line bundle on a normal toric variety")
    join(io, properties_string, ", ", " ")

end
