########################
# 1: The Julia type for toric line bundles
########################

abstract type ToricCoherentSheaf end

@attributes mutable struct ToricLineBundle <: ToricCoherentSheaf
    toric_variety::AbstractNormalToricVariety
    divisor_class::GrpAbFinGenElem
    ToricLineBundle(variety::AbstractNormalToricVariety, divisor_class::GrpAbFinGenElem) = new(variety, divisor_class)
end
export ToricLineBundle


########################
# 2: Generic constructors
########################

@doc Markdown.doc"""
    ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{T}) where {T <: IntegerUnion}

Construct the line bundle on the abstract normal toric variety `v` with class `c`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(v, [fmpz(2)])
A toric line bundle on a normal toric variety
```
"""
function ToricLineBundle(v::AbstractNormalToricVariety, input_class::Vector{T}) where {T <: IntegerUnion}
    class = picard_group(v)(input_class)
    return ToricLineBundle(v, class)
end


########################
# 3: Special constructor
########################

@doc Markdown.doc"""
    ToricLineBundle(v::AbstractNormalToricVariety, d::ToricDivisor)

Construct the toric variety associated to a (Cartier) torus-invariant divisor `d` on the normal toric variety `v`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(v, ToricDivisor(v, [1, 2, 3]))
A toric line bundle on a normal toric variety
```
"""
function ToricLineBundle(v::AbstractNormalToricVariety, d::ToricDivisor)
    if !is_cartier(d)
        throw(ArgumentError("The toric divisor must be Cartier to define a toric line bundle"))
    end
    f = map_from_torusinvariant_cartier_divisor_group_to_picard_group(v)
    class = f(sum(coefficients(d)[i] * gens(domain(f))[i] for i in 1:length(gens(domain(f)))))
    l = ToricLineBundle(v, class)
    set_attribute!(l, :toric_divisor, d)
    return l
end
ToricLineBundle(d::ToricDivisor) = ToricLineBundle(toric_variety(d), d)


########################
# 4: Tensor products
########################

function Base.:*(l1::ToricLineBundle, l2::ToricLineBundle)
    if toric_variety(l1) !== toric_variety(l2)
        throw(ArgumentError("The line bundles must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return ToricLineBundle(toric_variety(l1), divisor_class(l1) + divisor_class(l2))
end
Base.:inv(l::ToricLineBundle) = ToricLineBundle(toric_variety(l), (-1)*divisor_class(l))
Base.:^(l::ToricLineBundle, p::fmpz) = ToricLineBundle(toric_variety(l), p * divisor_class(l))
Base.:^(l::ToricLineBundle, p::Int) = l^fmpz(p)


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
    properties_string = ["A toric"]

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
