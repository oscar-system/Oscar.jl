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
    ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{fmpz})

Construct the line bundle on the abstract normal toric variety `v` with class `c`.
"""
function ToricLineBundle(v::AbstractNormalToricVariety, input_class::Vector{fmpz})
    class = picard_group(v)(input_class)
    return ToricLineBundle(v, class)
end

@doc Markdown.doc"""
    ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{Int})

Convenience method for ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{fmpz}).

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle( v, [ 2 ] )
A toric line bundle on a normal toric variety
```
"""
function ToricLineBundle(v::AbstractNormalToricVariety, input_class::Vector{Int})
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

julia> l = ToricLineBundle(v, ToricDivisor(v,[1,2,3]))
A toric line bundle on a normal toric variety
```
"""
function ToricLineBundle(v::AbstractNormalToricVariety, d::ToricDivisor)
    if !iscartier(d)
        throw(ArgumentError("The toric divisor must be Cartier to define a toric line bundle."))
    end
    f = map_from_torusinvariant_cartier_divisor_group_to_picard_group(v)
    class = f(sum([coefficients(d)[i] * gens(domain(f))[i] for i in 1:length(gens(domain(f)))]))
    return ToricLineBundle(v, class)
end
ToricLineBundle(d::ToricDivisor) = ToricLineBundle(toric_variety(d), d)


########################
# 4: Tensor products
########################

@doc Markdown.doc"""
    Base.:*(l1::ToricLineBundle, l2::ToricLineBundle)

Return the tensor product of two line bundles `l1` and `l2`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l1 = ToricLineBundle(P2, [1])
A toric line bundle on a normal toric variety

julia> l2 = ToricLineBundle(P2, [2])
A toric line bundle on a normal toric variety

julia> l1*l2
A toric line bundle on a normal toric variety
```
"""
function Base.:*(l1::ToricLineBundle, l2::ToricLineBundle)
    # check input
    if toric_variety(l1) !== toric_variety(l2)
        throw(ArgumentError("The line bundles must be defined on identically the same toric variety."))
    end
    
    # return the new divisor class
    return ToricLineBundle(toric_variety(l1), divisor_class(l1) + divisor_class(l2))
end


@doc Markdown.doc"""
    Base.:inv(l::ToricLineBundle)

Returns the inverse of the toric line bundle `l`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(P2, [1])
A toric line bundle on a normal toric variety

julia> inv(l)
A toric line bundle on a normal toric variety
```
"""
Base.:inv(l::ToricLineBundle) = ToricLineBundle(toric_variety(l), (-1)*divisor_class(l))


@doc Markdown.doc"""
    Base.:^(l::ToricLineBundle, p::fmpz)

Return the `p`-th tensor power of the toric line bundle `l`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(P2, [1])
A toric line bundle on a normal toric variety

julia> l^(-1)
A toric line bundle on a normal toric variety
```
"""
Base.:^(l::ToricLineBundle, p::fmpz) = ToricLineBundle(toric_variety(l), p * divisor_class(l))
Base.:^(l::ToricLineBundle, p::Int) = l^fmpz(p)


########################
# 4: Special constructors
########################

# The following is both a special constructor of a toric line bundle and an attribute of a toric variety.
# We hence provide snake case and Camel case methods to call it.

@doc Markdown.doc"""
    StructureSheaf(v::AbstractNormalToricVariety)

Construct the structure sheaf of a normal toric variety.
For convenience, we also support `structure_sheaf(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> StructureSheaf(v)
A toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle function StructureSheaf(v::AbstractNormalToricVariety)
    class = zero(picard_group(v))
    return ToricLineBundle(v, class)
end
structure_sheaf(v::AbstractNormalToricVariety) = StructureSheaf(v)
export StructureSheaf, structure_sheaf


@doc Markdown.doc"""
    AnticanonicalBundle(v::AbstractNormalToricVariety)

Construct the anticanonical bundle of a normal toric variety.
For convenience, we also support `anticanonical_bundle(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> AnticanonicalBundle(v)
A toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle function AnticanonicalBundle(v::AbstractNormalToricVariety)
    return ToricLineBundle(v, sum(cox_ring(v).d))
end
anticanonical_bundle(v::AbstractNormalToricVariety) = AnticanonicalBundle(v)
export AnticanonicalBundle, anticanonical_bundle


@doc Markdown.doc"""
    CanonicalBundle(v::AbstractNormalToricVariety)

Construct the canonical bundle of a normal toric variety.
For convenience, we also support `canonical_bundle(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> CanonicalBundle(v)
A toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle function CanonicalBundle(v::AbstractNormalToricVariety)
    return ToricLineBundle(v, (-1)*sum(cox_ring(v).d))
end
canonical_bundle(v::AbstractNormalToricVariety) = CanonicalBundle(v)
export CanonicalBundle, canonical_bundle


########################
# 5: Equality
########################

@doc Markdown.doc"""
    Base.:(==)(l1::ToricLineBundle, l2::ToricLineBundle)

Returns true if the toric line bundles `l1` and `l2` are isomorphic and false otherwise.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l1 = ToricLineBundle(P2, [0])
A toric line bundle on a normal toric variety

julia> l2 = l1^(-1)
A toric line bundle on a normal toric variety

julia> l1 == l2
true
```
"""
function Base.:(==)(l1::ToricLineBundle, l2::ToricLineBundle)
    if toric_variety(l1) !== toric_variety(l2)
        return false
    end
    return iszero(divisor_class(l1) - divisor_class(l2))
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
        push_attribute_if_exists!(properties_string, td, :isprincipal, "trivial")
        push_attribute_if_exists!(properties_string, td, :is_basepoint_free, "basepoint-free")
        ample_cb!(a,b) = push_attribute_if_exists!(a, b, :isample, "ample")
        push_attribute_if_exists!(properties_string, td, :is_very_ample, "very-ample"; callback=ample_cb!)
    end
    
    # print
    push!(properties_string, "line bundle on a normal toric variety")
    join(io, properties_string, ", ", " ")
    
end
