########################
# 1: The Julia type for toric line bundles
########################

abstract type ToricCoherentSheaf end

@attributes mutable struct ToricLineBundle <: ToricCoherentSheaf
    variety::AbstractNormalToricVariety
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
julia> v = toric_projective_space(2)
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
# 3: Display
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
