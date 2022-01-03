########################
# 1: The Julia type for toric line bundles
########################

abstract type ToricCoherentSheaf end

@attributes mutable struct ToricLineBundle <: ToricCoherentSheaf
    variety::AbstractNormalToricVariety
    divisor_class::GrpAbFinGenElem
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
    return ToricLineBundle(v, class, Dict())
end

@doc Markdown.doc"""
    ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{Int})

Convenience method for ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{fmpz}).

# Examples
```jldoctest
julia> v = toric_projective_space(2)
A normal, non-affine, smooth, projective, gorenstein, q-gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle( v, [ 2 ] )
A line bundle on a normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function ToricLineBundle(v::AbstractNormalToricVariety, input_class::Vector{Int})
    class = picard_group(v)(input_class)
    return ToricLineBundle(v, class, Dict())
end


########################
# 3: Display
########################

function Base.show(io::IO, line_bundle::ToricLineBundle)
    ambdim = pm_object(line_bundle.variety).FAN_AMBIENT_DIM
    print(io, "A line bundle on a normal toric variety corresponding to a polyhedral fan in ambient dimension $(ambdim)")
end
