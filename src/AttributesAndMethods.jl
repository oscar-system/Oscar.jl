######################
# 1: Attributes of ToricVarieties
######################

@doc Markdown.doc"""
    dim( v::AbstractNormalToricVariety )

Computes the dimension of the normal toric variety `v`.
"""
function dim( v::AbstractNormalToricVariety )
    return v.polymakeNTV.FAN_DIM::Int
end
export dim


@doc Markdown.doc"""
    dim_of_torusfactor( v::AbstractNormalToricVariety )

Computes the dimension of the torus factor of the normal toric variety `v`.
"""
function dim_of_torusfactor( v::AbstractNormalToricVariety )

    if has_torusfactor( v ) == false
        return 0
    end
    
    dimension_of_fan = v.polymakeNTV.FAN_DIM::Int
    ambient_dimension = v.polymakeNTV.FAN_AMBIENT_DIM::Int
    return ambient_dimension - dimension_of_fan
end
export dim_of_torusfactor


@doc Markdown.doc"""
    euler_characteristic( v::AbstractNormalToricVariety )

Computes the Euler characteristic of the normal toric variety `v`.
"""
function euler_characteristic( v::AbstractNormalToricVariety )
    f_vector = Vector{Int}(v.polymakeNTV.F_VECTOR)
    return f_vector[ dim( v ) ]
end
export euler_characteristic


struct NefCone
           polymakeNefCone::Polymake.BigObject
end
export NefCone

"""
    nef_cone( v::NormalToricVariety )

Computes the nef cone of the normal toric variety `v`.
"""
function nef_cone( v::NormalToricVariety )
    return NefCone( v.polymakeNTV.NEF_CONE )
end
export nef_cone


struct MoriCone
           polymakeNefCone::Polymake.BigObject
end
export MoriCone

"""
    mori_cone( v::NormalToricVariety )

Computes the mori cone of the normal toric variety `v`.
"""
function mori_cone( v::NormalToricVariety )
    return MoriCone( v.polymakeNTV.MORI_CONE )
end
export mori_cone


######################
# 2: Methods of ToricVarieties
######################


@doc Markdown.doc"""
    ith_betti_number( v::AbstractNormalToricVariety, i::Int )

Compute the i-th Betti number of the normal toric variety `v`.
"""
function ith_betti_number( v::AbstractNormalToricVariety, i::Int )
    if isodd(i)
        return 0
    end
    k = i÷2
    f_vector = Vector{Int}(v.polymakeNTV.F_VECTOR)
    pushfirst!(f_vector, 1)
    betti_number = sum( (-1)^(i-k) * binomial(i,k) * f_vector[ dim( v ) - i + 1 ] for i=k:dim( v ))
    return betti_number
    
end
export ith_betti_number


@doc Markdown.doc"""
    toric_ideal_binomial_generators(antv::AffineNormalToricVariety)

Get the exponent vectors corresponding to the generators of the toric ideal
associated to the affine normal toric variety `antv`.

# Examples
Take the cyclic quotient singularity corresponding to the pair of integers
`(2,5)`.
```jldoctest
julia> C = Oscar.positive_hull([-2 5; 1 0])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> toric_ideal_binomial_generators(antv)
pm::Matrix<long>
-1 -1 2 1
-1 0 3 -1
0 -1 -1 2
```
"""
function toric_ideal_binomial_generators(antv::AffineNormalToricVariety)
    pmntv = pm_ntv(antv)
    result = pmntv.TORIC_IDEAL.BINOMIAL_GENERATORS
    return result
end
export toric_ideal_binomial_generators
