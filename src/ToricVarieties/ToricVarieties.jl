################################################################################
################################################################################
##  Structs
################################################################################
################################################################################
struct AffineNormalToricVariety
    polymakeNTV::Polymake.BigObject
end

struct NormalToricVariety
    polymakeNTV::Polymake.BigObject
end

const NormalToricVarieties = Union{AffineNormalToricVariety, NormalToricVariety}


################################################################################
################################################################################
##  Constructors
################################################################################
################################################################################
@doc Markdown.doc"""
    NormalToricVariety(PF::PolyhedralFan)

Construct the normal toric variety corresponding to a polyhedral fan `PF`.

# Example
Take `PF` to be the normal fan of the square.
```jldoctest
julia> square = cube(2)
A polyhedron in ambient dimension 2

julia> nf = normal_fan(square)
A polyhedral fan in ambient dimension 2

julia> ntv = NormalToricVariety(nf)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""    
function NormalToricVariety(PF::PolyhedralFan)
    pmntv = Polymake.fulton.NormalToricVariety(pm_fan(PF))
    return NormalToricVariety(pmntv)
end


@doc Markdown.doc"""
    NormalToricVariety(P::Polyhedron)

Construct the normal toric variety corresponding to the normal fan of the given
polyhedron `P`.

Note that this only coincides with the projective variety associated to `P`, if
`P` is normal.

# Example
Set `P` to be a square.
```jldoctest
julia> square = cube(2)
A polyhedron in ambient dimension 2

julia> ntv = NormalToricVariety(square)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""    
function NormalToricVariety(P::Polyhedron)
    fan = normal_fan(P)
    return NormalToricVariety(fan)
end


@doc Markdown.doc"""
    AffineNormalToricVariety(C::Cone)

Construct the affine normal toric variety corresponding to a polyhedral cone
`C`.

# Example
Set `C` to be the positive orthant in two dimensions.
```jldoctest
julia> C = positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function AffineNormalToricVariety(C::Cone)
    pmc = pm_cone(C)
    fan = Polymake.fan.check_fan_objects(pmc)
    pmntv = Polymake.fulton.NormalToricVariety(fan)
    return AffineNormalToricVariety(pmntv)
end


@doc Markdown.doc"""
    NormalToricVariety(C::Cone)

Construct the (affine) normal toric variety corresponding to a polyhedral cone
`C`.

# Example
Set `C` to be the positive orthant in two dimensions.
```jldoctest
julia> C = positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> ntv = NormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function NormalToricVariety(C::Cone)
    pmc = pm_cone(C)
    fan = Polymake.fan.check_fan_objects(pmc)
    pmntv = Polymake.fulton.NormalToricVariety(fan)
    return NormalToricVariety(pmntv)
end

function pm_ntv(ntv::NormalToricVarieties)
    return ntv.polymakeNTV
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, ntv::NormalToricVarieties)
    # fan = get_polyhedral_fan(ntv)
    pmntv = pm_ntv(ntv)
    ambdim = pmntv.FAN_AMBIENT_DIM
    print(io, "A normal toric variety corresponding to a polyhedral fan in ambient dimension $(ambdim)")
end



################################################################################
################################################################################
##  Access properties
################################################################################
################################################################################

@doc Markdown.doc"""
    toric_ideal_binomial_generators(antv::AffineNormalToricVariety)

Get the exponent vectors corresponding to the generators of the toric ideal
associated to the affine normal toric variety `antv`.

# Example
Take the cyclic quotient singularity corresponding to the pair of integers
`(2,5)`.
```jldoctest
julia> C = positive_hull([-2 5; 1 0])
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
    pmntv = antv.polymakeNTV
    result = pmntv.TORIC_IDEAL.BINOMIAL_GENERATORS
    return result
end

# function get_polyhedral_fan(ntv::NormalToricVarieties)
#     pmntv = pm_ntv(ntv)
#     fan = Polymake.fan.PolyhedralFan(pmntv)
#     println("get_polyhedral_fan done")
#     return PolyhedralFan(fan)
# end


@doc Markdown.doc"""
    isaffine(ntv::NormalToricVarieties)

Check whether a normal toric variety is affine, i.e. the associated polyhedral
fan has exactly one maximal cone.

# Examples
```jldoctest
julia> C = positive_hull([-2 5; 1 0])
A polyhedral cone in ambient dimension 2

julia> ntv = NormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> isaffine(ntv)
true

julia> square = cube(2)
A polyhedron in ambient dimension 2

julia> ntv = NormalToricVariety(square)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> isaffine(ntv)
false
```
"""
function isaffine(ntv::NormalToricVarieties)
    pmntv = pm_ntv(ntv)
    return pmntv.AFFINE::Bool
end


@doc Markdown.doc"""
    iscomplete(ntv::NormalToricVarieties)

Check whether a normal toric variety is complete, i.e. the associated
polyhedral fan covers the whole space.

# Examples
```jldoctest
julia> C = positive_hull([-2 5; 1 0])
A polyhedral cone in ambient dimension 2

julia> ntv = NormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> iscomplete(ntv)
false

julia> square = cube(2)
A polyhedron in ambient dimension 2

julia> ntv = NormalToricVariety(square)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> iscomplete(ntv)
true
```
"""
function iscomplete(ntv::NormalToricVarieties)
    pmntv = pm_ntv(ntv)
    return pmntv.COMPLETE::Bool
end


@doc Markdown.doc"""
    isnormal_variety(ntv::NormalToricVarieties)

Check whether a normal toric variety is normal, i.e. the semigroups generated
by the Hilbert basis of the dual cones to the maximal cones of the associated
polyhedral fan are all saturated.

This is trivially true for toric varieties coming from fans or cones.
"""
function isnormal_variety(ntv::NormalToricVarieties)
    return true
end


@doc Markdown.doc"""
    isprojective(ntv::NormalToricVarieties)

Check whether a normal toric variety is projective.
"""
function isprojective(ntv::NormalToricVarieties)
    pmntv = pm_ntv(ntv)
    return pmntv.PROJECTIVE::Bool
end


@doc Markdown.doc"""
    issimplicial(ntv::NormalToricVarieties)

Check whether a normal toric variety is simplicial.

This is trivially true for all two-dimensional normal toric varieties.

# Examples
Take the cone over a square at height one.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
A polyhedral cone in ambient dimension 3

julia> antv = AffineNormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 3

julia> issimplicial(antv)
false
```
"""
function issimplicial(ntv::NormalToricVarieties)
    pmntv = pm_ntv(ntv)
    return pmntv.SIMPLICIAL::Bool
end


@doc Markdown.doc"""
    issmooth(ntv::NormalToricVarieties)

Check whether a normal toric variety is smooth, i.e. whether all the maximal
cones of the associated polyhedral fan are simplicial and generated by a
lattice basis.

# Examples
```jldoctest
julia> C = positive_hull([-2 5; 1 0])
A polyhedral cone in ambient dimension 2

julia> ntv = NormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> issmooth(ntv)
false

julia> square = cube(2)
A polyhedron in ambient dimension 2

julia> ntv = NormalToricVariety(square)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> issmooth(ntv)
true
```
"""
function issmooth(ntv::NormalToricVarieties)
    pmntv = pm_ntv(ntv)
    return pmntv.SMOOTH::Bool
end



export
    AffineNormalToricVariety,
    NormalToricVariety,
#     get_polyhedral_fan,
    isaffine,
    iscomplete,
    isnormal_variety,
    isprojective,
    issimplicial,
    issmooth,
    toric_ideal_binomial_generators
