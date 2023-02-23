#############################
# 1: The Julia type for toric cohomology classes
#############################

@attributes mutable struct CohomologyClass
    v::AbstractNormalToricVariety
    p::MPolyQuoRingElem
    function CohomologyClass(v::AbstractNormalToricVariety, p::MPolyQuoRingElem)
        if parent(p) != cohomology_ring(v)
            throw(ArgumentError("The polynomial must reside in the cohomology ring of the toric variety"))
        end
        return new(v, p)
    end
end
export CohomologyClass


######################
# 2: Generic constructors
######################

@doc Markdown.doc"""
    cohomology_class(v::AbstractNormalToricVariety, p::MPolyQuoRingElem)

Construct the toric cohomology class
on the toric variety `v` corresponding
to the polynomial `p`. Note that `p` must
reside in the cohomology ring of `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> c = cohomology_class(P2, gens(cohomology_ring(P2))[1])
A cohomology class on a normal toric variety given by x1
```
"""
cohomology_class(v::AbstractNormalToricVariety, p::MPolyQuoRingElem) = CohomologyClass(v, p)
export cohomology_class


@doc Markdown.doc"""
    cohomology_class(d::ToricDivisor)

Construct the toric cohomology class
corresponding to the toric divisor `d`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = toric_divisor(P2, [1, 2, 3])
A torus-invariant, non-prime divisor on a normal toric variety

julia> cohomology_class(d)
A cohomology class on a normal toric variety given by 6*x3
```
"""
function cohomology_class(d::ToricDivisor)
    indets = gens(cohomology_ring(toric_variety(d)))
    coeff_ring = coefficient_ring(toric_variety(d))
    poly = sum(coeff_ring(coefficients(d)[k]) * indets[k] for k in 1:length(indets))
    return CohomologyClass(toric_variety(d), poly)
end


@doc Markdown.doc"""
    cohomology_class(c::ToricDivisorClass)

Construct the toric cohomology class
corresponding to the toric divisor class `c`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, [2])
A divisor class on a normal toric variety

julia> cohomology_class(tdc)
A cohomology class on a normal toric variety given by 2*x3
```
"""
cohomology_class(c::ToricDivisorClass) = cohomology_class(toric_divisor(c))


@doc Markdown.doc"""
    cohomology_class(l::ToricLineBundle)

Construct the toric cohomology class
corresponding to the toric line bundle `l`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(P2, [2])
A toric line bundle on a normal toric variety

julia> polynomial(cohomology_class(l))
2*x3
```
"""
cohomology_class(l::ToricLineBundle) = cohomology_class(toric_divisor(l))


#################################
# 3: Addition, subtraction and scalar multiplication
#################################

function Base.:+(cc1::CohomologyClass, cc2::CohomologyClass)
    if toric_variety(cc1) !== toric_variety(cc2)
        throw(ArgumentError("The cohomology classes must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    ring = cohomology_ring(toric_variety(cc1))
    poly = polynomial(ring, cc1) + polynomial(ring, cc2)
    return CohomologyClass(toric_variety(cc1), poly)
end


function Base.:-(cc1::CohomologyClass, cc2::CohomologyClass)
    if toric_variety(cc1) !== toric_variety(cc2)
        throw(ArgumentError("The cohomology classes must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    ring = cohomology_ring(toric_variety(cc1))
    poly = polynomial(ring, cc1) - polynomial(ring, cc2)
    return CohomologyClass(toric_variety(cc1), poly)
end


Base.:*(c::fmpq, cc::CohomologyClass) = CohomologyClass(toric_variety(cc), coefficient_ring(toric_variety(cc))(c) * polynomial(cc))
Base.:*(c::Rational{Int64}, cc::CohomologyClass) = CohomologyClass(toric_variety(cc), coefficient_ring(toric_variety(cc))(c) * polynomial(cc))
Base.:*(c::T, cc::CohomologyClass) where {T <: IntegerUnion} = CohomologyClass(toric_variety(cc), coefficient_ring(toric_variety(cc))(c) * polynomial(cc))


#################################
# 4: Wedge product
#################################

function Base.:*(cc1::CohomologyClass, cc2::CohomologyClass)
    if toric_variety(cc1) !== toric_variety(cc2)
        throw(ArgumentError("The cohomology classes must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    ring = cohomology_ring(toric_variety(cc1))
    poly = polynomial(ring, cc1) * polynomial(ring, cc2)
    return CohomologyClass(toric_variety(cc1), poly)
end


Base.:^(cc::CohomologyClass, p::T) where {T <: IntegerUnion} = CohomologyClass(toric_variety(cc), polynomial(cc)^p)


########################
# 5: Equality
########################

function Base.:(==)(cc1::CohomologyClass, cc2::CohomologyClass)
    return toric_variety(cc1) === toric_variety(cc2) && iszero(polynomial(cc1-cc2))
end


######################
# 6: Display
######################s

function Base.show(io::IO, cc::CohomologyClass)
    join(io, "A cohomology class on a normal toric variety given by $(string(polynomial(cc)))")
end
