#############################
# 1: The Julia type for toric cohomology classes
#############################

@attributes mutable struct CohomologyClass
    v::NormalToricVarietyType
    p::MPolyQuoRingElem
    function CohomologyClass(v::NormalToricVarietyType, p::MPolyQuoRingElem; quick::Bool = false)
        @req parent(p) == cohomology_ring(v, check = quick) "The polynomial must reside in the cohomology ring of the toric variety"
        return new(v, p)
    end
end


######################
# 2: Generic constructors
######################

@doc raw"""
    cohomology_class(v::NormalToricVarietyType, p::MPolyQuoRingElem)

Construct the toric cohomology class
on the toric variety `v` corresponding
to the polynomial `p`. Note that `p` must
reside in the cohomology ring of `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> c = cohomology_class(P2, gens(cohomology_ring(P2))[1])
Cohomology class on a normal toric variety given by x1
```
"""
cohomology_class(v::NormalToricVarietyType, p::MPolyQuoRingElem; quick::Bool = false) = CohomologyClass(v, p; quick = quick)


@doc raw"""
    cohomology_class(d::ToricDivisor)

Construct the toric cohomology class
corresponding to the toric divisor `d`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> d = toric_divisor(P2, [1, 2, 3])
Torus-invariant, non-prime divisor on a normal toric variety

julia> cohomology_class(d)
Cohomology class on a normal toric variety given by 6*x3
```
"""
function cohomology_class(d::ToricDivisor; quick::Bool = false)
  if quick == false
    indets = gens(cohomology_ring(toric_variety(d)))
    coeff_ring = coefficient_ring(toric_variety(d))
    poly = sum(coeff_ring(coefficients(d)[k]) * indets[k] for k in 1:length(indets))
    return CohomologyClass(toric_variety(d), poly)
  end
  indets = gens(parent(cohomology_ring(toric_variety(d), check = false)))
  coeff_ring = coefficient_ring(toric_variety(d))
  poly = cohomology_ring(toric_variety(d))(sum(coeff_ring(coefficients(d)[k]) * indets[k] for k in 1:length(indets)))
  return CohomologyClass(toric_variety(d), poly, quick = quick)
end


@doc raw"""
    cohomology_class(c::ToricDivisorClass)

Construct the toric cohomology class
corresponding to the toric divisor class `c`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> tdc = toric_divisor_class(P2, [2])
Divisor class on a normal toric variety

julia> cohomology_class(tdc)
Cohomology class on a normal toric variety given by 2*x3
```
"""
cohomology_class(c::ToricDivisorClass; quick::Bool = false) = cohomology_class(toric_divisor(c), quick = quick)


@doc raw"""
    cohomology_class(l::ToricLineBundle)

Construct the toric cohomology class
corresponding to the toric line bundle `l`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(P2, [2])
Toric line bundle on a normal toric variety

julia> polynomial(cohomology_class(l))
2*x3
```
"""
cohomology_class(l::ToricLineBundle; quick::Bool = false) = cohomology_class(toric_divisor(l), quick = quick)


#################################
# 3: Addition, subtraction and scalar multiplication
#################################

function Base.:+(cc1::CohomologyClass, cc2::CohomologyClass)
    @req toric_variety(cc1) === toric_variety(cc2) "The cohomology classes must be defined on the same toric variety"
    ring = cohomology_ring(toric_variety(cc1))
    poly = polynomial(ring, cc1) + polynomial(ring, cc2)
    return CohomologyClass(toric_variety(cc1), poly)
end


function Base.:-(cc1::CohomologyClass, cc2::CohomologyClass)
    @req toric_variety(cc1) === toric_variety(cc2) "The cohomology classes must be defined on the same toric variety"
    ring = cohomology_ring(toric_variety(cc1))
    poly = polynomial(ring, cc1) - polynomial(ring, cc2)
    return CohomologyClass(toric_variety(cc1), poly)
end


Base.:*(c::QQFieldElem, cc::CohomologyClass) = CohomologyClass(toric_variety(cc), coefficient_ring(toric_variety(cc))(c) * polynomial(cc))
Base.:*(c::Rational{Int64}, cc::CohomologyClass) = CohomologyClass(toric_variety(cc), coefficient_ring(toric_variety(cc))(c) * polynomial(cc))
Base.:*(c::T, cc::CohomologyClass) where {T <: IntegerUnion} = CohomologyClass(toric_variety(cc), coefficient_ring(toric_variety(cc))(c) * polynomial(cc))


#################################
# 4: Wedge product
#################################

function Base.:*(cc1::CohomologyClass, cc2::CohomologyClass)
    @req toric_variety(cc1) === toric_variety(cc2) "The cohomology classes must be defined on the same toric variety"
    ring = cohomology_ring(toric_variety(cc1))
    poly = polynomial(ring, cc1) * polynomial(ring, cc2)
    return CohomologyClass(toric_variety(cc1), poly)
end


Base.:^(cc::CohomologyClass, p::T) where {T <: IntegerUnion} = CohomologyClass(toric_variety(cc), polynomial(cc)^p)


########################
# 5: Equality and hash
########################

function Base.:(==)(cc1::CohomologyClass, cc2::CohomologyClass) 
    toric_variety(cc1) === toric_variety(cc2) && polynomial(cc1) == polynomial(cc2)
end

function Base.hash(cc::CohomologyClass, h::UInt) 
    b = 0x4de32042e67d89c8 % UInt
    h = hash(toric_variety(cc), h)
    h = hash(polynomial(cc), h)
    return xor(h, b)
end

######################
# 6: Display
######################s

Base.show(io::IO, cc::CohomologyClass) = join(io, "Cohomology class on a normal toric variety given by $(string(polynomial(cc)))")
