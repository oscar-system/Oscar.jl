#############################
# 1: The Julia type for toric cohomology classes
#############################

@attributes mutable struct CohomologyClass
    v::NormalToricVarietyType
    p::MPolyQuoRingElem
    function CohomologyClass(v::NormalToricVarietyType, p::MPolyQuoRingElem, completeness_check::Bool)
        @req parent(p) == cohomology_ring(v; completeness_check) "The polynomial must reside in the cohomology ring of the toric variety"
        return new(v, p)
    end
end


######################
# 2: Generic constructors
######################

@doc raw"""
    cohomology_class(v::NormalToricVarietyType, p::MPolyQuoRingElem; completeness_check::Bool = true)

Construct the toric cohomology class on the toric variety `v` corresponding to the polynomial `p`.  
The polynomial `p` must lie in the cohomology ring of `v`.

!!! note "Simplicial and complete toric varieties"
    This function assumes that the toric variety is both **simplicial** and **complete**.
    Since completeness checks can be slow, you may skip them by passing
    the optional keyword argument `completeness_check = false`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> c = cohomology_class(P2, gens(cohomology_ring(P2))[1])
Cohomology class on a normal toric variety given by x1
```
"""
cohomology_class(v::NormalToricVarietyType, p::MPolyQuoRingElem; completeness_check::Bool = true) = CohomologyClass(v, p, completeness_check)


@doc raw"""
    cohomology_class(d::ToricDivisor; completeness_check::Bool = true)

Construct the toric cohomology class corresponding to the toric divisor `d`.

!!! note "Simplicial and complete toric varieties"
    This function assumes that the underlying toric variety is both **simplicial** and **complete**.
    Since completeness checks can be slow, you may skip them by passing
    the optional keyword argument `completeness_check = false`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> d = toric_divisor(P2, [1, 2, 3])
Torus-invariant, non-prime divisor on a normal toric variety

julia> cohomology_class(d)
Cohomology class on a normal toric variety given by x1 + 2*x2 + 3*x3
```
"""
function cohomology_class(d::ToricDivisor; completeness_check::Bool = true)
  R = cohomology_ring(toric_variety(d); completeness_check)
  indets = gens(base_ring(R))
  coeff_ring = coefficient_ring(toric_variety(d))
  poly = R(sum(coeff_ring(coefficients(d)[k]) * indets[k] for k in 1:length(indets)))
  return CohomologyClass(toric_variety(d), poly, completeness_check)
end


@doc raw"""
    cohomology_class(c::ToricDivisorClass; completeness_check::Bool = true)

Construct the toric cohomology class corresponding to the toric divisor class `c`.

!!! note "Simplicial and complete toric varieties"
    This function assumes that the underlying toric variety is both **simplicial** and **complete**.
    Since completeness checks can be slow, you may skip them by passing
    the optional keyword argument `completeness_check = false`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> tdc = toric_divisor_class(P2, [2])
Divisor class on a normal toric variety

julia> cohomology_class(tdc)
Cohomology class on a normal toric variety given by 2*x1
```
"""
cohomology_class(c::ToricDivisorClass; completeness_check::Bool = true) = cohomology_class(toric_divisor(c), completeness_check = completeness_check)


@doc raw"""
    cohomology_class(l::ToricLineBundle; completeness_check::Bool = true)

Construct the toric cohomology class corresponding to the toric line bundle `l`.

!!! note "Simplicial and complete toric varieties"
    This function assumes that the underlying toric variety is both **simplicial** and **complete**.
    Since completeness checks can be slow, you may skip them by passing
    the optional keyword argument `completeness_check = false`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(P2, [2])
Toric line bundle on a normal toric variety

julia> polynomial(cohomology_class(l))
2*x1
```
"""
cohomology_class(l::ToricLineBundle; completeness_check::Bool = true) = cohomology_class(toric_divisor(l), completeness_check = completeness_check)


#################################
# 3: Addition, subtraction and scalar multiplication
#################################

function Base.:+(cc1::CohomologyClass, cc2::CohomologyClass)
    @req toric_variety(cc1) === toric_variety(cc2) "The cohomology classes must be defined on the same toric variety"
    ring = cohomology_ring(toric_variety(cc1))
    poly = polynomial(ring, cc1) + polynomial(ring, cc2)
    return CohomologyClass(toric_variety(cc1), poly, true)
end


function Base.:-(cc1::CohomologyClass, cc2::CohomologyClass)
    @req toric_variety(cc1) === toric_variety(cc2) "The cohomology classes must be defined on the same toric variety"
    ring = cohomology_ring(toric_variety(cc1))
    poly = polynomial(ring, cc1) - polynomial(ring, cc2)
    return CohomologyClass(toric_variety(cc1), poly, true)
end


Base.:*(c::QQFieldElem, cc::CohomologyClass) = CohomologyClass(toric_variety(cc), coefficient_ring(toric_variety(cc))(c) * polynomial(cc), true)
Base.:*(c::Rational{Int64}, cc::CohomologyClass) = CohomologyClass(toric_variety(cc), coefficient_ring(toric_variety(cc))(c) * polynomial(cc), true)
Base.:*(c::T, cc::CohomologyClass) where {T <: IntegerUnion} = CohomologyClass(toric_variety(cc), coefficient_ring(toric_variety(cc))(c) * polynomial(cc), true)


#################################
# 4: Wedge product
#################################

function Base.:*(cc1::CohomologyClass, cc2::CohomologyClass)
    @req toric_variety(cc1) === toric_variety(cc2) "The cohomology classes must be defined on the same toric variety"
    ring = cohomology_ring(toric_variety(cc1), completeness_check = false)
    poly = polynomial(ring, cc1) * polynomial(ring, cc2)
    return CohomologyClass(toric_variety(cc1), poly, true)
end


Base.:^(cc::CohomologyClass, p::T) where {T <: IntegerUnion} = CohomologyClass(toric_variety(cc), polynomial(cc)^p, true)


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
