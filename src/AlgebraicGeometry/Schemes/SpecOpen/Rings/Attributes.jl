
########################################################################
# Attributes of SpecOpenRing                                           #
########################################################################

########################################################################
# Type getters                                                         #
########################################################################
SpecOpenRing(U::SpecOpen) = SpecOpenRing(ambient_scheme(U), U)

spec_open_ring_type(::Type{T}) where {T<:Spec} = SpecOpenRing{T, open_subset_type(T)}
spec_open_ring_type(X::AbsSpec) = spec_open_ring_type(typeof(X))

ring_type(::Type{SpecOpenType}) where {SpecOpenType<:SpecOpen} = SpecOpenRing{affine_patch_type(SpecOpenType), SpecOpenType}
ring_type(U::SpecOpen) = ring_type(typeof(U))

########################################################################
# Basic attributes                                                     #
########################################################################
@doc raw"""
    scheme(R::SpecOpenRing)

The ring ``R = 𝒪(X, U)`` belongs to a sheaf of rings ``𝒪(X, -)`` and this returns 
the scheme ``X`` on which ``𝒪`` is defined.
"""
scheme(R::SpecOpenRing) = R.scheme

gens(R::SpecOpenRing) = R.(gens(ambient_coordinate_ring(scheme(R))))
number_of_generators(R::SpecOpenRing) = number_of_generators(ambient_coordinate_ring(scheme(R)))
gen(R::SpecOpenRing, i::Int) = R(gen(ambient_coordinate_ring(scheme(R)), i))


@doc raw"""
    domain(R::SpecOpenRing)

For a ring ``R = 𝒪(X, U)``, return ``U``.
"""
domain(R::SpecOpenRing) = R.domain

########################################################################
# Attributes of SpecOpenRingElem                                       #
########################################################################

########################################################################
# Type getters                                                         #
########################################################################
elem_type(::Type{SpecOpenRing{S, T}}) where {S, T} = SpecOpenRingElem{SpecOpenRing{S, T}}
parent_type(::Type{SpecOpenRingElem{S}}) where {S} = S

########################################################################
# Basic getters                                                        #
########################################################################
parent(f::SpecOpenRingElem) = f.parent
scheme(f::SpecOpenRingElem) = scheme(parent(f))
domain(f::SpecOpenRingElem) = domain(parent(f))
restrictions(f::SpecOpenRingElem) = f.restrictions
affine_patches(f::SpecOpenRingElem) = affine_patches(domain(f))
number_of_patches(f::SpecOpenRingElem) = length(restrictions(f))
getindex(f::SpecOpenRingElem, i::Int) = getindex(restrictions(f), i)
getindex(f::SpecOpenRingElem, U::AbsSpec) = restrictions(f)[domain(f)[U]]

