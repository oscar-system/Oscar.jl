export spec_open_ring_type, ring_type 
export scheme, restrictions, affine_patches

########################################################################
# Attributes of SpecOpenRing                                           #
########################################################################

########################################################################
# Type getters                                                         #
########################################################################
SpecOpenRing(U::SpecOpen) = SpecOpenRing(ambient(U), U)

spec_open_ring_type(::Type{T}) where {T<:Spec} = SpecOpenRing{T, open_subset_type(T)}
spec_open_ring_type(X::AbsSpec) = spec_open_ring_type(typeof(X))

ring_type(::Type{SpecOpenType}) where {SpecOpenType<:SpecOpen} = SpecOpenRing{affine_patch_type(SpecOpenType), SpecOpenType}
ring_type(U::SpecOpen) = ring_type(typeof(U))

########################################################################
# Basic attributes                                                     #
########################################################################
@Markdown.doc """
    scheme(R::SpecOpenRing)

The ring ``R = 𝒪(X, U)`` belongs to a sheaf of rings ``𝒪(X, -)`` and this returns 
the scheme ``X`` on which ``𝒪`` is defined.
"""
scheme(R::SpecOpenRing) = R.scheme
gens(R::SpecOpenRing) = R.(gens(ambient_ring(scheme(R))))

@Markdown.doc """
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

elem_type(R::SpecOpenRing) = elem_type(typeof(R))

parent_type(::Type{SpecOpenRingElem{S}}) where {S} = S
parent_type(f::SpecOpenRingElem) = parent_type(typeof(f))

########################################################################
# Basic getters                                                        #
########################################################################
parent(f::SpecOpenRingElem) = f.parent
scheme(f::SpecOpenRingElem) = scheme(parent(f))
domain(f::SpecOpenRingElem) = domain(parent(f))
restrictions(f::SpecOpenRingElem) = f.restrictions
affine_patches(f::SpecOpenRingElem) = affine_patches(domain(f))
npatches(f::SpecOpenRingElem) = length(restrictions(f))
getindex(f::SpecOpenRingElem, i::Int) = getindex(restrictions(f), i)
getindex(f::SpecOpenRingElem, U::AbsSpec) = restrictions(f)[domain(f)[U]]

