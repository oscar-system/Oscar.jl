
########################################################################
# Attributes of AffineSchemeOpenSubschemeRing                                           #
########################################################################

########################################################################
# Type getters                                                         #
########################################################################
AffineSchemeOpenSubschemeRing(U::AffineSchemeOpenSubscheme) = AffineSchemeOpenSubschemeRing(ambient_scheme(U), U)

affine_scheme_open_subscheme_ring_type(::Type{T}) where {T<:AffineScheme} = AffineSchemeOpenSubschemeRing{T, open_subset_type(T)}
affine_scheme_open_subscheme_ring_type(X::AbsAffineScheme) = affine_scheme_open_subscheme_ring_type(typeof(X))

ring_type(::Type{AffineSchemeOpenSubschemeType}) where {AffineSchemeOpenSubschemeType<:AffineSchemeOpenSubscheme} = AffineSchemeOpenSubschemeRing{affine_patch_type(AffineSchemeOpenSubschemeType), AffineSchemeOpenSubschemeType}
ring_type(U::AffineSchemeOpenSubscheme) = ring_type(typeof(U))

########################################################################
# Basic attributes                                                     #
########################################################################
@doc raw"""
    scheme(R::AffineSchemeOpenSubschemeRing)

The ring ``R = 𝒪(X, U)`` belongs to a sheaf of rings ``𝒪(X, -)`` and this returns
the scheme ``X`` on which ``𝒪`` is defined.
"""
scheme(R::AffineSchemeOpenSubschemeRing) = R.scheme

gens(R::AffineSchemeOpenSubschemeRing) = R.(gens(ambient_coordinate_ring(scheme(R))))
number_of_generators(R::AffineSchemeOpenSubschemeRing) = number_of_generators(ambient_coordinate_ring(scheme(R)))
gen(R::AffineSchemeOpenSubschemeRing, i::Int) = R(gen(ambient_coordinate_ring(scheme(R)), i))


@doc raw"""
    domain(R::AffineSchemeOpenSubschemeRing)

For a ring ``R = 𝒪(X, U)``, return ``U``.
"""
domain(R::AffineSchemeOpenSubschemeRing) = R.domain

########################################################################
# Attributes of AffineSchemeOpenSubschemeRingElem                                       #
########################################################################

########################################################################
# Type getters                                                         #
########################################################################
elem_type(::Type{AffineSchemeOpenSubschemeRing{S, T}}) where {S, T} = AffineSchemeOpenSubschemeRingElem{AffineSchemeOpenSubschemeRing{S, T}}
parent_type(::Type{AffineSchemeOpenSubschemeRingElem{S}}) where {S} = S

########################################################################
# Basic getters                                                        #
########################################################################
parent(f::AffineSchemeOpenSubschemeRingElem) = f.parent
scheme(f::AffineSchemeOpenSubschemeRingElem) = scheme(parent(f))
domain(f::AffineSchemeOpenSubschemeRingElem) = domain(parent(f))
restrictions(f::AffineSchemeOpenSubschemeRingElem) = f.restrictions
affine_patches(f::AffineSchemeOpenSubschemeRingElem) = affine_patches(domain(f))
number_of_patches(f::AffineSchemeOpenSubschemeRingElem) = length(restrictions(f))
getindex(f::AffineSchemeOpenSubschemeRingElem, i::Int) = getindex(restrictions(f), i)
getindex(f::AffineSchemeOpenSubschemeRingElem, U::AbsAffineScheme) = restrictions(f)[domain(f)[U]]

