########################################################################
# Properties of AffineSchemeOpenSubschemeRingElem                                       #
########################################################################

########################################################################
# Required functionality for the ring interface                        #
########################################################################
is_domain_type(::Type{T}) where {T<:AffineSchemeOpenSubschemeRingElem} = true
is_domain_type(a::AffineSchemeOpenSubschemeRingElem) = is_domain_type(typeof(a))
is_exact_type(::Type{T}) where {T<:AffineSchemeOpenSubschemeRingElem} = true
is_exact_type(a::AffineSchemeOpenSubschemeRingElem) = is_exact_type(typeof(a))
is_domain_type(::Type{T}) where {T<:AffineSchemeOpenSubschemeRing} = true
is_domain_type(R::AffineSchemeOpenSubschemeRing) = is_domain_type(typeof(R))
is_exact_type(::Type{T}) where {T<:AffineSchemeOpenSubschemeRing} = true
is_exact_type(R::AffineSchemeOpenSubschemeRing) = is_exact_type(typeof(R))

