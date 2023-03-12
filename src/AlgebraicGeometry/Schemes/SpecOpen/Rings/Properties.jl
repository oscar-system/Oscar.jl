########################################################################
# Properties of SpecOpenRingElem                                       #
########################################################################

########################################################################
# Required functionality for the ring interface                        #
########################################################################
is_domain_type(::Type{T}) where {T<:SpecOpenRingElem} = true
is_domain_type(a::SpecOpenRingElem) = is_domain_type(typeof(a))
is_exact_type(::Type{T}) where {T<:SpecOpenRingElem} = true
is_exact_type(a::SpecOpenRingElem) = is_exact_type(typeof(a))
is_domain_type(::Type{T}) where {T<:SpecOpenRing} = true
is_domain_type(R::SpecOpenRing) = is_domain_type(typeof(R))
is_exact_type(::Type{T}) where {T<:SpecOpenRing} = true
is_exact_type(R::SpecOpenRing) = is_exact_type(typeof(R))

