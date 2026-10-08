########################################################################
# Properties of AffineSchemeOpenSubschemeRingElem                                       #
########################################################################

########################################################################
# Required functionality for the ring interface                        #
########################################################################
is_domain_type(::Type{T}) where {T<:AffineSchemeOpenSubschemeRingElem} = true
is_exact_type(::Type{T}) where {T<:AffineSchemeOpenSubschemeRingElem} = true
