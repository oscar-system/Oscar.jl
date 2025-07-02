########################################################################
# This file contains declarations of functions and types which need 
# to come before the later inclusions in Rings.jl
########################################################################

@doc raw"""
    is_local(R::Ring)

Return whether a given ring `R` is a local ring.
"""
function is_local(R::Ring)
  error("check whether $R is local is not implemented")
end

# In general we can not assume it to be known whether a given ring is local
is_known(::typeof(is_local), R::Ring) = false

is_local(::Field) = true

is_local(R::MPolyRing{<:FieldElem}) = is_zero(ngens(R))
is_known(::typeof(is_local), R::MPolyRing{<:FieldElem}) = true

