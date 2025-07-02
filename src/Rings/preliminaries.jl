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

is_local(::Field) = true

is_known(::typeof(is_local), R::Ring) = false

