################################################################################
#
#  Functionality for multivariate polynomials over tropical semirings
#  (alternatives for generic functions that do not work over tropical semirings)
#
################################################################################

# exponentiation has to be in multiple functions, otherwise calls become ambiguous
#   due to ^(a::AbstractAlgebra.Generic.MPoly{T}, b::Int64) where T<:RingElement
#   in AbstractAlgebra/src/generic/MPoly.jl:2257
function Base.:(^)(a::Generic.MPoly{<:TropicalSemiringElem}, n::Int)
    @req n>=0 "polynomial exponentiation with negative exponent"
    return Generic.pow_rmul(a,n)
end
function Base.:(^)(a::Generic.MPoly{<:TropicalSemiringElem}, n::ZZRingElem)
    return Generic.pow_rmul(a,Int(n))
end
