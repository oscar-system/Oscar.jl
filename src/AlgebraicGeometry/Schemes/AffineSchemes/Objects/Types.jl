########################################################################
#
# There are 4 types of affine schemes based on the following
#
# (1) (Polynomial) ring
# (2) Quotient of (polynoimal) ring
# (3) Localization of (polynomial) ring
# (4) Localization of quotient of (polynomial) ring
#
# Consequently, most method must be implemented 4 (or a multiple of 4) times.




###########################################
# (1) Abstract affine schemes
###########################################

#@doc raw"""
#    AbsAffineScheme{BaseRingType, RingType<:Ring}
#
#An affine scheme ``X = Spec(R)`` with ``R`` of type `RingType` over
#a ring ``𝕜`` of type `BaseRingType`.
#"""
#abstract type AbsAffineScheme{BaseRingType, RingType<:Ring} <: Scheme{BaseRingType} end
#
# Moved to src/forward_declarations.jl

############################################
# (2) Basic concrete type for affine schemes
############################################

@doc raw"""
    AffineScheme{BaseRingType, RingType}

An affine scheme ``X = Spec(R)`` with ``R`` a Noetherian ring of type `RingType`
over a base ring ``𝕜`` of type `BaseRingType`.
"""
@attributes mutable struct AffineScheme{BaseRingType, RingType<:Ring} <: AbsAffineScheme{BaseRingType, RingType}
  # the basic fields
  OO::RingType
  kk::BaseRingType

  function AffineScheme(OO::MPolyQuoLocRing)
    kk = coefficient_ring(base_ring(OO))
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end
  function AffineScheme(OO::MPolyLocRing)
    kk = coefficient_ring(base_ring(OO))
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end
  function AffineScheme(OO::MPolyRing)
    kk = coefficient_ring(OO)
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end
  function AffineScheme(OO::MPolyQuoRing)
    kk = coefficient_ring(base_ring(OO))
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end

  function AffineScheme(R::Ring)
    return new{typeof(ZZ), typeof(R)}(R, ZZ)
  end

  function AffineScheme(kk::Ring, R::Ring)
    return new{typeof(kk), typeof(R)}(R, kk)
  end

  function AffineScheme(kk::Field)
    return new{typeof(kk), typeof(kk)}(kk, kk)
  end
end



###########################################
# (3) Auxiliary types
###########################################

const StdAffineScheme = AbsAffineScheme{<:Ring, <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}

