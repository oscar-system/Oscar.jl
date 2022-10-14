export SpaceGerm

export representative, point

export germ_at_point, ring_of_germ, ideal_of_germ, ambient_germ

import AbstractAlgebra: Ring

@Markdown.doc """
    AbsSpaceGerm{BaseRingType<:Ring,RingType<:Ring}

A space germ, i.e. a ringed space ``(X,O_{(X,x)})`` with local ring ``O_{(X,x)}`` of type `RingType` over a ring ``k`` of type `BaseRingType`.
"""
abstract type AbsSpaceGerm{BaseRingType<:Ring, RingType<:Ring} <: AbsSpec{BaseRingType, RingType} end


####################################################################
## two short hand definitions for internal use only
####################################################################
GermAtClosedPoint = Spec{<:Ring, 
                         <:AbsLocalizedRing{<:Ring, <:Any, 
                                            <:MPolyComplementOfKPointIdeal}
                        }
GermAtGeometricPoint = Spec{<:Ring, 
                            <:AbsLocalizedRing{<:Ring, <:Any, 
                                               <:MPolyComplementOfPrimeIdeal}
                           }

###################################################################
## start of the definition of space germ functionality
###################################################################

@Markdown.doc """
    SpaceGerm{BaseRingType, RingType, SpecType}
A space germ ``(X,O_{(X,x)}``, i.e. a ringed space with underlying scheme ``X`` of type SpecType and local ring ``O_{(X,x)}`` of type `RingType` over some base ring ``k`` of type `BaseRingType`.
"""
@attributes mutable struct SpaceGerm{BaseRingType<:Ring, RingType<:Ring, SpecType<:Spec} <: AbsSpaceGerm{BaseRingType, RingType}
  X::SpecType

  function SpaceGerm(X::GermAtClosedPoint)
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X)
  end
  
  function SpaceGerm(X::GermAtGeometricPoint)
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X)
  end
end

### Getter functions

function underlying_scheme(X::SpaceGerm)
  return X.X
end

@attr Spec function representative(X::SpaceGerm{<:Ring, <:MPolyQuoLocalizedRing})
    R = base_ring(OO(X))
    I = modulus(OO(X))
    Q, _ = quo(R, I)
    return Spec(Q)
end

@attr Spec function representative(X::SpaceGerm{<:Ring, <:MPolyLocalizedRing})
    R = base_ring(OO(X))
    return Spec(R)
end

function point(X::SpaceGerm{<:Any, <:Any, <:GermAtClosedPoint})
  return point_coordinates(inverted_set(OO(X)))
end

function point(X::SpaceGerm{<:Any, <:Any, <:GermAtGeometricPoint})
  return prime_ideal(inverted_set(OO(X)))
end

ring_of_germ(X::AbsSpaceGerm) = OO(X)

function ideal_of_germ(X::AbsSpaceGerm{<:Ring,<:MPolyQuoLocalizedRing})
    return localized_ring(OO(X))(modulus(OO(X)))
end

function ideal_of_germ(X::AbsSpaceGerm{<:Ring,<:MPolyLocalizedRing})
    return ideal(OO(X),[zero(OO(X))])
end

function ambient_germ(X::AbsSpaceGerm{<:Ring,<:MPolyQuoLocalizedRing})
    Y,_ = germ_at_point(localized_ring(OO(X)))
    return Y
end

function ambient_germ(X::AbsSpaceGerm{<:Ring,<:MPolyLocalizedRing})
    return X
end

############################################################################################################
# allow user to specify point also as ideal
# currently broken -- suggestions welcome on how to do this Oscar-style
############################################################################################################
function _maxideal_to_point(I::MPolyIdeal)
  G = groebner_assure(I)
  dim(I)==0 || error("Ideal does not describe finite set of points")
  singular_assure(G)
  vdim(G)==1 || error("Ideal does not describe a single K-point")
  return [singular.reduce(v,G) for v in gens(base_ring(I))]
end

### constructors
function SpaceGerm(X::Spec, a::Vector)
  R = ambient_ring(X)
  kk = coefficient_ring(R)
  b = [kk.(v) for v in a]  ## throws and error, if vector entries are not compatible
  U = MPolyComplementOfKPointIdeal(R,b)
  Y = Spec(Localization(OO(X), U)[1])
  return SpaceGerm(Y)
end
  
function SpaceGerm(X::Spec, I::MPolyLocalizedIdeal)
  R = base_ring(base_ring(I))
  R == ambient_ring(X) || error("rings are not compatible")
  J = ideal(R, [numerator(p) for p in gens(I)]) # expected to be a single k-point, only tested in nex line
  a = _maxideal_to_point(J)
  Y = SpaceGerm(X,a)
  return Y
end

function SpaceGerm(X::Spec, I::MPolyQuoLocalizedIdeal)
  L = base_ring(I)
  R = base_ring(L)
  R == ambient_ring(X) || error("rings are not compatible")
  J = ideal(R, [numerator(p) for p in lift.(gens(I))]) + modulus(quotient_ring(L))
  a = _maxideal_to_point(J)
  Y = SpaceGerm(X,a)
  return Y
end  

function SpaceGerm(X::Spec, I::MPolyIdeal)
  R = base_ring(I)
  R == ambient_ring(X) || error("rings are not compatible")
  a = _maxideal_to_point(I)
  Y = SpaceGerm(X,a)
  return Y
end

function SpaceGerm(X::Spec, I::MPolyQuoIdeal)
  A = base_ring(I)
  A == OO(X) || error("rings are not compatible")
  R = base_ring(A)
  I = ideal(R, lift.(gens(I))) + modulus(A)
  a = _maxideal_to_point(I)
  Y = SpaceGerm(X,a)
  return Y
end

function germ_at_point(X::Spec, I::Ideal)
  Y = SpaceGerm(X, I)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

function germ_at_point(X::Spec, a::Vector)
  Y = SpaceGerm(X, a)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

function germ_at_point(A::MPolyRing, I::Ideal)
  X = Spec(A)
  Y = SpaceGerm(X, I)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

function germ_at_point(A::MPolyRing, a::Vector)
  X = Spec(A)
  Y = SpaceGerm(X, a)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

function germ_at_point(A::MPolyQuo, I::Ideal)
  X = Spec(A)
  Y = SpaceGerm(X, I)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

function germ_at_point(A::MPolyQuo, a::Vector)
  X = Spec(A)
  Y = SpaceGerm(X, a)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

LocalRing = Union{MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                        <:MPolyComplementOfKPointIdeal},
                  MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                     <:MPolyComplementOfKPointIdeal},
                  MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                        <:MPolyComplementOfPrimeIdeal},
                  MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                     <:MPolyComplementOfPrimeIdeal}
                 }

function germ_at_point(A::LocalRing; check::Bool=true)
  X = Spec(A)
  Y = SpaceGerm(X)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

### basic functionality for space germs

##############################################################################
# note: issubset, isempty, is_canonically_isomorphic, intersect are 
#       inherited from Spec
##############################################################################

function Base.union(X::AbsSpaceGerm, Y::AbsSpaceGerm)
  ambient_ring(X)==ambient_ring(Y) || error("not subgerms of a common space germ")
  point(X)==point(Y) || error("not the same point of the germ")
  # comparison of points implicitly also checks that localization was performed at points
  # otherwise 'point' is not implemented
  I=intersect(modulus(X),modulus(Y))
  Z,_ = germ_at_point(quo(localized_ring(OO(X)),I))
  return Z
end

##############################################################################
# note: singular_locus, is_smooth and is_regular are inherited from Spec
##############################################################################




