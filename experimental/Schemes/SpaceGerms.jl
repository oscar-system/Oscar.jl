export SpaceGerm, germ_at_point

export representative, point

export ambient_germ

import AbstractAlgebra: Ring

@Markdown.doc """
    AbsSpaceGerm{BaseRingType<:Ring,RingType<:Ring}

A space germ, i.e. a ringed space ``(X,O_{(X,x)})`` with local ring ``O_{(X,x)}`` of type `RingType` over a ring ``k`` of type `BaseRingType`.
"""
abstract type AbsSpaceGerm{BaseRingType<:Ring, RingType<:Ring} <: AbsSpec{BaseRingType, RingType} end


####################################################################
## two short hand definitions for internal use only
####################################################################
GermAtClosedPoint = Spec{<:Field, 
                         <:AbsLocalizedRing{<:Ring, <:Any, 
                                            <:MPolyComplementOfKPointIdeal}
                        }
GermAtGeometricPoint = Spec{<:Field, 
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
      return Spec(quotient_ring(OO(X)))
end

@attr Spec function representative(X::SpaceGerm{<:Ring, <:MPolyLocalizedRing})
    R = ambient_ring(X)
    return Spec(R)
end

function point(X::SpaceGerm{<:Any, <:Any, <:GermAtClosedPoint})
  return point_coordinates(inverted_set(OO(X)))
end

function point(X::SpaceGerm{<:Any, <:Any, <:GermAtGeometricPoint})
  return prime_ideal(inverted_set(OO(X)))
end

Oscar.ring(X::AbsSpaceGerm) = OO(X)

function Oscar.ideal(X::AbsSpaceGerm{<:Ring,<:MPolyQuoLocalizedRing})
    return modulus(OO(X))
end

function Oscar.ideal(X::AbsSpaceGerm{<:Ring,<:MPolyLocalizedRing})
    return ideal(OO(X),[zero(OO(X))])
end

@attr SpaceGerm function ambient_germ(X::AbsSpaceGerm{<:Ring,<:MPolyQuoLocalizedRing})
    Y,_ = germ_at_point(localized_ring(OO(X)))
    return Y
end

@attr SpaceGerm function ambient_germ(X::AbsSpaceGerm{<:Ring,<:MPolyLocalizedRing})
    return X
end

############################################################################################################
# allow user to specify point also as ideal
############################################################################################################

function _maxideal_to_point(I::MPolyIdeal)
  R=base_ring(I)
  o = degrevlex(gens(R))
  G=groebner_basis(I)
  LG = leading_ideal(I;ordering=o)
  dim(LG)==0 || error("Ideal does not describe finite set of points")
  vd,vbasis = _vdim_hack(LG)
  vd ==1 || error("Ideal does not describe a single K-point")
  return [leading_coefficient(normal_form(v,G)) for v in gens(R)]
end

### _vdim_hack only to be used until vdim of 0-dimensional ideals is implemented properly
### _vdim_hack assumes for simplicity that 
###         1. dim(I)=0 has been tested, 
###         2.I is a monomial ideal
function _vdim_hack(I::MPolyIdeal)
  leer=elem_type(base_ring(I))[]
  R = base_ring(I)
  if one(R) in I 
    return 0,leer
  end
  G=groebner_basis(I)
  M = ideal(R,gens(R))
  result=[R(1)]
  J = ideal(R,normal_form(gens(M),I))
  while dim(J) != length(gens(R))
    Jtemp = leer
    JN=gens(J)
    for i in 1:length(JN)
      if(JN[i]!=R(0))
        push!(result,JN[i])
        push!(Jtemp, JN[i])
      end
    end
    J = ideal(R,Jtemp)
    J = ideal(R,normal_form(gens(J*M),I))
  end
  return length(result),result
end

############################################################################################################
### constructors
############################################################################################################
function SpaceGerm(X::AbsSpec, a::Vector)
  R = ambient_ring(X)
  kk = coefficient_ring(R)
  b = [kk.(v) for v in a]  ## throws an error, if vector entries are not compatible
  U = MPolyComplementOfKPointIdeal(R,b)
  Y = Spec(Localization(OO(X), U)[1])
  return SpaceGerm(Y)
end

function SpaceGerm(X::AbsSpec, I::MPolyIdeal)
  R = base_ring(I)
  R == ambient_ring(X) || error("rings are not compatible")
  a = _maxideal_to_point(I)
  Y = SpaceGerm(X,a)
  return Y
end

function SpaceGerm(X::AbsSpec, I::MPolyQuoIdeal)
  A = base_ring(I)
  A == OO(X) || error("rings seem incompatible")
  R = base_ring(A)
  I = ideal(R, lift.(gens(I))) + modulus(A)
  a = _maxideal_to_point(I)
  Y = SpaceGerm(X,a)
  return Y
end

##### still missing: case of MPolyLocalizedIdeal and quo
#####               - for inverted set = powers of elements: point either outside of excluded hypersurface, in 
#####                 which case this should be referred to the global case, or nonsense
#####               - for complements of ideals: either I is the unique maximal ideal, 
#####                 i.e.the complement of the inverted set, or nonsense

function germ_at_point(X::AbsSpec, I::Ideal)
  Y = SpaceGerm(X, I)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

function germ_at_point(X::AbsSpec, a::Vector)
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

#########################################################################################
## for convenience of users thinking in terms of local rings
#########################################################################################

LocalRing = Union{MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                        <:MPolyComplementOfKPointIdeal},
                  MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                     <:MPolyComplementOfKPointIdeal},
                  MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                        <:MPolyComplementOfPrimeIdeal},
                  MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                     <:MPolyComplementOfPrimeIdeal}
                 }


function SpaceGerm(A::LocalRing)
  return SpaceGerm(Spec(A))
end

## and with identity map to keep usage consistent
function germ_at_point(A::LocalRing)
  X = SpaceGerm(A)
  restr_map = SpecMor(X, X, hom(OO(X), OO(X), gens(OO(X)), check=false), check=false)
  return X, restr_map
end

### basic functionality for space germs

##############################################################################
# note: ==, intersect are inherited from Spec
#       intersect with explicit fallback to Spec and change of return type
##############################################################################

function issubset(X::AbsSpaceGerm{<:Any, <:MPolyQuoLocalizedRing}, Y::AbsSpaceGerm{<:Any, <:MPolyQuoLocalizedRing})
  R = ambient_ring(X)
  R == ambient_ring(Y) || return false
  point(X) == point(Y) || return false
  IY=ideal(localized_ring(OO(X)),modulus(quotient_ring(OO(Y))))
  gens(IY)
  gens(modulus(quotient_ring(OO(X))))
  return issubset(IY,modulus(OO(X)))
end

function Base.intersect(X::AbsSpaceGerm, Y::AbsSpaceGerm)
  Z = intersect(underlying_scheme(X),underlying_scheme(Y))
  return SpaceGerm(Z)
end

function Base.union(X::AbsSpaceGerm, Y::AbsSpaceGerm)
  R = ambient_ring(X)
  R == ambient_ring(Y) || error("not subgerms of a common space germ")
  point(X) == point(Y) || error("not the same point of the germ")
  # comparison of points implicitly also checks that localization was performed at points
  # otherwise 'point' is not implemented
  I = intersect(modulus(quotient_ring(OO(X))),modulus(quotient_ring(OO(Y))))
  Z,_ = germ_at_point(MPolyQuoLocalizedRing(R, I ,inverted_set(OO(X))))
  return Z
end

##############################################################################
# note: singular_locus, is_smooth and is_regular are inherited from Spec
##############################################################################

function singular_locus(X::AbsSpaceGerm)
  return SpaceGerm(singular_locus(underlying_scheme(X)))
end