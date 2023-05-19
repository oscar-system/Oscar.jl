export SpaceGerm
export ambient_germ
export germ_at_point
export point
export rational_point_coordinates
export representative

import AbstractAlgebra: Ring

@doc raw"""
    AbsSpaceGerm{BaseRingType<:Ring,RingType<:Ring}

A space germ, i.e. a ringed space ``(X,O_{(X,x)})`` with local ring ``O_{(X,x)}`` of type `RingType` over a coefficient field ``k`` of type `BaseRingType`.
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

@doc raw"""
    SpaceGerm{BaseRingType, RingType, SpecType}
A space germ ``(X,O_{(X,x)}``, i.e. a ringed space with underlying scheme ``X`` of type SpecType and local ring ``O_{(X,x)}`` of type `RingType` over some base ring ``k`` of type `BaseRingType`.
"""
@attributes mutable struct SpaceGerm{BaseRingType<:Ring, RingType<:Ring, SpecType<:Spec} <: AbsSpaceGerm{BaseRingType, RingType}
  X::SpecType

  function SpaceGerm(X::GermAtClosedPoint)
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X)
  end
  
## the following one is currently unused..
  function SpaceGerm(X::GermAtGeometricPoint)
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X)
  end
end

### Getter functions

function underlying_scheme(X::SpaceGerm)
  return X.X
end

@attr Spec function representative(X::SpaceGerm{<:Ring, <:MPolyQuoLocRing})
      return Spec(underlying_quotient(OO(X)))
end

@attr Spec function representative(X::SpaceGerm{<:Ring, <:MPolyLocRing})
    R = ambient_coordinate_ring(X)
    return Spec(R)
end

function point(X::SpaceGerm{<:Any, <:Any, <:GermAtClosedPoint})
  return point_coordinates(inverted_set(OO(X)))
end

## currently unused case
function point(X::SpaceGerm{<:Any, <:Any, <:GermAtGeometricPoint})
  return prime_ideal(inverted_set(OO(X)))
end

Oscar.ring(X::AbsSpaceGerm) = OO(X)

function Oscar.ideal(X::AbsSpaceGerm{<:Ring,<:MPolyQuoLocRing})
    return modulus(OO(X))
end

function Oscar.ideal(X::AbsSpaceGerm{<:Ring,<:MPolyLocRing})
    return ideal(OO(X),[zero(OO(X))])
end

@attr SpaceGerm function ambient_germ(X::AbsSpaceGerm{<:Ring,<:MPolyQuoLocRing})
    Y,_ = germ_at_point(localized_ring(OO(X)))
    return Y
end

@attr SpaceGerm function ambient_germ(X::AbsSpaceGerm{<:Ring,<:MPolyLocRing})
    return X
end

############################################################################################################
# allow user to specify point also as ideal
############################################################################################################

@doc raw"""
    rational_point_coordinates(I::MPolyIdeal)

Return the $k$-coordinates of the point corresponding to a maximal ideal
$I \in k[x_1,\dots,x_n]$, which describes a $k$-point. If $I$ is not maximal
or does not describe a point with coordinates in the field $k$, an error 
exception results.

# Examples
```jldoctest
julia> R, (x, y) = QQ["x","y"]
(Multivariate Polynomial Ring in x, y over Rational Field, QQMPolyRingElem[x, y])

julia> I = ideal(R, [x-1,y-3])
ideal(x - 1, y - 3)

julia> rational_point_coordinates(I)
2-element Vector{QQFieldElem}:
 1
 3

```
"""
function rational_point_coordinates(I::MPolyIdeal)
  R=base_ring(I)
  o = degrevlex(gens(R))
  G=groebner_basis(I)
  LG = leading_ideal(I;ordering=o)
  dim(LG)==0 || error("Ideal does not describe finite set of points")
  vd,vbasis = _vdim_hack(LG)
  vd ==1 || error("Ideal does not describe a single K-point")
  return [AbstractAlgebra.leading_coefficient(normal_form(v,I)) for v in gens(R)] # TODO does the ordering matter?
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
  M = ideal(R,gens(R))
  result=[R(1)]
  J = ideal(R,normal_form(gens(M),I))
  while dim(J) != ngens(R)
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
  R = ambient_coordinate_ring(X)
  kk = coefficient_ring(R)
  b = [kk.(v) for v in a]  ## throws an error, if vector entries are not compatible
  U = MPolyComplementOfKPointIdeal(R,b)
  Y = Spec(Localization(OO(X), U)[1])
  Z = SpaceGerm(Y)
  set_attribute!(Z,:representative,X)
  return SpaceGerm(Y)
end

function SpaceGerm(X::AbsSpec, I::MPolyIdeal)
  R = base_ring(I)
  R === ambient_coordinate_ring(X) || error("rings are not compatible")
  a = rational_point_coordinates(I)
  Y = SpaceGerm(X,a)
  return Y
end

to_poly_ideal(I::MPolyQuoIdeal) = ideal(base_ring(base_ring(I)),lift.(gens(I))) + modulus(base_ring(I))
to_poly_ideal(I::MPolyLocalizedIdeal) = ideal(base_ring(base_ring(I)), gens(saturated_ideal(I)))
to_poly_ideal(I::MPolyQuoLocalizedIdeal) = ideal(base_ring(base_ring(I)),gens(saturated_ideal(I))) + modulus(underlying_quotient(base_ring(I)))

function SpaceGerm(X::AbsSpec, I::Ideal)
  A = base_ring(I)
  A === OO(X) || error("rings are incompatible")
  J = to_poly_ideal(I)
  a = rational_point_coordinates(J)
  Y = SpaceGerm(X,a)
  return Y
end

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

function germ_at_point(A::MPolyQuoRing, I::Ideal)
  X = Spec(A)
  Y = SpaceGerm(X, I)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

function germ_at_point(A::MPolyQuoRing, a::Vector)
  X = Spec(A)
  Y = SpaceGerm(X, a)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

#########################################################################################
## for convenience of users thinking in terms of local rings
#########################################################################################

LocalRing = Union{MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, 
                                        <:MPolyComplementOfKPointIdeal},
                  MPolyLocRing{<:Any, <:Any, <:Any, <:Any, 
                                     <:MPolyComplementOfKPointIdeal},
                  MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, 
                                        <:MPolyComplementOfPrimeIdeal},
                  MPolyLocRing{<:Any, <:Any, <:Any, <:Any, 
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

function issubset(X::AbsSpaceGerm{<:Any, <:MPolyQuoLocRing}, Y::AbsSpaceGerm{<:Any, <:MPolyQuoLocRing})
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  point(X) == point(Y) || return false
  IY=ideal(localized_ring(OO(X)), gens(modulus(underlying_quotient(OO(Y)))))
  return issubset(IY,modulus(OO(X)))
end

function issubset(X::AbsSpaceGerm{<:Any, <:MPolyLocRing}, Y::AbsSpaceGerm{<:Any, <:MPolyQuoLocRing})
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  point(X) == point(Y) || return false
  return iszero(modulus(OO(Y)))
end

function issubset(X::AbsSpaceGerm{<:Any, <:MPolyLocRing}, Y::AbsSpaceGerm{<:Any, <:MPolyLocRing})
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  return point(X) == point(Y)
end

function issubset(X::AbsSpaceGerm{<:Any, <:MPolyQuoLocRing}, Y::AbsSpaceGerm{<:Any, <:MPolyLocRing})
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  return point(X) == point(Y)
end

function Base.intersect(X::AbsSpaceGerm, Y::AbsSpaceGerm)
  point(X) == point(Y) || error("not the same point of the germ")
  Z = intersect(underlying_scheme(X),underlying_scheme(Y))
  return SpaceGerm(Z)
end

function Base.union(X::AbsSpaceGerm, Y::AbsSpaceGerm)
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("not subgerms of a common space germ")
  point(X) == point(Y) || error("not the same point of the germ")
  # comparison of points implicitly also checks that localization was performed at points
  # otherwise 'point' is not implemented
  I = intersect(modulus(underlying_quotient(OO(X))),modulus(underlying_quotient(OO(Y))))
  Z,_ = germ_at_point(MPolyQuoLocRing(R, I ,inverted_set(OO(X))))
  return Z
end

##############################################################################
# note: singular_locus, is_smooth and is_regular are inherited from Spec
##############################################################################

# We want the singular locus of a `SpaceGerm` to be a `SpaceGerm` again and 
# not a plain `Spec`.
function singular_locus(X::AbsSpaceGerm)
  S, inc = singular_locus(underlying_scheme(X))
  Sgerm = SpaceGerm(S)
  return Sgerm, ClosedEmbedding(SpecMor(Sgerm, X, pullback(inc), check=false), image_ideal(inc), check=false)
end

function subscheme(X::SpaceGerm, I::Ideal)
  base_ring(I) === OO(X) || error("ideal does not belong to the correct ring")
  Y = subscheme(underlying_scheme(X), I)
  return SpaceGerm(Y)
end

