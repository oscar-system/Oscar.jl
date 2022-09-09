export SpaceGerm

export representative, point

export germ_at_point

import AbstractAlgebra: Ring

abstract type AbsSpaceGerm{BaseRingType<:Ring, RingType<:Ring} <: AbsSpec{BaseRingType, RingType} end

GermAtClosedPoint = Spec{<:Ring, 
                         <:AbsLocalizedRing{<:Ring, <:Any, 
                                            <:MPolyComplementOfKPointIdeal}
                        }
GermAtGeometricPoint = Spec{<:Ring, 
                            <:AbsLocalizedRing{<:Ring, <:Any, 
                                               <:MPolyComplementOfPrimeIdeal}
                           }

@attributes mutable struct SpaceGerm{BaseRingType<:Ring, RingType<:Ring, SpecType<:Spec} <: AbsSpaceGerm{BaseRingType, RingType}
  X::SpecType

  function SpaceGerm(X::GermAtClosedPoint)
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X)
  end
  
  function SpaceGerm(X::GermAtGeometricPoint)
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X)
  end
end

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

function SpaceGerm(X::Spec, I::MPolyLocalizedIdeal; check::Bool=true)
  R = base_ring(base_ring(I))
  R == ambient_ring(X) || error("rings are not compatible")
  P = saturated_ideal(I)
  if check
    is_prime(P) || error("the given ideal is not prime")
  end
  U = complement_of_ideal(P, check=check)
  Y = Spec(Localization(OO(X), U)[1])
  return SpaceGerm(Y)
end

function SpaceGerm(X::Spec, I::MPolyQuoLocalizedIdeal; check::Bool=true)
  R = base_ring(base_ring(I))
  R == ambient_ring(X) || error("rings are not compatible")
  P = saturated_ideal(I)
  if check
    is_prime(P) || error("the given ideal is not prime")
  end
  U = complement_of_ideal(P, check=check)
  Y = Spec(Localization(OO(X), U)[1])
  return SpaceGerm(Y)
end

function SpaceGerm(X::Spec, I::MPolyIdeal; check::Bool=true)
  if check
    is_prime(I) || error("the given ideal is not prime")
  end
  R = base_ring(I)
  R == ambient_ring(X) || error("rings are not compatible")

  U = complement_of_ideal(I, check=check)
  Y = Spec(Localization(OO(X), U)[1])

  return SpaceGerm(Y)
end

function SpaceGerm(X::Spec, I::MPolyQuoIdeal; check::Bool=true)
  A = base_ring(I)
  A == OO(X) || error("rings are not compatible")
  R = base_ring(A)
  P = ideal(R, lift.(gens(I))) + modulus(A)
  if check
    is_prime(P) || error("the given ideal is not prime")
  end
  U = complement_of_ideal(P, check=check)
  Y = Spec(Localization(OO(X), U)[1])

  return SpaceGerm(Y)
end

function germ_at_point(X::Spec, I::Ideal; check::Bool=true)
  Y = SpaceGerm(X, I, check=check)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

function germ_at_point(A::MPolyRing, I::Ideal; check::Bool=true)
  X = Spec(A)
  Y = SpaceGerm(X, I, check=check)
  restr_map = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

function germ_at_point(A::MPolyQuo, I::Ideal; check::Bool=true)
  X = Spec(A)
  Y = SpaceGerm(X, I, check=check)
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

function is_canonically_isomorphic(X::AbsSpec, Y::AbsSpec) 
  return is_canonically_isomorphic(underlying_scheme(X), underlying_scheme(Y))
end
