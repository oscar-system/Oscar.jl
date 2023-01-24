########################################################################
# Computation of kernels of homomorphisms of polynomial rings 
# over coefficient rings which are not fields.
#
# In the admissible cases, this flattens the structure of the rings 
# and computes the kernel on the flat counterpart. 
########################################################################
function kernel(
    f::MPolyAnyMap{<:AbstractAlgebra.Generic.MPolyRing{<:RingElemType}, 
                   <:AbstractAlgebra.Generic.MPolyRing{<:RingElemType}
                  }
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                 <:MPolyQuoLocalizedRingElem
                                }
          }
  f_flat = _flatten(f)
  K_flat = kernel(f_flat)
  K = ideal(domain(f), _flat_to_orig(domain(f)).(gens(K_flat)))
  #set_attribute!(K, :_flat_counterpart, K_flat) # TODO: Support attribute storage?
  return K
end

@attr function _flatten(
    S::AbstractAlgebra.Generic.MPolyRing{RingElemType}
  ) where {RingElemType <: MPolyElem}
  R = base_ring(S)
  kk = coefficient_ring(R)
  S_flat, _ = PolynomialRing(kk, vcat(symbols(S), symbols(R)))
  R_to_S_flat = hom(R, S_flat, gens(S_flat)[ngens(S)+1:end])
  S_to_S_flat = hom(S, S_flat, R_to_S_flat, gens(S_flat)[1:ngens(S)])
  S_flat_to_S = hom(S_flat, S, vcat(gens(S), S.(gens(R))))
  
  set_attribute!(S, :_flat_to_orig, S_flat_to_S)
  set_attribute!(S, :_orig_to_flat, S_to_S_flat)
  set_attribute!(S, :_base_ring_to_flat, R_to_S_flat)

  return S_flat
end

@attr function _flatten(
    S::AbstractAlgebra.Generic.MPolyRing{RingElemType}
  ) where {RingElemType <: MPolyQuoElem}
  A = base_ring(S)::MPolyQuo
  R = base_ring(A)::MPolyRing
  I = modulus(A)::MPolyIdeal
  kk = coefficient_ring(R)::Field

  # Before building S_flat, we have to create a polynomial 
  # ring T_flat from which we can pass to the quotient.
  T_flat, _ = PolynomialRing(kk, vcat(symbols(S), symbols(R)))
  R_to_T_flat = hom(R, T_flat, gens(T_flat)[ngens(S)+1:end])
  
  I_flat = ideal(T_flat, R_to_T_flat.(gens(I)))
  S_flat, pr = quo(T_flat, I_flat)

  A_to_S_flat = hom(A, S_flat, gens(S_flat)[ngens(S)+1:end])
  
  S_to_S_flat = hom(S, S_flat, A_to_S_flat, gens(S_flat)[1:ngens(S)])
  S_flat_to_S = hom(S_flat, S, vcat(gens(S), S.(gens(A))))

  set_attribute!(S, :_flat_to_orig, S_flat_to_S)
  set_attribute!(S, :_orig_to_flat, S_to_S_flat)
  set_attribute!(S, :_base_ring_to_flat, A_to_S_flat)

  return S_flat
end

@attr function _flatten(
    S::AbstractAlgebra.Generic.MPolyRing{RingElemType}
  ) where {RingElemType <: MPolyQuoLocalizedRingElem}
  Q = base_ring(S)::MPolyQuoLocalizedRing
  R = base_ring(Q)::MPolyRing
  L = localized_ring(Q)::MPolyLocalizedRing
  A = underlying_quotient(Q)::MPolyQuo
  U = inverted_set(Q)::AbsMultSet
  I = modulus(Q)::MPolyLocalizedIdeal
  kk = coefficient_ring(R)::Field

  # Before building S_flat, we have to create a polynomial 
  # ring T_flat from which we can pass to the localized quotient.
  T_flat, _ = PolynomialRing(kk, vcat(symbols(S), symbols(R)))
  R_to_T_flat = hom(R, T_flat, gens(T_flat)[ngens(S)+1:end])

  T_flat_loc, _ = localization(T_flat, R_to_T_flat(U)) # Will throw if multiplicative set can not be transfered. 
  L_to_T_flat_loc = hom(L, T_flat_loc, gens(T_flat_loc)[ngens(S)+1:end])
  
  I_flat = ideal(T_flat_loc, L_to_T_flat_loc.(gens(I)))
  S_flat, pr = quo(T_flat_loc, I_flat)

  Q_to_S_flat = hom(Q, S_flat, gens(S_flat)[ngens(S)+1:end])
  
  S_to_S_flat = hom(S, S_flat, Q_to_S_flat, gens(S_flat)[1:ngens(S)])
  S_flat_to_S = hom(S_flat, S, vcat(gens(S), S.(gens(Q))))

  set_attribute!(S, :_flat_to_orig, S_flat_to_S)
  set_attribute!(S, :_orig_to_flat, S_to_S_flat)
  set_attribute!(S, :_base_ring_to_flat, Q_to_S_flat)

  return S_flat
end

@attr function _flatten(
    S::AbstractAlgebra.Generic.MPolyRing{RingElemType}
  ) where {RingElemType <: MPolyLocalizedRingElem}
  L = base_ring(S)::MPolyLocalizedRing
  R = base_ring(L)::MPolyRing
  U = inverted_set(L)::AbsMultSet
  kk = coefficient_ring(R)::Field

  # Before building S_flat, we have to create a polynomial 
  # ring T_flat from which we can pass to the localized quotient.
  T_flat, _ = PolynomialRing(kk, vcat(symbols(S), symbols(R)))
  R_to_T_flat = hom(R, T_flat, gens(T_flat)[ngens(S)+1:end])

  S_flat, _ = localization(T_flat, R_to_T_flat(U)) # Will throw if multiplicative set can not be transfered. 
  L_to_S_flat = hom(L, S_flat, gens(S_flat)[ngens(S)+1:end])
  
  R_to_S_flat = hom(R, S_flat, gens(S_flat)[ngens(S)+1:end])
  
  S_to_S_flat = hom(S, S_flat, L_to_S_flat, gens(S_flat)[1:ngens(S)])
  S_flat_to_S = hom(S_flat, S, vcat(gens(S), S.(gens(R))))

  set_attribute!(S, :_flat_to_orig, S_flat_to_S)
  set_attribute!(S, :_orig_to_flat, S_to_S_flat)
  set_attribute!(S, :_base_ring_to_flat, L_to_S_flat)

  return S_flat
end

### Transfer of multiplicative sets along ring homomorphisms 
function (phi::MPolyAnyMap{<:MPolyRing, <:MPolyRing, Nothing})(U::MPolyPowersOfElement;
                                                               check::Bool=true
                                                              )
  ambient_ring(U) === domain(phi) || error("multiplicative set does not lay in the domain of the morphism")
  S = codomain(phi) 
  SU = MPolyPowersOfElement(S, phi.(denominators(U)))
  return SU
end

function (phi::MPolyAnyMap{<:MPolyRing, <:MPolyRing, Nothing})(U::MPolyComplementOfPrimeIdeal;
                                                               check::Bool=true
                                                              )
  ambient_ring(U) === domain(phi) || error("multiplicative set does not lay in the domain of the morphism")
  S = codomain(phi) 
  Q = ideal(S, phi.(gens(prime_ideal(U))))
  SU = MPolyComplementOfPrimeIdeal(S, Q, check=check)
  return SU
end

### Canonical maps involved in flattening of towers of polynomial rings
@attr function _flat_to_orig(
    S::AbstractAlgebra.Generic.MPolyRing{RingElemType}
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                 <:MPolyQuoLocalizedRingElem
                                }
          }
  _flatten(S) # Will set the attribute
  return get_attribute(S, :_flat_to_orig)
end

@attr function _orig_to_flat(
    S::AbstractAlgebra.Generic.MPolyRing{RingElemType}
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                 <:MPolyQuoLocalizedRingElem
                                }
          }
  _flatten(S) # Will set the attribute
  return get_attribute(S, :_orig_to_flat)
end

@attr function _base_ring_to_flat(
    S::AbstractAlgebra.Generic.MPolyRing{RingElemType}
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                 <:MPolyQuoLocalizedRingElem
                                }
          }
  _flatten(S) # Will set the attribute
  return get_attribute(S, :_base_ring_to_flat)
end

### Computation of induced morphisms on flattened towers of polynomial rings
@attr function _flatten(
    f::MPolyAnyMap{AbstractAlgebra.Generic.MPolyRing{RingElemType}, 
                   AbstractAlgebra.Generic.MPolyRing{RingElemType},
                   Nothing
                  }
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                 <:MPolyQuoLocalizedRingElem
                                }
          }
  S = domain(f)
  T = codomain(f)
  RS = base_ring(domain(f))

  S_flat = _flatten(S)
  T_flat = _flatten(T)
  imgs = _orig_to_flat(T).(f.(gens(S))) # The first half of the images
  imgs = vcat(imgs, _base_ring_to_flat(T).(gens(RS))) # the ones from the base ring
  return hom(S_flat, T_flat, imgs)
end

@attr function _flatten(
    f::MPolyAnyMap{AbstractAlgebra.Generic.MPolyRing{RingElemType}, 
                   AbstractAlgebra.Generic.MPolyRing{RingElemType}
                   # Note the missing requirement here: It allows for a non-trivial coefficient map
                  }
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                 <:MPolyQuoLocalizedRingElem
                                }
          }
  S = domain(f)
  T = codomain(f)
  RS = base_ring(domain(f))

  S_flat = _flatten(S)
  T_flat = _flatten(T)
  imgs = _orig_to_flat(T).(f.(gens(S))) # The first half of the images
  imgs = vcat(imgs, _orig_to_flat(T).(T.(coefficient_map(f).(gens(RS))))) # the ones from the base ring
  return hom(S_flat, T_flat, imgs)
end

### Deflecting some of the basic methods for ideals in towers of polynomial rings
@attr function _flat_counterpart(
    I::MPolyIdeal{RingElemType}
  ) where { RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                  <:MPolyQuoLocalizedRingElem
                                 }
          }
  S = base_ring(I)
  S_flat = _flatten(S)
  return ideal(S_flat, _orig_to_flat(S).(gens(I)))
end

function ideal_membership(
    x::AbstractAlgebra.Generic.MPoly{RingElemType},
    I::MPolyIdeal{AbstractAlgebra.Generic.MPoly{RingElemType}}
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                  <:MPolyQuoLocalizedRingElem
                                 }
          }
  parent(x) === base_ring(I) || return base_ring(I)(x) in I

  S = parent(x)
  S_flat = _flatten(S)
  I_flat = _flat_counterpart(I)
  return _orig_to_flat(S)(x) in I_flat
end
