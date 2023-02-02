########################################################################
# RingFlattening
#
# A struct to accomodate all the information necessary for using 
# flattenings of polynomial rings efficiently to render them effectively
# computable. 
#
# To the outside, this is primarily the map identifying S with its 
# flattened version.
########################################################################
mutable struct RingFlattening{TowerRingType<:MPolyRing, FlatRingType<:Ring, 
                              CoeffRingType<:Ring
                             } <: Hecke.Map{TowerRingType, FlatRingType, 
                                            SetMap, RingFlattening
                                           }
  orig::TowerRingType
  flat::FlatRingType
  R::CoeffRingType
  orig_to_flat::MPolyAnyMap{TowerRingType, FlatRingType}
  flat_to_orig::Hecke.Map{FlatRingType, TowerRingType}
  base_ring_to_flat::Hecke.Map{CoeffRingType, FlatRingType}

  function RingFlattening(
      S::AbstractAlgebra.Generic.MPolyRing{RingElemType}
    ) where {RingElemType <: MPolyElem}
    R = base_ring(S)
    kk = coefficient_ring(R)
    S_flat, _ = PolynomialRing(kk, vcat(symbols(S), symbols(R)))
    R_to_S_flat = hom(R, S_flat, gens(S_flat)[ngens(S)+1:end])
    S_to_S_flat = hom(S, S_flat, R_to_S_flat, gens(S_flat)[1:ngens(S)])
    S_flat_to_S = hom(S_flat, S, vcat(gens(S), S.(gens(R))))
    return new{typeof(S), typeof(S_flat), typeof(R)}(S, S_flat, R, 
                                                     S_to_S_flat, S_flat_to_S,
                                                     R_to_S_flat
                                                    )
  end

  function RingFlattening(
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

    return new{typeof(S), typeof(S_flat), typeof(A)}(S, S_flat, A, 
                                                     S_to_S_flat, S_flat_to_S,
                                                     A_to_S_flat
                                                    )
  end

  function RingFlattening(
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

    return new{typeof(S), typeof(S_flat), typeof(Q)}(S, S_flat, Q, 
                                                     S_to_S_flat, S_flat_to_S,
                                                     Q_to_S_flat
                                                    )
  end

  function RingFlattening(
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

    return new{typeof(S), typeof(S_flat), typeof(L)}(S, S_flat, L, 
                                                     S_to_S_flat, S_flat_to_S,
                                                     L_to_S_flat
                                                    )
  end
end

### Getters 
function (phi::RingFlattening)(x::RingElem)
  return phi.orig_to_flat(x)
end

function inverse(phi::RingFlattening)
  return phi.flat_to_orig
end

function codomain(phi::RingFlattening)
  return phi.flat
end

function domain(phi::RingFlattening)
  return phi.orig
end

function map_from_coefficient_ring_to_flattening(phi::RingFlattening)
  return phi.base_ring_to_flat
end

### Some basic functionality
@attr RingFlattening function flatten(R::MPolyRing)
  return RingFlattening(R)
end

function (phi::RingFlattening)(I::MPolyIdeal)
  if !has_attribute(I, :flat_counterpart)
    I_flat = ideal(codomain(phi), phi.(gens(I)))
    set_attribute!(I, :flat_counterpart, I_flat)
    return I_flat
  end
  return get_attribute(I, :flat_counterpart)
end

function preimage(phi::RingFlattening, x::RingElem)
  parent(x) === codomain(phi) || error("element does not belong to the correct ring")
  return inverse(phi)(x)
end

### Computation of induced morphisms on flattened towers of polynomial rings
@attr function flatten(
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
  flat_S = flatten(S)
  flat_T = flatten(T)
  imgs = flat_T.(f.(gens(S))) # The first half of the images
  imgs = vcat(imgs, flat_T.(gens(coefficient_ring(S)))) # the ones from the base ring
  return hom(codomain(flat_S), codomain(flat_T), imgs)
end

@attr function flatten(
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
  flat_S = flatten(S)
  flat_T = flatten(T)
  RS = coefficient_ring(domain(f))
  
  imgs = flat_T.(f.(gens(S))) # The first half of the images
  imgs = vcat(imgs, flat_T.(coefficient_map(f).(gens(RS)))) # the ones from the base ring
  return hom(codomain(flat_S), codomain(flat_T), imgs)
end

### Deflecting some of the basic methods for ideals in towers of polynomial rings
function ideal_membership(
    x::AbstractAlgebra.Generic.MPoly{RingElemType},
    I::MPolyIdeal{AbstractAlgebra.Generic.MPoly{RingElemType}}
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                  <:MPolyQuoLocalizedRingElem
                                 }
          }
  parent(x) === base_ring(I) || return base_ring(I)(x) in I
  phi = flatten(parent(x))
  return ideal_membership(phi(x), phi(I))
end

function coordinates(
    x::AbstractAlgebra.Generic.MPoly{RingElemType},
    I::MPolyIdeal{AbstractAlgebra.Generic.MPoly{RingElemType}}
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                  <:MPolyQuoLocalizedRingElem
                                 }
          }
  phi = flatten(parent(x))
  return change_base_ring(inverse(phi), (coordinates(phi(x), phi(I))))
end

########################################################################
# Computation of kernels of homomorphisms of polynomial rings 
# over coefficient rings which are themselves polynomial rings and 
# variants thereof.
#
# In the admissible cases, this flattens the structure of the rings 
# and computes the kernel on the flat counterpart. 
########################################################################

# This first signature has to stay to avoid ambiguities.
function kernel(
    f::MPolyAnyMap{<:AbstractAlgebra.Generic.MPolyRing{<:RingElemType}, 
                   <:AbstractAlgebra.Generic.MPolyRing{<:RingElemType}
                  }
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                 <:MPolyQuoLocalizedRingElem
                                }
          }
  f_flat = flatten(f)
  K_flat = kernel(f_flat)
  phi = flatten(domain(f))
  K = ideal(domain(f), inverse(phi).(gens(K_flat)))
  set_attribute!(K, :flat_counterpart, K_flat) # TODO: Support attribute storage?
  return K
end

function kernel(
    f::MPolyAnyMap{<:AbstractAlgebra.Generic.MPolyRing{<:RingElemType}, 
                   <:Ring
                  }
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                 <:MPolyQuoLocalizedRingElem
                                }
          }
  phi = flatten(domain(f))
  f_flat = hom(codomain(phi), codomain(f), f.(inverse(phi).(gens(codomain(phi)))))
  K_flat = kernel(f_flat)
  K = ideal(domain(f), inverse(phi).(gens(K_flat)))
  set_attribute!(K, :flat_counterpart, K_flat) # TODO: Support attribute storage?
  return K
end

function kernel(
    f::MPolyAnyMap{<:Union{<:MPolyRing, <:MPolyQuo, 
                           <:MPolyLocalizedRing, <:MPolyQuoLocalizedRing
                          }, # Restriction necessary to avoid ambiguities
                   <:AbstractAlgebra.Generic.MPolyRing{<:RingElemType}
                  }
  ) where {RingElemType <: Union{<:MPolyElem, <:MPolyQuoElem, <:MPolyLocalizedRingElem, 
                                 <:MPolyQuoLocalizedRingElem
                                }
          }
  phi = flatten(codomain(f))
  f_flat = hom(domain(f), codomain(phi), phi.(f.(gens(domain(f)))))
  return kernel(f_flat)
end
