########################################################################
# RingFlattening
#
# A struct to accommodate all the information necessary for using 
# flattenings of polynomial rings efficiently to render them effectively
# computable. 
#
# To the outside, this is primarily the map identifying S with its 
# flattened version.
########################################################################
@attributes mutable struct RingFlattening{TowerRingType<:Union{MPolyRing, MPolyQuoRing}, 
                              FlatRingType<:Ring, CoeffRingType<:Ring
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
      S::MPolyRing{RingElemType}
    ) where {RingElemType <: MPolyRingElem}
    R = base_ring(S)
    kk = coefficient_ring(R)
    S_flat, _ = polynomial_ring(kk, vcat(symbols(S), symbols(R)))
    R_to_S_flat = hom(R, S_flat, gens(S_flat)[ngens(S)+1:end])
    S_to_S_flat = hom(S, S_flat, R_to_S_flat, gens(S_flat)[1:ngens(S)])
    S_flat_to_S = hom(S_flat, S, vcat(gens(S), S.(gens(R))))
    return new{typeof(S), typeof(S_flat), typeof(R)}(S, S_flat, R, 
                                                     S_to_S_flat, S_flat_to_S,
                                                     R_to_S_flat
                                                    )
  end

  function RingFlattening(
      S::MPolyRing{RingElemType}
    ) where {RingElemType <: MPolyQuoRingElem}
    A = base_ring(S)::MPolyQuoRing
    R = base_ring(A)::MPolyRing
    I = modulus(A)::MPolyIdeal
    kk = coefficient_ring(R)::Field

    # Before building S_flat, we have to create a polynomial 
    # ring T_flat from which we can pass to the quotient.
    T_flat, _ = polynomial_ring(kk, vcat(symbols(S), symbols(R)))
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
      S::MPolyRing{RingElemType}
    ) where {RingElemType <: MPolyQuoLocRingElem}
    Q = base_ring(S)::MPolyQuoLocRing
    R = base_ring(Q)::MPolyRing
    L = localized_ring(Q)::MPolyLocRing
    A = underlying_quotient(Q)::MPolyQuoRing
    U = inverted_set(Q)::AbsMultSet
    I = modulus(Q)::MPolyLocalizedIdeal
    kk = coefficient_ring(R)::Field

    # Before building S_flat, we have to create a polynomial 
    # ring T_flat from which we can pass to the localized quotient.
    T_flat, _ = polynomial_ring(kk, vcat(symbols(S), symbols(R)))
    R_to_T_flat = hom(R, T_flat, gens(T_flat)[ngens(S)+1:end])

    T_flat_loc, _ = localization(T_flat, R_to_T_flat(U)) # Will throw if multiplicative set can not be transferred. 
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
      S::MPolyRing{RingElemType}
    ) where {RingElemType <: MPolyLocRingElem}
    L = base_ring(S)::MPolyLocRing
    R = base_ring(L)::MPolyRing
    U = inverted_set(L)::AbsMultSet
    kk = coefficient_ring(R)::Field

    # Before building S_flat, we have to create a polynomial 
    # ring T_flat from which we can pass to the localized quotient.
    T_flat, _ = polynomial_ring(kk, vcat(symbols(S), symbols(R)))
    R_to_T_flat = hom(R, T_flat, gens(T_flat)[ngens(S)+1:end])

    S_flat, _ = localization(T_flat, R_to_T_flat(U)) # Will throw if multiplicative set can not be transferred. 
    L_to_S_flat = hom(L, S_flat, gens(S_flat)[ngens(S)+1:end])

    R_to_S_flat = hom(R, S_flat, gens(S_flat)[ngens(S)+1:end])

    S_to_S_flat = hom(S, S_flat, L_to_S_flat, gens(S_flat)[1:ngens(S)])
    S_flat_to_S = hom(S_flat, S, vcat(gens(S), S.(gens(R))))

    return new{typeof(S), typeof(S_flat), typeof(L)}(S, S_flat, L, 
                                                     S_to_S_flat, S_flat_to_S,
                                                     L_to_S_flat
                                                    )
  end

  # Flattenings of quotient rings of the form (𝕜[x][u])/J → 𝕜[x, u]/J'
  function RingFlattening(
      S::MPolyQuoRing{RingElemType}
    ) where {RingElemType <: MPolyRingElem{<:MPolyRingElem}}
    P = base_ring(S) # the free polynomial ring
    P_flattening = flatten(P)
    R = base_ring(P) # the coefficient ring of S
    kk = coefficient_ring(R)

    P_flat = codomain(P_flattening)
    S_flat, _ = quo(P_flat, P_flattening(modulus(S)))

    R_to_S_flat = hom(R, S_flat, gens(S_flat)[ngens(S)+1:end])
    S_to_S_flat = hom(S, S_flat, R_to_S_flat, gens(S_flat)[1:ngens(S)])
    S_flat_to_S = hom(S_flat, S, vcat(gens(S), S.(gens(R))))
    return new{typeof(S), typeof(S_flat), typeof(R)}(S, S_flat, R, 
                                                     S_to_S_flat, S_flat_to_S,
                                                     R_to_S_flat
                                                    )
  end

  # Flattenings of quotient rings of the form ((𝕜[x]/I)[u])/J → 𝕜[x, u]/(I' + J')
  function RingFlattening(
      S::MPolyQuoRing{RingElemType}
    ) where {RingElemType <: MPolyRingElem{<:MPolyQuoRingElem}}
    P = base_ring(S)::MPolyRing # the polynomial ring behind S
    A = base_ring(P)::MPolyQuoRing # the coefficient ring of S

    P_flattening = flatten(P)
    P_flat = codomain(P_flattening)
    mod_flat = P_flattening(modulus(S))
    S_flat, _ = quo(P_flat, mod_flat)

    A_to_S_flat = hom(A, S_flat, gens(S_flat)[ngens(S)+1:end])

    S_to_S_flat = hom(S, S_flat, A_to_S_flat, gens(S_flat)[1:ngens(S)])
    S_flat_to_S = hom(S_flat, S, vcat(gens(S), S.(P.(gens(A)))))

    return new{typeof(S), typeof(S_flat), typeof(A)}(S, S_flat, A, 
                                                     S_to_S_flat, S_flat_to_S,
                                                     A_to_S_flat
                                                    )
  end

  # Flattenings of quotient rings of the form (((𝕜[x]/I)[U⁻¹])[u])/J → (𝕜[x, u]/(I' + J'))[U'⁻¹]
  function RingFlattening(
      S::MPolyQuoRing{RingElemType}
    ) where {RingElemType <: MPolyRingElem{<:MPolyQuoLocRingElem}}
    P = base_ring(S)::MPolyRing
    Q = base_ring(P)::MPolyQuoLocRing
    R = base_ring(Q)::MPolyRing

    P_flattening = flatten(P)
    P_flat = codomain(P_flattening)
    mod_flat = P_flattening(modulus(S))
    S_flat, _ = quo(P_flat, mod_flat)

    Q_to_S_flat = hom(Q, S_flat, gens(S_flat)[ngens(S)+1:end])
    S_to_S_flat = hom(S, S_flat, Q_to_S_flat, gens(S_flat)[1:ngens(S)])

    # The common type of homomorphisms from a localization to S
    # involves computation of the inverse of elements in S. 
    # This has no computational backend. Hence, we need to cheat here. 
    B = base_ring(S_flat)
    psi = hom(B, S, vcat(gens(S), S.(gens(Q))))
    v = vcat([zero(Q) for i in 1:ngens(S)], gens(Q))
    function my_map(a::RingElem) 
      p = lifted_numerator(a)
      q = lifted_denominator(a)
      pp = psi(p)
      qq = S(P(inv(evaluate(q, v)))) # We know that all admissible representatives of denominators must come from elements in L.
      return pp*qq
    end

    #S_flat_to_S = hom(S_flat, S, vcat(gens(S), S.(gens(Q))))
    S_flat_to_S = MapFromFunc(my_map, S_flat, S)

    return new{typeof(S), typeof(S_flat), typeof(Q)}(S, S_flat, Q, 
                                                     S_to_S_flat, S_flat_to_S,
                                                     Q_to_S_flat
                                                    )
  end

  # Flattenings of quotient rings of the form ((𝕜[x][U⁻¹])[u])/J → (𝕜[x, u])[U'⁻¹]/J'
  function RingFlattening(
      S::MPolyQuoRing{RingElemType}
    ) where {RingElemType <: MPolyRingElem{<:MPolyLocRingElem}}
    P = base_ring(S)::MPolyRing
    L = base_ring(P)::MPolyLocRing

    P_flattening = flatten(P)
    P_flat = codomain(P_flattening)
    mod_flat = P_flattening(modulus(S))
    S_flat, _ = quo(P_flat, mod_flat)

    L_to_S_flat = hom(L, S_flat, gens(S_flat)[ngens(S)+1:end])
    S_to_S_flat = hom(S, S_flat, L_to_S_flat, gens(S_flat)[1:ngens(S)])

    # The common type of homomorphisms from a localization to S
    # involves computation of the inverse of elements in S. 
    # This has no computational backend. Hence, we need to cheat here. 
    B = base_ring(S_flat)
    psi = hom(B, S, vcat(gens(S), S.(gens(L))))
    v = vcat([zero(L) for i in 1:ngens(S)], gens(L))
    function my_map(a::RingElem) 
      p = lifted_numerator(a)
      q = lifted_denominator(a)
      pp = psi(p)
      qq = S(P(inv(evaluate(q, v)))) # We know that all admissible representatives of denominators must come from elements in L.
      return pp*qq
    end

    S_flat_to_S = MapFromFunc(my_map, S_flat, S)
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

@attr RingFlattening function flatten(R::MPolyQuoRing)
  return RingFlattening(R)
end

function (phi::RingFlattening)(I::MPolyIdeal)
  if !has_attribute(I, :flat_counterpart)
    I_flat = ideal(codomain(phi), elem_type(codomain(phi))[phi(g) for g in gens(I)])
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
    f::MPolyAnyMap{<:MPolyRing{RingElemType}, 
                   <:MPolyRing{RingElemType},
                   Nothing
                  }
  ) where {RingElemType <: Union{<:MPolyRingElem, <:MPolyQuoRingElem, <:MPolyLocRingElem, 
                                 <:MPolyQuoLocRingElem
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
    f::MPolyAnyMap{<:MPolyRing{RingElemType}, 
                   <:MPolyRing{RingElemType}
                   # Note the missing requirement here: It allows for a non-trivial coefficient map
                  }
  ) where {RingElemType <: Union{<:MPolyRingElem, <:MPolyQuoRingElem, <:MPolyLocRingElem, 
                                 <:MPolyQuoLocRingElem
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
    x::MPolyRingElem{RingElemType},
    I::MPolyIdeal{<:MPolyRingElem{RingElemType}}
  ) where {RingElemType <: Union{<:MPolyRingElem, <:MPolyQuoRingElem, <:MPolyLocRingElem, 
                                  <:MPolyQuoLocRingElem
                                 }
          }
  parent(x) === base_ring(I) || return base_ring(I)(x) in I
  phi = flatten(parent(x))
  return ideal_membership(phi(x), phi(I))
end

function coordinates(
    x::MPolyRingElem{RingElemType},
    I::MPolyIdeal{<:MPolyRingElem{RingElemType}}
  ) where {RingElemType <: Union{<:MPolyRingElem, <:MPolyQuoRingElem, <:MPolyLocRingElem, 
                                  <:MPolyQuoLocRingElem
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
    f::MPolyAnyMap{<:MPolyRing{<:RingElemType}, 
                   <:MPolyRing{<:RingElemType}
                  }
  ) where {RingElemType <: Union{<:MPolyRingElem, <:MPolyQuoRingElem, <:MPolyLocRingElem, 
                                 <:MPolyQuoLocRingElem
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
    f::MPolyAnyMap{<:MPolyRing{<:RingElemType}, 
                   <:Ring
                  }
  ) where {RingElemType <: Union{<:MPolyRingElem, <:MPolyQuoRingElem, <:MPolyLocRingElem, 
                                 <:MPolyQuoLocRingElem
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
    f::MPolyAnyMap{<:Union{<:MPolyRing, <:MPolyQuoRing, 
                           <:MPolyLocRing, <:MPolyQuoLocRing
                          }, # Restriction necessary to avoid ambiguities
                   <:MPolyRing{<:RingElemType}
                  }
  ) where {RingElemType <: Union{<:MPolyRingElem, <:MPolyQuoRingElem, <:MPolyLocRingElem, 
                                 <:MPolyQuoLocRingElem
                                }
          }
  phi = flatten(codomain(f))
  f_flat = hom(domain(f), codomain(phi), phi.(f.(gens(domain(f)))))
  return kernel(f_flat)
end

### saturations
function saturation(
    I::MPolyIdeal{T}, J::MPolyIdeal{T}
  ) where {T<:MPolyRingElem{<:Union{<:MPolyRingElem, MPolyQuoRingElem, MPolyLocRingElem, 
                                <:MPolyQuoLocRingElem
                               }}}
  S = base_ring(I)
  phi = flatten(S)
  pre_res = saturation(phi(I), phi(J))
  return ideal(S, inverse(phi).(gens(pre_res)))
end

### transferred functionality for quotient rings
function is_invertible_with_inverse(
    a::MPolyQuoRingElem{<:MPolyRingElem{<:Union{MPolyRingElem, MPolyQuoRingElem, 
                                            MPolyQuoLocRingElem, MPolyLocRingElem}
                                   }})
  phi = flatten(parent(a))
  a_flat = phi(a)
  is_unit(a_flat) || return false, a
  return true, inverse(phi)(inv(a_flat))
end

function is_unit(
    a::MPolyQuoRingElem{<:MPolyRingElem{<:Union{MPolyRingElem, MPolyQuoRingElem, 
                                            MPolyQuoLocRingElem, MPolyLocRingElem}
                                   }})
  return is_unit(flatten(parent(a))(a))
end

function inv(
    a::MPolyQuoRingElem{<:MPolyRingElem{<:Union{MPolyRingElem, MPolyQuoRingElem, 
                                            MPolyQuoLocRingElem, MPolyLocRingElem}
                                   }})
  phi = flatten(parent(a))
  a_flat = phi(a)
  return inverse(phi)(inv(a_flat))
end

