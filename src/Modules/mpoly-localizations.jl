########################################################################
#
# Localizations of finitely generated modules over multivariate 
# polynomial rings.
#
# This file implements the required functionality for the 
# Posur-interface; see [1].
#
# [1] Posur: Linear systems over localizations of rings, arXiv:1709.08180v2


function has_nonempty_intersection(U::MPolyPowersOfElement, I::MPolyIdeal; check::Bool=true)
  R = ring(U)
  R == base_ring(I) || error("the multiplicative set and the ideal must be defined over the same ring")

  d = prod(denominators(U); init=one(R))
  if check
    inradical(d, I) || return false, zero(R), zero_matrix(R, 1, ngens(I))
  end
  (k, f) = _minimal_power_such_that(d, in(I))
  return true, f, coordinates(f, I)
end

function has_nonempty_intersection(U::MPolyComplementOfPrimeIdeal, I::MPolyIdeal; check::Bool=true)
  R = ring(U)
  R == base_ring(I) || error("the multiplicative set and the ideal must be defined over the same ring")
  P = prime_ideal(U)
  candidates = [(f, i) for (f, i) in zip(gens(I), 1:ngens(I)) if !(f in P)]
  length(candidates) == 0 && return false, zero(R), zero_matrix(R, 1, ngens(I))
  d = maximum([total_degree(f) for (f, i) in candidates])
  (g, j) = candidates[1]
  for (h, k) in candidates
    if total_degree(h) < total_degree(g) 
      (g, j) = (h, k)
    end
  end
  A = zero_matrix(R, 1, ngens(I))
  A[1, j] = 1
  return true, g, A
end

function has_nonempty_intersection(U::MPolyComplementOfKPointIdeal, I::MPolyIdeal; check::Bool=true)
  R = ring(U)
  R == base_ring(I) || error("the multiplicative set and the ideal must be defined over the same ring")
  a = point_coordinates(U)
  candidates = [(f, i) for (f, i) in zip(gens(I), 1:ngens(I)) if !(iszero(evaluate(f, a)))]
  length(candidates) == 0 && return false, zero(R), zero_matrix(R, 1, ngens(I))
  d = maximum([total_degree(f) for (f, i) in candidates])
  (g, j) = candidates[1]
  for (h, k) in candidates
    if total_degree(h) < total_degree(g) 
      (g, j) = (h, k)
    end
  end
  A = zero_matrix(R, 1, ngens(I))
  A[1, j] = 1
  return true, g, A
end

function has_nonempty_intersection(U::MPolyProductOfMultSets, I::MPolyIdeal; check::Bool=true)
  J = I
  R = ring(U) 
  R == base_ring(I) || error("rings not compatible")
  Usets = sets(U)
  if length(Usets) == 1 
    return Oscar.has_nonempty_intersection(Usets[1], I, check=check)
  end

  V = pop!(Usets)
  Iloc = MPolyLocRing(R, V)(I)
  saturated_ideal(Iloc, with_generator_transition=true)
  J = pre_saturated_ideal(Iloc)
  (success, g, A) = has_nonempty_intersection(MPolyProductOfMultSets(R, Usets), J, check=check)
  if !success 
    return false, zero(R), zero_matrix(R, 1, ngens(I))
  end
  T = pre_saturation_data(Iloc)
  Bext = transpose(T * transpose(A))
  #Bext = A*T
  u = lcm(_vec(denominator.(Bext)))
  B = map_entries(x->preimage(map_from_base_ring(Iloc), x), u*Bext)
  return true, u*g, B
end

# For a `RingElem` f this computes a pair (k, h) where h = f^k and 
# k is the minimal natural number such that the property P(f^k) is 
# satisfied. 
#
# The algorithm is a simple implementation of logarithmic bisection.
function _minimal_power_such_that(f::RingElemType, P::PropertyType) where {RingElemType<:RingElem, PropertyType}
  P(one(parent(f))) && return (0, one(f))
  P(f) && return (1, f)
  f_powers = [(1,f)]

  while !P(last(f_powers)[2])
    push!(f_powers, (last(f_powers)[1]*2, last(f_powers)[2]^2))
  end
  upper = pop!(f_powers)
  lower = pop!(f_powers)
  while upper[1]!=lower[1]+1
    middle = pop!(f_powers)
    middle = (lower[1]+middle[1], lower[2]*middle[2])
    if P(middle[2])
      upper = middle
    else
      lower = middle
    end
  end
  return upper
end

########################################################################
# Background material to make the use of local orderings available     #
########################################################################
@doc raw"""
    shifted_module(
        F::FreeMod{T}
      ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                                   <:MPolyComplementOfKPointIdeal}}

For a free module ``F`` over a localized polynomial ring ``Rₘ`` at a maximal ideal ``𝔪`` 
of a rational point this returns a triple ``(F♭, φ, φ⁻¹)`` where ``F♭`` is the 
corresponding free module over ``R`` and ``φ : F♭ → F♭`` is the isomorphism over 
the shift map ``Φ : R → R`` which is moving the point of ``𝔪`` to the origin.
"""
@attr Any function shifted_module(
    F::FreeMod{T}
  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}

  (a, b) = base_ring_shifts(F)
  return base_ring_module(F), a, b
end

@attr Any function base_ring_shifts(
    F::FreeMod{T}
  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}
  Fb = base_ring_module(F)
  L = base_ring(F)
  R = base_ring(L)
  shift, back_shift = base_ring_shifts(L)
  F_shift = hom(Fb, Fb, gens(Fb), shift)
  F_backshift = hom(Fb, Fb, gens(Fb), back_shift)
  return F_shift, F_backshift
end

@doc raw"""
    shifted_module(
        M::SubquoModule{T}
      ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                                   <:MPolyComplementOfKPointIdeal}}

For a subquotient module ``M`` over a localized polynomial ring ``Rₘ`` at a maximal ideal ``𝔪`` 
of a rational point and a `pre_saturated_module` ``N`` over ``R``, this returns a triple 
``(N', φ, φ⁻¹)`` where ``N'`` is a module over ``R``, and ``φ : N → N'`` is an isomorphism over 
the shift map ``Φ : R → R`` which is moving the point of ``𝔪`` to the origin.
"""
@attr Any function shifted_module(
    M::SubquoModule{T}
  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}

  R = base_ring(M)::MPolyLocRing
  P = base_ring(R)::MPolyRing

  Mp = pre_saturated_module(M)
  F = ambient_free_module(M)
  shift, backshift = base_ring_shifts(R)
  Fb, F_shift, F_backshift = shifted_module(ambient_free_module(M))
  Mp_sub = F_shift.(ambient_representatives_generators(Mp))
  Mp_rel = F_shift.(relations(Mp))
  result = quo_object(sub_object(Fb, Mp_sub), Mp_rel)
  a = hom(Mp, result, gens(result), shift)
  b = hom(result, Mp, gens(Mp), backshift)
  return result, a, b
end

@doc raw"""
    shifted_module(
        M::SubModuleOfFreeModule{T}
                  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem,
                                <:MPolyRing, <:MPolyRingElem,
                                <:MPolyComplementOfKPointIdeal}}
For a submodule ``M`` of a free module over a localized polynomial ring ``Rₘ`` at a maximal
ideal ``𝔪`` of a rational point and a `pre_saturated_module` ``N`` over ``R``, this returns a triple
``(N', φ, φ⁻¹)`` where ``N'`` is a module over ``R``, and ``φ : N → N'`` is an isomorphism over
the shift map ``Φ : R → R`` which is moving the point of ``𝔪`` to the origin.
"""
@attr Any function shifted_module(
    M::SubModuleOfFreeModule{T}
  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem,
                               <:MPolyComplementOfKPointIdeal}}

  R = base_ring(M)::MPolyLocRing
  P = base_ring(R)::MPolyRing
  FL = ambient_free_module(M)
  F = base_ring_module(FL)
  (A,D) = clear_denominators(M.matrix)
  Mp = SubModuleOfFreeModule(F,A)
  F_shifted, shift, shift_back = shifted_module(FL)
  Mp_gens_shift = shift.(gens(Mp))
  result = SubModuleOfFreeModule(F_shifted,Mp_gens_shift)
  return result
end
