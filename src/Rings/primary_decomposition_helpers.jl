########################################################################
# Expansion of coefficient fields
#
# For primary decomposition and related operations in a polynomial 
# ring L[x₁,…,xₙ] over a number field L = ℚ[α₁,…,αᵣ] it is better to 
# first translate the respective problem to one in a polynomial ring 
# ℚ [θ₁,…,θᵣ,x₁,…,xₙ] over ℚ by adding further variables for the 
# algebraic elemnts αᵢ and dividing by their minimum polynomials.
#
# This might seem counter intuitive, but it has proved to provide 
# significant speedup in most cases. The reason is that Singular does 
# not have a native data type for most algebraic extension fields 
# used is Oscar and reducing a number field L as above via 
# `absolute_simple_field` destroys sparseness. 
########################################################################

# Since these need the declaration of MPolyQuoRing, they have to come later here 
# in an extra file.
function _expand_coefficient_field(R::MPolyRing{T}; rec_depth::Int=0) where {T<:QQFieldElem}
  return R, identity_map(R), identity_map(R)
end

function _expand_coefficient_field(Q::MPolyQuoRing{<:MPolyRingElem{T}}; rec_depth::Int=0) where {T<:QQFieldElem}
  return Q, identity_map(Q), identity_map(Q)
end

function _expand_coefficient_field(R::MPolyRing{T}; rec_depth=0) where {T<:Union{nf_elem, <:Hecke.NfRelElem}}
  K = coefficient_ring(R)
  alpha = first(gens(K))
  kk = base_field(K)
  P, _ = polynomial_ring(kk, vcat([Symbol("θ_$(rec_depth)")], symbols(R)), cached=false)
  theta = first(gens(P))
  f = defining_polynomial(K)
  d = degree(f)
  powers = elem_type(P)[theta^k for k in 0:d-1]
  function help_map(a)
    return sum(c*powers[i] for (i, c) in enumerate(coefficients(a)); init=zero(P))
  end
  R_flat, pr = quo(P, ideal(P, evaluate(f, theta)))
  to_R_flat = hom(R, R_flat, x->pr(help_map(x)), gens(R_flat)[2:end])
  to_R = hom(R_flat, R, vcat([R(alpha)], gens(R)))
  return R_flat, to_R, to_R_flat
end

# Special dispatch for graded rings to preserve gradings
function _expand_coefficient_field(R::MPolyDecRing{T}; rec_depth=0) where {T<:Union{nf_elem, <:Hecke.NfRelElem}}
  RR = forget_grading(R)
  # We have to do the expansion for RR and then rewrap everything as graded rings/algebras 
  # with appropriate weights
  G = grading_group(R)
  RR_exp, iso, iso_inv = _expand_coefficient_field(RR; rec_depth)
  # RR_exp is now an MPolyQuo because of the modulus given by the minimum polynomial
  PP_exp = base_ring(RR_exp)::MPolyRing
  # Grade the new variable with zero since it belongs to the coefficient field really.
  P_exp, _ = grade(PP_exp, vcat([zero(G)], degree.(gens(R))))
  # Add the modulus again
  gr_mod = ideal(P_exp, P_exp.(gens(modulus(RR_exp))))
  R_exp, _ = quo(P_exp, gr_mod)
  iso_gr = hom(R_exp, R, R.(iso.(gens(RR_exp))))
  coeff_map = coefficient_map(iso_inv)
  @assert coeff_map(one(coefficient_ring(R))) == one(RR_exp)
  iso_inv_gr = hom(R, R_exp, x->R_exp(lift(coeff_map(x))), R_exp.(lift.(iso_inv.(gens(RR)))))
  return R_exp, iso_gr, iso_inv_gr
end

function _expand_coefficient_field(Q::MPolyQuoRing{<:MPolyRingElem{T}}; rec_depth::Int=0) where {T<:Union{nf_elem, <:Hecke.NfRelElem}}
  R = base_ring(Q)
  R_flat, iso, iso_inv = _expand_coefficient_field(R; rec_depth)
  I = modulus(Q)
  I_flat = ideal(R_flat, iso_inv.(gens(I)))
  Q_flat, pr = quo(R_flat, I_flat)
  alpha = first(gens(coefficient_ring(Q)))
  to_Q = hom(Q_flat, Q, vcat([Q(alpha)], gens(Q)))
  theta = Q_flat[1]
  f = defining_polynomial(coefficient_ring(Q))
  d = degree(f)
  powers = [theta^k for k in 0:d-1]
  function help_map(a)
    return sum(c*powers[i] for (i, c) in enumerate(coefficients(a)); init=zero(Q_flat))
  end
  to_Q_flat = hom(Q, Q_flat, help_map, gens(Q_flat)[2:end])
  return Q_flat, to_Q, to_Q_flat
end

function _expand_coefficient_field_to_QQ(R::Union{<:MPolyRing, <:MPolyQuoRing}; rec_depth::Int=0)
  R_flat, to_R, to_R_flat = _expand_coefficient_field(R; rec_depth = rec_depth + 1)
  res, a, b = _expand_coefficient_field_to_QQ(R_flat; rec_depth = rec_depth + 1)
  return res, compose(a, to_R), compose(to_R_flat, b)
end

function _expand_coefficient_field_to_QQ(R::MPolyRing{T}; rec_depth=0) where {T<:QQFieldElem}
  return _expand_coefficient_field(R; rec_depth = rec_depth + 1)
end

function _expand_coefficient_field_to_QQ(R::MPolyQuoRing{<:MPolyRingElem{T}}; rec_depth::Int=0) where {T<:QQFieldElem}
  return _expand_coefficient_field(R; rec_depth = rec_depth + 1)
end

@attr function equidimensional_decomposition_weak(I::MPolyQuoIdeal)
  A = base_ring(I)::MPolyQuoRing
  R = base_ring(A)::MPolyRing
  J = saturated_ideal(I)
  res = equidimensional_decomposition_weak(J)
  return typeof(I)[ideal(A, unique!([x for x in A.(gens(K)) if !iszero(x)])) for K in res]
end


@attr function equidimensional_decomposition_radical(I::MPolyQuoIdeal)
  A = base_ring(I)::MPolyQuoRing
  R = base_ring(A)::MPolyRing
  J = saturated_ideal(I)
  res = equidimensional_decomposition_radical(J)
  return typeof(I)[ideal(A, unique!([x for x in A.(gens(K)) if !iszero(x)])) for K in res]
end

@attr function equidimensional_hull(I::MPolyQuoIdeal)
  A = base_ring(I)::MPolyQuoRing
  R = base_ring(A)::MPolyRing
  J = saturated_ideal(I)
  res = equidimensional_hull(J)
  return ideal(A, unique!([x for x in A.(gens(res)) if !iszero(x)]))
end

@attr function equidimensional_hull_radical(I::MPolyQuoIdeal)
  A = base_ring(I)::MPolyQuoRing
  R = base_ring(A)::MPolyRing
  J = saturated_ideal(I)
  res = equidimensional_hull_radical(J)
  return ideal(A, unique!([x for x in A.(gens(res)) if !iszero(x)]))
end

