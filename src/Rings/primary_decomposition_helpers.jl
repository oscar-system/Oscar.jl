########################################################################
# Expansion of coefficient fields
#
# For primary decomposition and related operations in a polynomial 
# ring L[x₁,…,xₙ] over a number field L = ℚ[α₁,…,αᵣ] it is better to 
# first translate the respective problem to one in a polynomial ring 
# ℚ [θ₁,…,θᵣ,x₁,…,xₙ] over ℚ by adding further variables for the 
# algebraic elements αᵢ and dividing by their minimum polynomials.
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
  return R, id_hom(R), id_hom(R)
end

function _expand_coefficient_field(Q::MPolyQuoRing{<:MPolyRingElem{T}}; rec_depth::Int=0) where {T<:QQFieldElem}
  return Q, id_hom(Q), id_hom(Q)
end

function _expand_coefficient_field(R::MPolyRing{T}; rec_depth=0) where {T<:Union{AbsSimpleNumFieldElem, <:Hecke.RelSimpleNumFieldElem}}
  K = coefficient_ring(R)
  alpha = first(gens(K))
  kk = base_field(K)
  P, _ = polynomial_ring(kk, vcat([Symbol("θ_$(rec_depth)")], symbols(R)); cached = false)
  theta = first(gens(P))
  f = defining_polynomial(K)
  d = degree(f)
  powers = elem_type(P)[theta^k for k in 0:d-1]
  R_flat, pr = quo(P, ideal(P, evaluate(f, theta)))
  to_R_flat = hom(R, R_flat, hom(K, R_flat, pr(theta); check=false), gens(R_flat)[2:end]; check=false)
  to_R = hom(R_flat, R, vcat([R(alpha)], gens(R)); check=false)
  return R_flat, to_R, to_R_flat
end

function _expand_coefficient_field(
    R::MPolyRing{T}; rec_depth=0
  ) where {T<:Union{<:AbsNonSimpleNumFieldElem, <:Hecke.RelNonSimpleNumFieldElem}}
  K = coefficient_ring(R)
  alpha = gens(K)
  r = length(alpha)
  kk = base_field(K)
  P, _ = polynomial_ring(kk, vcat([Symbol("θ_$(rec_depth)_$i") for i in 1:r], symbols(R)); cached = false)
  theta = gens(P)[1:r]
  f = defining_polynomials(K)
  d = degree.(f)
  R_flat, pr = quo(P, ideal(P, [evaluate(a, b) for (a, b) in zip(f, theta)]))
  to_R_flat = hom(R, R_flat, hom(K, R_flat, pr.(theta); check=false), gens(R_flat)[r+1:end]; check=false)
  to_R = hom(R_flat, R, vcat(R.(alpha), gens(R)); check=false)
  return R_flat, to_R, to_R_flat
end

function _expand_coefficient_field(
    A::MPolyQuoRing{S}; rec_depth::Int=0
  ) where {T<:Union{<:AbsNonSimpleNumFieldElem, <:Hecke.RelNonSimpleNumFieldElem}, S<:MPolyRingElem{T}}
  R = base_ring(A)
  R_exp, iso, iso_inv = _expand_coefficient_field(R; rec_depth)
  I = ideal(R_exp, iso_inv.(gens(modulus(A))))
  A_exp, pr = quo(R_exp, I)
  r = ngens(R_exp) - ngens(R) # The first r variables have been added for the field extension
  theta = gens(A_exp)[1:r]
  alpha = gens(coefficient_ring(A))
  to_A = hom(A_exp, A, vcat(A.(alpha), gens(A)); check=false)
  to_A_exp = hom(A, A_exp, hom(coefficient_ring(A), A_exp, theta; check=false), gens(A_exp)[r+1:end]; check=false)
  return A_exp, to_A, to_A_exp
end

# Special dispatch for graded rings to preserve gradings
function _expand_coefficient_field(R::MPolyDecRing{T}; rec_depth::Int=0) where {T<:Union{AbsSimpleNumFieldElem, <:Hecke.RelSimpleNumFieldElem, <:AbsNonSimpleNumFieldElem, <:Hecke.RelNonSimpleNumFieldElem}}
  RR = forget_grading(R)
  # We have to do the expansion for RR and then rewrap everything as graded rings/algebras 
  # with appropriate weights
  G = grading_group(R)
  RR_exp, iso, iso_inv = _expand_coefficient_field(RR; rec_depth)
  # RR_exp is now an MPolyQuo because of the modulus given by the minimum polynomial
  PP_exp = base_ring(RR_exp)::MPolyRing
  r = ngens(PP_exp) - ngens(RR) # The number of variables added for the field extension
  # Grade the new variable with zero since it belongs to the coefficient field really.
  P_exp, _ = grade(PP_exp, vcat([zero(G) for i in 1:r], degree.(gens(R))))
  # Add the modulus again
  gr_mod = ideal(P_exp, P_exp.(gens(modulus(RR_exp))))
  R_exp, _ = quo(P_exp, gr_mod)
  iso_gr = hom(R_exp, R, R.(iso.(gens(RR_exp))); check=false)
  coeff_map = coefficient_map(iso_inv)
  @assert coeff_map(one(coefficient_ring(R))) == one(RR_exp) "coefficient map is incorrect"
  iso_inv_gr = hom(R, R_exp, x->R_exp(lift(coeff_map(x))), R_exp.(lift.(iso_inv.(gens(RR)))); check=false)
  return R_exp, iso_gr, iso_inv_gr
end

function _expand_coefficient_field(Q::MPolyQuoRing{<:MPolyRingElem{T}}; rec_depth::Int=0) where {T<:Union{<:AbsSimpleNumFieldElem, <:Hecke.RelSimpleNumFieldElem}}
  R = base_ring(Q)
  R_flat, iso, iso_inv = _expand_coefficient_field(R; rec_depth)
  I = modulus(Q)
  I_flat = ideal(R_flat, iso_inv.(gens(I)))
  Q_flat, pr = quo(R_flat, I_flat)
  alpha = first(gens(coefficient_ring(Q)))
  to_Q = hom(Q_flat, Q, vcat([Q(alpha)], gens(Q)); check=false)
  theta = Q_flat[1]
  f = defining_polynomial(coefficient_ring(Q))
  d = degree(f)
  powers = [theta^k for k in 0:d-1]
  to_Q_flat = hom(Q, Q_flat, hom(coefficient_ring(Q), Q_flat, theta; check=false), gens(Q_flat)[2:end]; check=false)
  return Q_flat, to_Q, to_Q_flat
end

function _expand_coefficient_field_to_QQ(R::Union{<:MPolyRing, <:MPolyQuoRing}; rec_depth::Int=0)
  get_attribute!(R, :coefficient_field_expansion) do
    R_flat, to_R, to_R_flat = _expand_coefficient_field(R; rec_depth = rec_depth + 1)
    res, a, b = _expand_coefficient_field_to_QQ(R_flat; rec_depth = rec_depth + 1)
    return res, compose(a, to_R), compose(to_R_flat, b)
  end::Tuple{<:Ring, <:Map, <:Map}
end

function _expand_coefficient_field_to_QQ(R::MPolyRing{T}; rec_depth=0) where {T<:QQFieldElem}
  get_attribute!(R, :coefficient_field_expansion) do
    return _expand_coefficient_field(R; rec_depth = rec_depth + 1)
  end::Tuple{<:Ring, <:Map, <:Map}
end

function _expand_coefficient_field_to_QQ(R::MPolyQuoRing{<:MPolyRingElem{T}}; rec_depth::Int=0) where {T<:QQFieldElem}
  get_attribute!(R, :coefficient_field_expansion) do
    return _expand_coefficient_field(R; rec_depth = rec_depth + 1)
  end::Tuple{<:Ring, <:Map, <:Map}
end

@attr Any function equidimensional_decomposition_weak(I::MPolyQuoIdeal)
  A = base_ring(I)::MPolyQuoRing
  R = base_ring(A)::MPolyRing
  J = saturated_ideal(I)
  res = equidimensional_decomposition_weak(J)
  return typeof(I)[ideal(A, unique!([x for x in A.(gens(K)) if !iszero(x)])) for K in res]
end


@attr Any function equidimensional_decomposition_radical(I::MPolyQuoIdeal)
  A = base_ring(I)::MPolyQuoRing
  R = base_ring(A)::MPolyRing
  J = saturated_ideal(I)
  res = equidimensional_decomposition_radical(J)
  return typeof(I)[ideal(A, unique!([x for x in A.(gens(K)) if !iszero(x)])) for K in res]
end

@attr Any function equidimensional_hull(I::MPolyQuoIdeal)
  A = base_ring(I)::MPolyQuoRing
  R = base_ring(A)::MPolyRing
  J = saturated_ideal(I)
  res = equidimensional_hull(J)
  return ideal(A, unique!([x for x in A.(gens(res)) if !iszero(x)]))
end

@attr Any function equidimensional_hull_radical(I::MPolyQuoIdeal)
  A = base_ring(I)::MPolyQuoRing
  R = base_ring(A)::MPolyRing
  J = saturated_ideal(I)
  res = equidimensional_hull_radical(J)
  return ideal(A, unique!([x for x in A.(gens(res)) if !iszero(x)]))
end

@attr Any function absolute_primary_decomposition(I::MPolyQuoIdeal)
  A = base_ring(I)::MPolyQuoRing
  R = base_ring(A)::MPolyRing
  J = saturated_ideal(I)
  res = absolute_primary_decomposition(J)
  result = []
  if is_empty(res)
    U = ideal(A, one(A))
    return [(U, U, U, 0)]
  end
  
  # Create the ring for the return values
  for (P, Q, P_prime, d) in res
    R_prime = base_ring(P_prime) # the new polynomial ring
    L = coefficient_ring(R_prime) # the new field for the result
    A_ext, ext_map = change_base_ring(L, A) # recreate the quo-ring over that field
    @assert coefficient_ring(base_ring(A_ext)) === L
    help_map = hom(R_prime, A_ext, gens(A_ext); check=false)
    PP = ideal(A, A.(gens(P)))
    QQ = ideal(A, A.(gens(Q)))
    trans_gens = help_map.(gens(P_prime))
    PP_prime = ideal(A_ext, trans_gens)
    push!(result, (PP, QQ, PP_prime, d))
  end
  return result
end

#=
function change_base_ring(phi::Any, A::MPolyQuoRing)
  R = base_ring(A)
  RR, map_R = change_base_ring(phi, R)
  I = ideal(RR, map_R.(gens(modulus(A))))
  AA, pr = quo(RR, I)
  psi = hom(A, AA, phi, gens(AA), check=false)
  return AA, psi
end
=#

function change_base_ring(phi::Any, R::MPolyRing)
  kk = coefficient_ring(R)
  L = parent(phi(zero(kk)))
  RR, _ = polynomial_ring(L, symbols(R); cached = false)
  psi = hom(R, RR, phi, gens(RR); check=false)
  return RR, psi
end

function change_base_ring(phi::Any, R::MPolyDecRing)
  kk = coefficient_ring(R)
  L = parent(phi(zero(kk)))
  RR, _ = graded_polynomial_ring(L, symbols(R); cached = false, weights=weights(R))
  psi = hom(R, RR, phi, gens(RR); check=false)
  return RR, psi
end
