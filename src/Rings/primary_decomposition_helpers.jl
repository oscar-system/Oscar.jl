# Since these need the declaration of MPolyQuoRing, they have to come later here 
# in an extra file.
function _flatten_ring(R::MPolyRing{T}; rec_depth::Int=0) where {T<:QQFieldElem}
  return R, identity_map(R), identity_map(R)
end

function _flatten_ring(Q::MPolyQuoRing{<:MPolyRingElem{T}}; rec_depth::Int=0) where {T<:QQFieldElem}
  return Q, identity_map(Q), identity_map(Q)
end

function _flatten_ring(R::MPolyRing{T}; rec_depth=0) where {T<:Union{nf_elem, <:Hecke.NfRelElem}}
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

function _flatten_ring(Q::MPolyQuoRing{<:MPolyRingElem{T}}; rec_depth::Int=0) where {T<:Union{nf_elem, <:Hecke.NfRelElem}}
  R = base_ring(Q)
  R_flat, iso, iso_inv = _flatten_ring(R; rec_depth)
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

function _flatten_to_QQ(R::Union{<:MPolyRing, <:MPolyQuoRing}; rec_depth::Int=0)
  R_flat, to_R, to_R_flat = _flatten_ring(R; rec_depth = rec_depth + 1)
  res, a, b = _flatten_to_QQ(R_flat; rec_depth = rec_depth + 1)
  return res, compose(a, to_R), compose(to_R_flat, b)
end

function _flatten_to_QQ(R::MPolyRing{T}; rec_depth=0) where {T<:QQFieldElem}
  return _flatten_ring(R; rec_depth = rec_depth + 1)
end

function _flatten_to_QQ(R::MPolyQuoRing{<:MPolyRingElem{T}}; rec_depth::Int=0) where {T<:QQFieldElem}
  return _flatten_ring(R; rec_depth = rec_depth + 1)
end
