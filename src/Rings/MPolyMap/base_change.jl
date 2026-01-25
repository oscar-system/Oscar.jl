function change_base_ring(A::Ring, P::MPolyRing)
  kk = coefficient_ring(P)
  res, _ = polynomial_ring(A, symbols(P); cached = false)
  return res, hom(P, res, A, gens(res), check=false)
end

function change_base_ring(A::Ring, Q::MPolyQuoRing)
  P = base_ring(Q)
  PA, phi = change_base_ring(A, P)
  I = phi(modulus(Q))
  res, pr = quo(PA, I)
  return res, hom(Q, res, A, gens(res), check=false)
end

# Necessary duplicate because of method ambiguity
function change_base_ring(A::Ring, P::MPolyDecRing)
  kk = coefficient_ring(P)
  res, _ = graded_polynomial_ring(A, symbols(P); weights=generator_degrees(P), cached = false)
  return res, hom(P, res, A, gens(res), check=false)
end

function change_base_ring(phi::Any, Q::MPolyQuoRing)
  P = base_ring(Q)::MPolyRing
  PP, psi = change_base_ring(phi, P)
  I = psi(modulus(Q))
  res, pr = quo(PP, I)
  return res, hom(Q, res, phi, gens(res), check=false)
end

