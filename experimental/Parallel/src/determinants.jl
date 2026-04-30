function modular_det(A::MatrixElem{T};
    wp::Union{OscarWorkerPool, Nothing}=nothing,
    primes_bound=1000000000
  ) where {T<:ZZPolyRingElem}
  m, n = size(A)
  m == n || error("matrix must be square")
  is_zero(m) && return one(ZZ)
  is_one(m) && return A[1, 1]
  bound = factorial(ZZ(m))*ZZ(maximum([maximum(abs.(coefficients(a)); init=0) for a in A]; init=0))^m
  row_degs = [maximum([degree(a) for a in A[i, :]]; init=0) for i in 1:m]
  max_deg = ZZ(sum(row_degs; init=0))
  bound *= binomial(m-1+max_deg, max_deg)

  primes = ZZRingElem[]
  primes_prod = one(ZZ)
  while 2*primes_prod < bound
    p = next_prime(rand(ZZ, 0:primes_bound))
    p in primes && continue
    push!(primes, p)
    primes_prod *= p
  end
  
  Bs = [map_entries(a->map_coefficients(GF(p), a), A) for p in primes]
  dets = isnothing(wp) ? [det(B) for B in Bs] : pmap(det, wp, Tuple(Bs))
  # reconstruct the coefficients
  max_deg = maximum(degree(b) for b in dets; init=0)
  q = prod(primes; init=one(ZZ))
  new_coeffs = sizehint!(ZZRingElem[], max_deg)
  for i in 0:max_deg
    c = crt([lift(ZZ, coeff(b, i)) for b in dets], primes; check=false)
    if 2*c > q
      c -= q
    end
    push!(new_coeffs, c)
  end
  P = base_ring(A)
  return P(new_coeffs)
end


