function modular_det(A::MatrixElem{T};
    wp::Union{OscarWorkerPool, Nothing}=nothing,
    primes_bound=ZZ(1000)^6 # Experiments suggest that beyond this bound
                            # modular arithmetic becomes slow again.
  ) where {T<:ZZPolyRingElem}
  m, n = size(A)
  m == n || error("matrix must be square")
  is_zero(m) && return one(ZZ)
  is_one(m) && return A[1, 1]

  # Determine a bound for the size of coefficients of the result.
  row_degs = [maximum([degree(a) for a in A[i, :]]; init=0) for i in 1:m]
  max_deg = ZZ(sum(row_degs; init=0))

  # This first bound is for the size of the determinant over ZZ
  # which arises for a single matrix `B` whose rows are the 
  # homogeneous parts of the rows of `A` for fixed degrees `dᵢ` 
  # per row. 
  bound = factorial(ZZ(m))*ZZ(maximum([maximum(abs.(coefficients(a)); init=0) for a in A]; init=0))^m
  # Such determinants add up to the coefficient of tᵈwith as many 
  # summands as there are possibilities for `Σᵢ dᵢ = d`. The following 
  # binomial coefficient is a significant overcount, but reliable 
  # bound of this.
  bound *= binomial(m-1+max_deg, max_deg)

  # Select random primes until the above bound is reached. 
  primes = ZZRingElem[]
  primes_prod = one(ZZ)
  while 2*primes_prod < bound
    p = next_prime(rand(ZZ, 0:primes_bound))
    p in primes && continue
    push!(primes, p)
    primes_prod *= p
  end
  
  # Compute determinants of the reduced matrices.
  Bs = [map_entries(a->map_coefficients(GF(p), a)::FqPolyRingElem, A) for p in primes]
  dets = isnothing(wp) ? [det(B) for B in Bs] : pmap(det, wp, Tuple(Bs))

  # Reconstruct the coefficients from chinese remaindering
  max_deg = maximum(degree(b) for b in dets; init=0)
  new_coeffs = sizehint!(ZZRingElem[], max_deg)
  for i in 0:max_deg
    c = crt([lift(ZZ, coeff(b, i)) for b in dets], primes; check=false)
    if 2*c > primes_prod
      c -= primes_prod
    end
    push!(new_coeffs, c)
  end
  P = base_ring(A)
  return P(new_coeffs)
end


