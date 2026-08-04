########################################################################
# Code for computing determinants of matrices over ZZ[x].
#
# At the moment this is in experimental because of the option to use 
# parallel methods. For the future we might want to port an 
# implementation of this with the parallel methods stripped off to 
# AbstractAlgebra. Comparison of timings and discussions where this 
# can favourably be applied can be found on PR #5971. 
########################################################################
@doc raw"""
    det_modular_coeff_rec(A::MatElem{T};
        wp::Union{OscarWorkerPool, Nothing}=nothing,
        bound::Symbol=:GoldsteinGraham
      ) where {T<:ZZPolyRingElem}

Compute the determinant of `A` via modular coefficient reconstruction, i.e. 
by reducing `A` modulo primes and using the Chinese remainder theorem on the coefficients. 

When `wp` is an `OscarWorkerPool`, the computation of determinants modulo primes is done in parallel on the workers
The symbol `bound` can be set either to `:GoldsteinGraham`, or `:Naive` to select a way to compute a certified bound 
for the selection of primes for reduction.
"""
function det_modular_coeff_rec(A::MatElem{T};
    wp::Union{OscarWorkerPool, Nothing}=nothing,
    bound::Symbol=:GoldsteinGraham
  ) where {T<:ZZPolyRingElem}
  m, n = size(A)
  @req is_square(A) "matrix must be square"
  is_zero(m) && return one(base_ring(A))
  is_one(m) && return A[1, 1]

  bound = _det_height_bound(Val{bound}, A)
  # Select primes until the above bound is reached. 
  primes = Int[next_prime(2^60)]
  primes_prod = ZZ(only(primes))
  bound *= 2 # makes enough space on positive and negative axis
  while primes_prod < bound
    p = next_prime(last(primes)+1)
    push!(primes, p)
    primes_prod *= p
  end
  
  # Compute determinants of the reduced matrices.
  GFs = [Native.GF(p; cached=false) for p in primes]
  Bs = [map_entries(a->map_coefficients(kk, a)::fpPolyRingElem, A) for kk in GFs]
  dets = isnothing(wp) ? [det(B) for B in Bs] : pmap(det, wp, Tuple(Bs))

  # Reconstruct the coefficients from chinese remaindering
  max_deg = maximum(degree(b) for b in dets; init=0)
  max_deg < 0 && return zero(base_ring(A))
  new_coeffs = sizehint!(ZZRingElem[], max_deg+1)
  env = crt_env(ZZ.(primes))
  pp2 = div(primes_prod, 2)
  for i in 0:max_deg
    c = crt([lift(ZZ, coeff(b, i)) for b in dets], env)
    # TODO: We might want to use the `signed` option for `crt` here. 
    # But that is not yet available in Hecke for this signature with 
    # an `crt_env`, so we will postpone this until we have a new 
    # version of Hecke to support this. 
    if c > pp2
      c -= primes_prod
    end
    push!(new_coeffs, c)
  end
  P = base_ring(A)
  return P(new_coeffs)
end

# Functionality for implementing the bound given in
# https://doi.org/10.1137/1016065
#
# This could be more general but hadamard_bound2 is restricted to ZZPolyRingElem
# To apply to a matrix containing QQPolyRingElem, one must clear denoms & convert to ZZPolyRingElem.
function _det_height_bound(::Type{Val{:GoldsteinGraham}}, M::MatElem{ZZPolyRingElem})
  is_square(M) || error("Matrix must be square")
  is_empty(M)  &&  return one(base_ring(M))
  # compute the matrix of 1-norms of the entries
  M2 = map(f->sum(abs(c) for c in coefficients(f); init=zero(ZZ))::ZZRingElem, M)
  B2 = min(hadamard_bound2(M2), hadamard_bound2(transpose(M2)))
  return 1+isqrt(B2)
end

# A naive bound for the absolute value of the coefficients of `det(A)` for a matrix
# over the ring of univariate polynomials
function _det_height_bound(::Type{Val{:Naive}}, A::MatElem{ZZPolyRingElem})
  row_degs = [maximum([degree(a) for a in A[i, :]]; init=0) for i in 1:m]
  max_deg = ZZ(sum(row_degs; init=0))
  # This first bound is for the size of the determinant over ZZ
  # which arises for a single matrix `B` whose rows are the 
  # homogeneous parts of the rows of `A` for fixed degrees `dᵢ` 
  # per row. 
  bound = factorial(ZZ(m))*ZZ(maximum([maximum(abs.(coefficients(a)); init=0) for a in A]; init=0))^m
  # Such determinants add up to the coefficient of tᵈ with as many 
  # summands as there are possibilities for `Σᵢ dᵢ = d`. The following 
  # binomial coefficient is a significant overcount, but reliable 
  # bound of this.
  bound *= binomial(m-1+max_deg, max_deg)
  return bound
end

