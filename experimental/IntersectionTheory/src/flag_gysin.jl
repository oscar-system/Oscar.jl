# Direct type-A Gysin pushforward, Darondeau--Pragacz [DP17], Theorem 1.1.
# There is no flag Groebner basis, normal form, or Vandermonde expansion here.
# We use the formula on E or the annihilator flag on E^vee, whichever omits
# the larger endpoint block. Ties use E^vee, hence ordinary quotient roots.

struct _FlagBundleGysin{D,A,S,P,C,R}
  coefficient_ring::A
  relative_ring::R
  root_ring::S
  chern_images::Vector{P}
  variable_weights::Vector{Int}
  exponents::NTuple{D,Int}
  relative_dimension::Int
  orientation::Symbol
  signed_chern::Vector{C}
  segre::Vector{C}
  cache_lock::ReentrantLock
end

function _FlagBundleGysin(F)
  A, ranks, b = F.coefficient_ring, F.ranks, F.chern
  n, l = sum(ranks), length(ranks)
  use_dual = first(ranks) >= last(ranks)
  omitted = use_dual ? 1 : l
  blocks = use_dual ? (2:l) : ((l - 1):-1:1)
  sign = use_dual ? 1 : -1
  d = n - ranks[omitted]
  # The trivial flag needs neither root variables nor a determinant.
  if d == 0
    T, t = A, elem_type(A)[]
  else
    T, t = polynomial_ring(A, ["t$i" for i in 1:d]; cached=false)
  end
  P, C = elem_type(T), elem_type(A)
  images = Vector{P}(undef, n)
  exponents = Int[]
  offset = 0
  for i in blocks
    r = ranks[i]
    # Elementary symmetric polynomials, without subset enumeration.
    elementary = P[one(T); [zero(T) for _ in 1:r]]
    for k in 1:r
      root = sign*t[offset + k]
      for j in k:-1:1
        elementary[j + 1] += root*elementary[j]
      end
    end
    for j in 1:r
      images[F.starts[i] + j] = elementary[j + 1]
      push!(exponents, n - j)
    end
    offset += r
  end
 
  # c(Q_omitted) = c(E) / prod(1 + sign*t_i*z), truncated to its rank.
  # Ascending k multiplies by the inverse of ONE linear factor at a time.
  # No division in the coefficient ring, and no elimination Groebner basis.
  r = ranks[omitted]
  reciprocal = P[one(T); [zero(T) for _ in 1:r]]
  for root in t
    for k in 1:r
      reciprocal[k + 1] -= sign*root*reciprocal[k]
    end
  end
  for k in 1:r
    value = reciprocal[k + 1]
    for j in 1:k
      value += T(b[j])*reciprocal[k - j + 1]
    end
    images[F.starts[omitted] + k] = value
  end
  signed_chern = C[_flag_coefficient((-sign)^j*b[j]) for j in 1:n]
  weights = reduce(vcat, [collect(1:r) for r in ranks])
  return _FlagBundleGysin{d,typeof(A),typeof(T),P,C,typeof(F.ring)}(
    A, F.ring, T, images, weights, Tuple(exponents),
    sum(ranks .* F.starts), use_dual ? :quotient : :subspace,
    signed_chern, C[one(A)], ReentrantLock(),
  )
end

# The polynomial passed here is relative, with actual base-ring coefficients.
# Powers are cached only within a call: no unbounded cache of user inputs.
function _flag_gysin_roots(G::_FlagBundleGysin{D,A,S,P}, f;
                           prune::Bool=true) where {D,A,S,P}
  @req parent(f) === G.relative_ring "Polynomial belongs to a different relative ring"
  result = zero(G.root_ring)
  powers = Dict{Tuple{Int,Int},P}()
  for (a, e) in zip(coefficients(f), AbstractAlgebra.exponent_vectors(f))
    if prune && sum(e .* G.variable_weights) < G.relative_dimension
      continue
    end
    term = G.root_ring(_flag_coefficient(a))
    iszero(term) && continue
    for j in eachindex(e)
      e[j] == 0 && continue
      power = get!(powers, (j, e[j])) do
        G.chern_images[j]^e[j]
      end
      term *= power
      iszero(term) && break
    end
    result += term
  end
  return result
end

# Return a snapshot. In particular, no reader races with resize!/push! on the
# cache vector. The lock is not a general thread-safety claim about Oscar rings.
function _flag_gysin_segre(G::_FlagBundleGysin, upto::Int)
  return lock(G.cache_lock) do
    n = length(G.signed_chern)
    while length(G.segre) <= upto
      k = length(G.segre)
      value = zero(G.coefficient_ring)
      for j in 1:min(k, n)
        value -= G.signed_chern[j]*G.segre[k - j + 1]
      end
      push!(G.segre, _flag_coefficient(value))
    end
    copy(G.segre[1:(upto + 1)])
  end
end

# Sort rows decreasingly, retaining the determinant sign. Repeated shifts
# give identical rows, so these monomials contribute zero over ANY ring.
function _flag_gysin_rows(alpha, exponents::NTuple{D,Int}) where {D}
  rows = [alpha[i] + D - exponents[i] for i in 1:D]
  sign = 1
  for i in 2:D
    value, j = rows[i], i
    while j > 1 && rows[j - 1] < value
      rows[j] = rows[j - 1]
      j -= 1
      sign = -sign
    end
    rows[j] = value
  end
  for i in 2:D
    rows[i - 1] == rows[i] && return 0, Tuple(rows)
  end
  return sign, Tuple(rows)
end

function _flag_gysin_evaluate(G::_FlagBundleGysin{0}, f)
  @req parent(f) === G.coefficient_ring "Unexpected trivial-flag coefficient ring"
  return _flag_coefficient(f)
end

function _flag_gysin_evaluate(G::_FlagBundleGysin{D,A,S,P,C}, f) where {D,A,S,P,C}
  @req parent(f) === G.root_ring "Polynomial belongs to a different root ring"
  weights = Dict{NTuple{D,Int},C}()
  triangle = div(D*(D + 1), 2)
  for (a, alpha) in zip(coefficients(f), AbstractAlgebra.exponent_vectors(f))
    sum(alpha) < G.relative_dimension && continue
    sign, rows = _flag_gysin_rows(alpha, G.exponents)
    (sign == 0 || rows[end] <= 0) && continue
    # Combine equal determinants BEFORE evaluating any of them. The cache is
    # per call and stores coefficients, not all previous input polynomials.
    value = _flag_coefficient(get(weights, rows, zero(G.coefficient_ring)) + sign*a)
    if iszero(value)
      delete!(weights, rows)
    else
      weights[rows] = value
    end
  end
  isempty(weights) && return zero(G.coefficient_ring)
  upto = maximum(sum(rows) - triangle for rows in keys(weights))
  segre = _flag_gysin_segre(G, upto)
  result = zero(G.coefficient_ring)
  for (rows, a) in weights
    degree = sum(rows) - triangle
    M = zero_matrix(G.coefficient_ring, D, D)
    for i in 1:D, j in 1:D
      k = rows[i] - j
      # Every nonzero determinant product has nonnegative indices summing to
      # degree. Entries above degree cannot contribute, so no larger Segre
      # classes need to be computed, even if individual row shifts are large.
      if 0 <= k <= degree
        M[i, j] = segre[k + 1]
      end
    end
    result += a*Oscar.AbstractAlgebra.det_clow(M)
  end
  return _flag_coefficient(result)
end

function _flag_pushforward(G::_FlagBundleGysin, f)
  return _flag_gysin_evaluate(G, _flag_gysin_roots(G, f))
end
