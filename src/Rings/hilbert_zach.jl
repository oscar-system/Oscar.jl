########################################################################
# Backend functionality for computation of multivariate Hilbert series 
# due to Matthias Zach
########################################################################


#######################################################################
# 06.10.2022
# The following internal routine can probably still be tuned.
#
# Compared to the implementation in Singular, it is still about 10 
# times slower. We spotted the following possible deficits:
#
#   1) The returned polynomials are not yet built with MPolyBuildCtx.
#      We tried that, but as for now, the code for the build context 
#      is even slower than the direct implementation that is in place 
#      now. Once MPolyBuildCtx is tuned, we should try that again; 
#      the code snippets are still there. In total, building the return 
#      values takes up more than one third of the computation time 
#      at this moment.
#
#   2) Numerous allocations for integer vectors. The code below 
#      performs lots of iterative allocations for lists of integer 
#      vectors. If we constructed a container data structure to 
#      maintain this list internally and do the allocations at once  
#      (for instance in a big matrix), this could significantly 
#      speed up the code. However, that is too much work for 
#      the time being, as long as a high-performing Hilbert series 
#      computation is not of technical importance.
#
#   3) Singular uses bitmasking for exponent vectors to decide on 
#      divisibility more quickly. This is particularly important in 
#      the method _divide_by. For singular, this leads to a speedup 
#      of factor 5. Here, only less than one third of the time is 
#      spent in `_divide_by`, so it does not seem overly important 
#      at this point, but might become relevant in the future. 
#      A particular modification of the singular version of bitmasking 
#      is to compute means for the exponents occurring for each variable 
#      and set the bits depending on whether a given exponent is greater 
#      or less than that mean value. 

function _hilbert_numerator_from_leading_exponents(
    a::Vector{Vector{Int}},
    weight_matrix::Matrix{Int},
    return_ring::Ring,
    #algorithm=:generator
    #algorithm=:custom
    #algorithm=:gcd
    #algorithm=:indeterminate
    #algorithm=:cocoa
    algorithm::Symbol # =:BayerStillmanA, # This is by far the fastest strategy. Should be used.
    # short exponent vectors where the k-th bit indicates that the k-th 
    # exponent is non-zero.
  )
  t = gens(return_ring)
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, return_ring, t)
  ret !== nothing && return ret

  if algorithm == :BayerStillmanA
    return _hilbert_numerator_bayer_stillman(a, weight_matrix, return_ring, t)
  elseif algorithm == :custom
    return _hilbert_numerator_custom(a, weight_matrix, return_ring, t)
  elseif algorithm == :gcd # see Remark 5.3.11
    return _hilbert_numerator_gcd(a, weight_matrix, return_ring, t)
  elseif algorithm == :generator # just choosing on random generator, cf. Remark 5.3.8
    return _hilbert_numerator_generator(a, weight_matrix, return_ring, t)
  elseif algorithm == :indeterminate # see Remark 5.3.8
    return _hilbert_numerator_indeterminate(a, weight_matrix, return_ring, t)
  elseif algorithm == :cocoa # see Remark 5.3.14
    return _hilbert_numerator_cocoa(a, weight_matrix, return_ring, t)
  end
  error("invalid algorithm")
end


########################################################################
# Implementations of the different algorithms below
#
# Compute the Hilbert series from the monomial leading ideal using 
# different bisection strategies along the lines of 
# [Kreuzer, Robbiano: Computational Commutative Algebra 2, Springer]
# Section 5.3.
########################################################################

function _hilbert_numerator_trivial_cases(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector, oneS = one(S)
  )
  length(a) == 0 && return oneS

  # See Proposition 5.3.6
  if _are_pairwise_coprime(a)
    return prod(oneS - _expvec_to_poly(S, t, weight_matrix, e) for e in a)
  end

  return nothing
end

function _hilbert_numerator_cocoa(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector
  )
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t)
  ret !== nothing && return ret

  n = length(a)
  m = length(a[1])
  counters = [0 for i in 1:m]
  for e in a
    counters += [iszero(k) ? 0 : 1 for k in e]
  end
  j = _find_maximum(counters)

  p = a[rand(1:n)]
  while p[j] == 0
    p = a[rand(1:n)]
  end

  q = a[rand(1:n)]
  while q[j] == 0 || p == q
    q = a[rand(1:n)]
  end

  pivot = [0 for i in 1:m]
  pivot[j] = minimum([p[j], q[j]])


  ### Assembly of the quotient ideal with less generators
  rhs = [e for e in a if !_divides(e, pivot)]
  push!(rhs, pivot)

  ### Assembly of the division ideal with less total degree
  lhs = _divide_by_monomial_power(a, j, pivot[j])

  f = one(S)
  for i in 1:nvars(S)
    z = t[i]
    f *= z^(sum([pivot[j]*weight_matrix[i, j] for j in 1:length(pivot)]))
  end

  return _hilbert_numerator_cocoa(rhs, weight_matrix, S, t) + f*_hilbert_numerator_cocoa(lhs, weight_matrix, S, t)
end

function _hilbert_numerator_indeterminate(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector
  )
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t)
  ret !== nothing && return ret

  e = first(a)
  found_at = findfirst(!iszero, e)::Int64
  pivot = zero(e)
  pivot[found_at] = 1

  ### Assembly of the quotient ideal with less generators
  rhs = [e for e in a if e[found_at] == 0]
  push!(rhs, pivot)

  ### Assembly of the division ideal with less total degree
  lhs = _divide_by(a, pivot)

  f = one(S)
  for i in 1:nvars(S)
    z = t[i]
    f *= z^(sum([pivot[j]*weight_matrix[i, j] for j in 1:length(pivot)]))
  end

  return _hilbert_numerator_indeterminate(rhs, weight_matrix, S, t) + f*_hilbert_numerator_indeterminate(lhs, weight_matrix, S, t)
end

function _hilbert_numerator_generator(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector
  )
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t)
  ret !== nothing && return ret

  b = copy(a)
  pivot = pop!(b)

  f = one(S)
  for i in 1:nvars(S)
    z = t[i]
    f *= z^(sum([pivot[j]*weight_matrix[i, j] for j in 1:length(pivot)]))
  end

  c = _divide_by(b, pivot)
  p1 = _hilbert_numerator_generator(b, weight_matrix, S, t)
  p2 = _hilbert_numerator_generator(c, weight_matrix, S, t)

  return p1 - f * p2
end

function _hilbert_numerator_gcd(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector
  )
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t)
  ret !== nothing && return ret

  n = length(a)
  counters = [0 for i in 1:length(a[1])]
  for e in a
    counters += [iszero(k) ? 0 : 1 for k in e]
  end
  j = _find_maximum(counters)

  p = a[rand(1:n)]
  while p[j] == 0
    p = a[rand(1:n)]
  end

  q = a[rand(1:n)]
  while q[j] == 0 || p == q
    q = a[rand(1:n)]
  end

  pivot = _gcd(p, q)

  ### Assembly of the quotient ideal with less generators
  rhs = [e for e in a if !_divides(e, pivot)]
  push!(rhs, pivot)

  ### Assembly of the division ideal with less total degree
  lhs = _divide_by(a, pivot)

  f = one(S)
  for i in 1:nvars(S)
    z = t[i]
    f *= z^(sum([pivot[j]*weight_matrix[i, j] for j in 1:length(pivot)]))
  end

  return _hilbert_numerator_gcd(rhs, weight_matrix, S, t) + f*_hilbert_numerator_gcd(lhs, weight_matrix, S, t)
end

function _hilbert_numerator_custom(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector
  )
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t)
  ret !== nothing && return ret

  p = Vector{Int}()
  q = Vector{Int}()
  max_deg = 0
  for i in 1:length(a)
    b = a[i]
    for j in i+1:length(a)
      c = a[j]
      r = _gcd(b, c)
      if sum(r) > max_deg
        max_deg = sum(r)
        p = b
        q = c
      end
    end
  end

  ### Assembly of the quotient ideal with less generators
  pivot = _gcd(p, q)
  rhs = [e for e in a if !_divides(e, pivot)]
  push!(rhs, pivot)

  ### Assembly of the division ideal with less total degree
  lhs = _divide_by(a, pivot)

  f = one(S)
  for i in 1:nvars(S)
    z = t[i]
    f *= z^(sum([pivot[j]*weight_matrix[i, j] for j in 1:length(pivot)]))
  end

  return _hilbert_numerator_custom(rhs, weight_matrix, S, t) + f*_hilbert_numerator_custom(lhs, weight_matrix, S, t)
end

function _hilbert_numerator_bayer_stillman(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector, oneS = one(S)
  )
  ###########################################################################
  # For this strategy see
  #
  # Bayer, Stillman: Computation of Hilber Series
  # J. Symbolic Computation (1992) No. 14, pp. 31--50
  #
  # Algorithm 2.6, page 35
  ###########################################################################
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t, oneS)
  ret !== nothing && return ret

  # make sure we have lexicographically ordered monomials
  sort!(a, alg=QuickSort)

  # initialize the result 
  h = oneS - _expvec_to_poly(S, t, weight_matrix, a[1])
  linear_mons = Vector{Int}()
  for i in 2:length(a)
    J = _divide_by(a[1:i-1], a[i])
    empty!(linear_mons)
    J1 = filter!(J) do m
      k = findfirst(!iszero, m)
      k === nothing && return false # filter out zero vector
      if m[k] == 1 && findnext(!iszero, m, k + 1) === nothing
        push!(linear_mons, k)
        return false
      end
      return true
    end
    q = _hilbert_numerator_bayer_stillman(J1, weight_matrix, S, t, oneS)
    for k in linear_mons
      q *= (oneS - prod(t[i]^weight_matrix[i, k] for i in 1:length(t)))
    end

    h -= q * _expvec_to_poly(S, t, weight_matrix, a[i])
  end
  return h
end


########################################################################
# Auxiliary helper functions below
########################################################################

### compute t ^ (weight_matrix * expvec), where t == gens(S)
function _expvec_to_poly(S::Ring, t::Vector, weight_matrix::Matrix{Int}, expvec::Vector{Int})
  o = one(coefficient_ring(S))  # TODO: cache this ?!
  return S([o], [weight_matrix * expvec])
end

### special case for univariate polynomial ring
function _expvec_to_poly(S::PolyRing, t::Vector, weight_matrix::Matrix{Int}, expvec::Vector{Int})
  @assert length(t) == 1
  @assert size(weight_matrix) == (1, length(expvec))

  # compute the dot-product of weight_matrix[1,:] and expvec, but faster than dot
  # TODO: what about overflows?
  s = weight_matrix[1,1] * expvec[1]
  for i in 2:length(expvec)
    @inbounds s += weight_matrix[1,i] * expvec[i]
  end
  return t[1]^s
end

function _find_maximum(a::Vector{Int})
  m = a[1]
  j = 1
  for i in 1:length(a)
    if a[i] > m 
      m = a[i]
      j = i
    end
  end
  return j
end

### This implements the speedup from the CoCoA strategy, 
# see Exercise 4 in Section 5.3
function _divide_by_monomial_power(a::Vector{Vector{Int}}, j::Int, k::Int)
  divides(a::Vector{Int}, b::Vector{Int}) = all(>=(0), b[1:j-1] - a[1:j-1]) && all(>=(0), b[j+1:end] - a[j+1:end])
  kept = [e for e in a if e[j] == 0]
  for i in 1:k
    next_slice = [e for e in a if e[j] == i]
    still_kept = next_slice
    for e in kept 
      all(v->!divides(v, e), next_slice) && push!(still_kept, e)
    end
    kept = still_kept
  end
  still_kept = vcat(still_kept, [e for e in a if e[j] > k])

  result = Vector{Vector{Int}}()
  for e in still_kept
    v = copy(e)
    v[j] = (e[j] > k ? e[j] - k : 0)
    push!(result, v)
  end
  return result
end

_divides(a::Vector{Int}, b::Vector{Int}) = all(a .>= b)

function _gcd(a::Vector{Int}, b::Vector{Int})
  return [a[k] > b[k] ? b[k] : a[k] for k in 1:length(a)]
end

function _are_pairwise_coprime(a::Vector{Vector{Int}})
  length(a) <= 1 && return true
  n = length(a)
  m = length(a[1])
  for i in 1:n-1
    for j in i+1:n
      all(k -> a[i][k] == 0 || a[j][k] == 0, 1:m) || return false
    end
  end
  return true
end

### Assume that a is minimal. We divide by `pivot` and return 
# a minimal set of generators for the resulting monomial ideal.
function _divide_by(a::Vector{Vector{Int}}, pivot::Vector{Int})
  # The good ones will contribute to a minimal generating 
  # set of the lhs ideal.
  #
  # The bad monomials come from those which hop over the boundaries of 
  # the monomial diagram by the shift. Their span has a new 
  # generator which is collected in `bad` a priori. It is checked 
  # whether they become superfluous and if not, they are added to 
  # the good ones.
  good = sizehint!(Vector{Vector{Int}}(), length(a))
  bad = sizehint!(Vector{Vector{Int}}(), length(a))
  for e in a
    if _divides(e, pivot)
      push!(good, e - pivot)
    else
      push!(bad, e)
    end
  end

  # pre-allocate m so that don't need to allocate it again each loop iteration
  m = similar(pivot)
  for e in bad
    # the next line computers   m = [k < 0 ? 0 : k for k in e]
    # but without allocations
    for i in 1:length(e)
      m[i] = max(e[i] - pivot[i], 0)
    end

    # check whether the new monomial m is already in the span 
    # of the good ones. If yes, discard it. If not, discard those 
    # elements of the good ones that are in the span of m and put 
    # m in the list of good ones instead. 
    if all(x->!_divides(m, x), good)
      # Remove those 'good' elements which are multiples of m
      filter!(x -> !_divides(x, m), good)
      push!(good, copy(m))
    end
  end

  return good
end

