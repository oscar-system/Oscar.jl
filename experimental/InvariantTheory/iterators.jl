################################################################################
#
#  All Monomials
#
################################################################################

# Iterate over all vectors in Z_{\geq 0}^n of weight (sum of the entries d)
WeightedIntegerVectors(n::T, d::T) where T = WeightedIntegerVectors{T}(n, d)

Base.eltype(WIV::WeightedIntegerVectors{T}) where T = Vector{T}

# We are basically doing d-multicombinations of n here and those are "the same"
# as d-combinations of n + d - 1 (according to Knuth).
function Base.length(WIV::WeightedIntegerVectors{T}) where T
  if iszero(WIV.d)
    return 0
  end
  return binomial(Int(WIV.n + WIV.d - 1), Int(WIV.d))
end

function Base.iterate(WIV::WeightedIntegerVectors{T}, state::Nothing = nothing) where T
  if WIV.n == 0 || WIV.d == 0
    return nothing
  end

  c = zeros(T, WIV.n)
  c[1] = WIV.d
  s = zeros(T, WIV.n)
  if isone(WIV.n)
    s[1] = WIV.d + 1
    return (c, s)
  end

  s[1] = WIV.d - 1
  s[2] = T(1)
  return (c, s)
end

function Base.iterate(WIV::WeightedIntegerVectors{T}, s::Vector{T}) where T
  n = WIV.n
  d = WIV.d
  if s[n] == d + 1
    return nothing
  end
  c = copy(s)
  if s[n] == d
    s[n] += 1
    return (c, s)
  end

  for i = n - 1:-1:1
    if !iszero(s[i])
      s[i] -= 1
      if i + 1 == n
        s[n] += 1
      else
        s[i + 1] = 1
        if !iszero(s[n])
          s[i + 1] += s[n]
          s[n] = 0
        end
      end
      return (c, s)
    end
  end
end

# Returns an iterator over all monomials of R of degree d
function all_monomials(R::MPolyRing, d::Int)
  @assert d >= 0

  if d == 0
    return ( one(R) for i = 1:1)
  end
  return ( set_exponent_vector!(one(R), 1, w) for w in WeightedIntegerVectors(nvars(R), d) )
end

################################################################################
#
#  Bases of Invariant Rings
#
################################################################################

# Returns the dimension of the graded component of degree d.
# If we cannot compute the Molien series (so far in the modular case), we return
# -1.
function dimension_via_molien_series(R::InvRing, d::Int)
  if ismodular(R)
    return -1
  end

  Qt, t = PowerSeriesRing(QQ, d + 1, "t")
  F = molien_series(R)
  k = coeff(numerator(F)(t)*inv(denominator(F)(t)), d)
  @assert isintegral(k)
  k = Int(numerator(k))
  return k
end

# Iterate over the basis of the degree d component of an invariant ring.
# algo can be either :default, :reynolds or :linear_algebra.
function iterate_basis(R::InvRing, d::Int, algo::Symbol = :default)
  @assert d >= 0 "Degree must be non-negativ"

  reynolds = false
  if algo == :reynolds
    reynolds = true
  elseif algo == :linear_algebra
    reynolds = false
  elseif algo == :default
    # TODO: Fine tune this: Depending on d and the group order it is better
    # to use "linear_algebra" also in the non-modular case.
    if ismodular(R)
      reynolds = false
    else
      reynolds = true
    end
  else
    error("Unsupported argument :$(algo) for algo.")
  end

  if reynolds
    return iterate_basis_reynolds(R, d)
  end
  return iterate_basis_linear_algebra(R, d)
end

function iterate_basis_reynolds(R::InvRing, d::Int)
  @assert d >= 0 "Degree must be non-negativ"

  monomials = all_monomials(polynomial_ring(R), d)

  k = dimension_via_molien_series(R, d)
  @assert k != -1

  N = zero_matrix(base_ring(polynomial_ring(R)), 0, 0)
  return InvRingBasisIterator{typeof(R), typeof(monomials), typeof(N)}(R, d, k, true, monomials, N)
end

# Sadly, we can't really do much iteratively here.
function iterate_basis_linear_algebra(IR::InvRing, d::Int)
  @assert d >= 0 "Degree must be non-negativ"

  R = polynomial_ring(IR)

  k = dimension_via_molien_series(IR, d)
  if k == 0
    N = zero_matrix(base_ring(R), 0, 0)
    mons = elem_type(R)[]
    return InvRingBasisIterator{typeof(IR), typeof(mons), typeof(N)}(IR, d, k, false, mons, N)
  end

  mons = collect(all_monomials(R, d))
  if d == 0
    N = identity_matrix(base_ring(R), 1, 1)
    return InvRingBasisIterator{typeof(IR), typeof(mons), typeof(N)}(IR, d, k, false, mons, N)
  end

  mons_to_rows = Dict{elem_type(R), Int}(mons .=> 1:length(mons))

  K = base_ring(R)

  group_gens = action(IR)
  M = zero_matrix(K, length(group_gens)*length(mons), length(mons))
  for i = 1:length(group_gens)
    offset = (i - 1)*length(mons)
    phi = right_action(R, group_gens[i])
    for j = 1:length(mons)
      f = mons[j]
      g = phi(f) - f
      for (c, m) in zip(coefficients(g), monomials(g))
        M[offset + mons_to_rows[m], j] = c
      end
    end
  end
  n, N = right_kernel(M)

  return InvRingBasisIterator{typeof(IR), typeof(mons), typeof(N)}(IR, d, n, false, mons, N)
end

Base.eltype(BI::InvRingBasisIterator) = elem_type(polynomial_ring(BI.R))

Base.length(BI::InvRingBasisIterator) = BI.dim

function Base.iterate(BI::InvRingBasisIterator)
  if BI.reynolds
    return iterate_reynolds(BI)
  end
  return iterate_linear_algebra(BI)
end

function Base.iterate(BI::InvRingBasisIterator, state::InvRingBasisIteratorState)
  if BI.reynolds
    return iterate_reynolds(BI, state)
  end
  return iterate_linear_algebra(BI, state)
end

function iterate_reynolds(BI::InvRingBasisIterator)
  @assert BI.reynolds
  if BI.dim == 0
    return nothing
  end

  K = base_ring(polynomial_ring(BI.R))

  monomial_to_basis = Dict{elem_type(polynomial_ring(BI.R)), Int}()

  if BI.degree == 0
    M = zero_matrix(K, 1, 0)
    state = 1 # Just a dummy
    return one(polynomial_ring(BI.R)), InvRingBasisIteratorState{typeof(M), elem_type(polynomial_ring(BI.R)), Int}(M, monomial_to_basis, state)
  end

  state = nothing
  while true
    fstate = iterate(BI.monomials, state)
    if fstate === nothing
      error("No monomials left")
    end
    f = fstate[1]
    state = fstate[2]

    g = reynolds_operator(BI.R, f)
    if iszero(g)
      continue
    end
    M = zero_matrix(K, 1, length(g))
    i = 1
    for (c, m) in zip(coefficients(g), monomials(g))
      monomial_to_basis[m] = i
      M[1, i] = c
      i += 1
    end
    return inv(leading_coefficient(g))*g, InvRingBasisIteratorState{typeof(M), typeof(g), typeof(state)}(M, monomial_to_basis, state)
  end
end

function iterate_reynolds(BI::InvRingBasisIterator, state::InvRingBasisIteratorState)
  @assert BI.reynolds

  if nrows(state.M) == BI.dim
    return nothing
  end

  K = base_ring(polynomial_ring(BI.R))

  M = state.M
  monomial_state = state.monomial_iter_state
  while true
    fmonomial_state = iterate(BI.monomials, monomial_state)
    if fmonomial_state === nothing
      error("No monomials left")
    end
    f = fmonomial_state[1]
    monomial_state = fmonomial_state[2]

    g = reynolds_operator(BI.R, f)
    if iszero(g)
      continue
    end
    N = vcat(M, zero_matrix(K, 1, ncols(M)))
    new_basis_mon = false
    for (c, m) in zip(coefficients(g), monomials(g))
      if !haskey(state.monomial_to_basis, m)
        new_basis_mon = true
        state.monomial_to_basis[m] = ncols(N) + 1
        N = hcat(N, zero_matrix(K, nrows(N), 1))
      end
      col = state.monomial_to_basis[m]
      N[nrows(N), col] = c
    end
    if new_basis_mon || rref!(N) == nrows(M) + 1
      return inv(leading_coefficient(g))*g, InvRingBasisIteratorState{typeof(N), typeof(g), typeof(monomial_state)}(N, state.monomial_to_basis, monomial_state)
    end
  end
end

function iterate_linear_algebra(BI::InvRingBasisIterator)
  @assert !BI.reynolds
  if BI.dim == 0
    return nothing
  end

  f = polynomial_ring(BI.R)()
  N = BI.kernel
  for i = 1:nrows(N)
    if iszero(N[i, 1])
      continue
    end
    f += N[i, 1]*BI.monomials[i]
  end
  return f, InvRingBasisIteratorState{typeof(N), typeof(f), Int}(2)
end

function iterate_linear_algebra(BI::InvRingBasisIterator, state::InvRingBasisIteratorState)
  @assert !BI.reynolds
  if state.i > BI.dim
    return nothing
  end

  f = polynomial_ring(BI.R)()
  N = BI.kernel
  for i = 1:nrows(N)
    if iszero(N[i, state.i])
      continue
    end
    f += N[i, state.i]*BI.monomials[i]
  end
  return f, InvRingBasisIteratorState{typeof(N), typeof(f), Int}(state.i + 1)
end
