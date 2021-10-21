################################################################################
#
#  All Monomials
#
################################################################################

# Returns an iterator over all monomials of R of degree d
all_monomials(R::MPolyRing, d::Int) = AllMonomials(R, d)

AllMonomials(R::MPolyRing, d::Int) = AllMonomials{typeof(R)}(R, d)

Base.eltype(AM::AllMonomials) = elem_type(AM.R)

# We are basically doing d-multicombinations of n here and those are "the same"
# as d-combinations of n + d - 1 (according to Knuth).
Base.length(AM::AllMonomials) = binomial(nvars(AM.R) + AM.d - 1, AM.d)

function Base.iterate(AM::AllMonomials, state::Nothing = nothing)
  n = nvars(AM.R)
  if AM.d == 0
    s = zeros(Int, n)
    s[n] = AM.d + 1
    return one(AM.R), s
  end

  c = zeros(Int, n)
  c[1] = AM.d
  s = zeros(Int, n)
  if isone(n)
    s[1] = AM.d + 1
    return (set_exponent_vector!(one(AM.R), 1, c), s)
  end

  s[1] = AM.d - 1
  s[2] = 1
  return (set_exponent_vector!(one(AM.R), 1, c), s)
end

function Base.iterate(AM::AllMonomials, s::Vector{Int})
  n = nvars(AM.R)
  d = AM.d
  if s[n] == d + 1
    return nothing
  end
  c = copy(s)
  if s[n] == d
    s[n] += 1
    return (set_exponent_vector!(one(AM.R), 1, c), s)
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
      return (set_exponent_vector!(one(AM.R), 1, c), s)
    end
  end
end

function Base.show(io::IO, AM::AllMonomials)
  println(io, "Iterator over the monomials of degree $(AM.d) of")
  print(io, AM.R)
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
  return InvRingBasisIterator{typeof(R), typeof(monomials), eltype(monomials), typeof(N)}(R, d, k, true, monomials, Vector{eltype(monomials)}(), N)
end

# Sadly, we can't really do much iteratively here.
function iterate_basis_linear_algebra(IR::InvRing, d::Int)
  @assert d >= 0 "Degree must be non-negativ"

  R = polynomial_ring(IR)

  k = dimension_via_molien_series(IR, d)
  if k == 0
    N = zero_matrix(base_ring(R), 0, 0)
    mons = elem_type(R)[]
    dummy_mons = all_monomials(R, 0)
    return InvRingBasisIterator{typeof(IR), typeof(dummy_mons), eltype(mons), typeof(N)}(IR, d, k, false, dummy_mons, mons, N)
  end

  mons_iterator = all_monomials(R, d)
  mons = collect(mons_iterator)
  if d == 0
    N = identity_matrix(base_ring(R), 1)
    dummy_mons = all_monomials(R, 0)
    return InvRingBasisIterator{typeof(IR), typeof(mons_iterator), eltype(mons), typeof(N)}(IR, d, k, false, mons_iterator, mons, N)
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

  return InvRingBasisIterator{typeof(IR), typeof(mons_iterator), eltype(mons), typeof(N)}(IR, d, n, false, mons_iterator, mons, N)
end

Base.eltype(BI::InvRingBasisIterator) = elem_type(polynomial_ring(BI.R))

Base.length(BI::InvRingBasisIterator) = BI.dim

function Base.show(io::IO, BI::InvRingBasisIterator)
  println(io, "Iterator over a basis of the component of degree $(BI.degree) of")
  print(io, BI.R)
end

function Base.iterate(BI::InvRingBasisIterator)
  if BI.reynolds
    return iterate_reynolds(BI)
  end
  return iterate_linear_algebra(BI)
end

function Base.iterate(BI::InvRingBasisIterator, state)
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
    return one(polynomial_ring(BI.R)), (zero_matrix(K, 1, 0), monomial_to_basis, Vector{Int}())
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
    return inv(leading_coefficient(g))*g, (M, monomial_to_basis, state)
  end
end

function iterate_reynolds(BI::InvRingBasisIterator, state)
  @assert BI.reynolds

  M = state[1]
  if nrows(M) == BI.dim
    return nothing
  end

  K = base_ring(polynomial_ring(BI.R))

  monomial_to_basis = state[2]
  monomial_state = state[3]
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
      if !haskey(monomial_to_basis, m)
        new_basis_mon = true
        monomial_to_basis[m] = ncols(N) + 1
        N = hcat(N, zero_matrix(K, nrows(N), 1))
      end
      col = monomial_to_basis[m]
      N[nrows(N), col] = c
    end
    if new_basis_mon || rref!(N) == nrows(M) + 1
      return inv(leading_coefficient(g))*g, (N, monomial_to_basis, monomial_state)
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
    f += N[i, 1]*BI.monomials_collected[i]
  end
  return f, 2
end

function iterate_linear_algebra(BI::InvRingBasisIterator, state::Int)
  @assert !BI.reynolds
  if state > BI.dim
    return nothing
  end

  f = polynomial_ring(BI.R)()
  N = BI.kernel
  for i = 1:nrows(N)
    if iszero(N[i, state])
      continue
    end
    f += N[i, state]*BI.monomials_collected[i]
  end
  return f, state + 1
end
