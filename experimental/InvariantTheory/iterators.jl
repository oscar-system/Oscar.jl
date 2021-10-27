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

################################################################################
#
#  "Iterate" a vector space
#
################################################################################

function vector_space_iterator(K::FieldT, basis_iterator::IteratorT) where {FieldT <: Union{Nemo.GaloisField, Nemo.GaloisFmpzField, FqNmodFiniteField, FqFiniteField}, IteratorT}
  return VectorSpaceIteratorFiniteField(K, basis_iterator)
end

vector_space_iterator(K::FieldT, basis_iterator::IteratorT, bound::Int = 10^5) where {FieldT, IteratorT} = VectorSpaceIteratorRand(K, basis_iterator, bound)

Base.eltype(VSI::VectorSpaceIterator{FieldT, IteratorT, ElemT}) where {FieldT, IteratorT, ElemT} = ElemT

Base.length(VSI::VectorSpaceIteratorFiniteField) = BigInt(order(VSI.field))^length(VSI.basis_iterator) - 1

# The "generic" iterate for all subtypes of VectorSpaceIterator
function _iterate(VSI::VectorSpaceIterator)
  if isempty(VSI.basis_iterator)
    return nothing
  end

  if length(VSI.basis_collected) > 0
    b = VSI.basis_collected[1]
  else
    b, s = iterate(VSI.basis_iterator)
    VSI.basis_collected = [ b ]
    VSI.basis_iterator_state = s
  end
  phase = length(VSI.basis_iterator) != 1 ? 1 : 3
  return b, (2, Int[ 1, 2 ], phase)
end

Base.iterate(VSI::VectorSpaceIteratorRand) = _iterate(VSI)

function Base.iterate(VSI::VectorSpaceIteratorFiniteField)
  if isempty(VSI.basis_iterator)
    return nothing
  end

  b, state = _iterate(VSI)
  e, s = iterate(VSI.field)
  elts = fill(e, length(VSI.basis_iterator))
  states = [ deepcopy(s) for i = 1:length(VSI.basis_iterator) ]

  return b, (state..., elts, states)
end

# state is supposed to be a tuple of length 3 containing:
# - the next basis element to be returned (relevant for phase 1)
# - a Vector{Int} being the state for phase 2
# - an Int: the number of the phase
function _iterate(VSI::VectorSpaceIterator, state)
  phase = state[3]
  @assert phase == 1 || phase == 2

  if phase == 1
    # Check whether we already computed this basis element
    if length(VSI.basis_collected) >= state[1]
      b = VSI.basis_collected[state[1]]
    else
      b, s = iterate(VSI.basis_iterator, VSI.basis_iterator_state)
      push!(VSI.basis_collected, b)
      VSI.basis_iterator_state = s
    end
    new_phase = state[1] != length(VSI.basis_iterator) ? 1 : 2
    return b, (state[1] + 1, state[2], new_phase)
  end

  # Iterate all possible sums of basis elements
  @assert length(VSI.basis_iterator) > 1
  s = state[2]
  b = sum([ VSI.basis_collected[i] for i in s ])

  expand = true
  if s[end] < length(VSI.basis_collected)
    s[end] += 1
    expand = false
  else
    for i = length(s) - 1:-1:1
      if s[i] + 1 < s[i + 1]
        s[i] += 1
        for j = i + 1:length(s)
          s[j] = s[j - 1] + 1
        end
        expand = false
        break
      end
    end
  end

  if expand
    if length(s) < length(VSI.basis_collected)
      s = collect(1:length(s) + 1)
    else
      phase = 3
    end
  end
  return b, (state[1], s, phase)
end

function Base.iterate(VSI::VectorSpaceIteratorRand, state)
  if state[3] != 3
    return _iterate(VSI, state)
  end

  # Phase 3: Random linear combinations
  coeffs = rand(-VSI.rand_bound:VSI.rand_bound, length(VSI.basis_collected))
  return sum([ coeffs[i]*VSI.basis_collected[i] for i = 1:length(VSI.basis_collected) ]), (state[1], state[2], 3)
end

function Base.iterate(VSI::VectorSpaceIteratorFiniteField, state)
  if state[3] != 3
    b, s = _iterate(VSI, state[1:3])
    return b, (s..., state[4], state[5])
  end

  # Phase 3: Iterate over all possible coefficients in the field

  a = state[4]
  b = state[5]

  n = length(VSI.basis_collected)
  j = n
  ab = iterate(VSI.field, b[j])
  while ab == nothing
    a[j], b[j] = iterate(VSI.field)
    j -= 1
    if j == 0
      return nothing
    end
    ab = iterate(VSI.field, b[j])
  end
  a[j], b[j] = ab[1], ab[2]

  if all( x -> iszero(x) || isone(x), a)
    # We already visited this element in phase 2
    return iterate(VSI, (state[1:3]..., a, b))
  end
  return dot(a, VSI.basis_collected), (state[1:3]..., a, b)
end

function collect_basis(VSI::VectorSpaceIterator)
  while length(VSI.basis_collected) != length(VSI.basis_iterator)
    if !isdefined(VSI, :basis_iterator_state)
      @assert length(VSI.basis_collected) == 0
      b, s = iterate(VSI.basis_iterator)
    else
      b, s = iterate(VSI.basis_iterator, VSI.basis_iterator_state)
    end
    push!(VSI.basis_collected, b)
    VSI.basis_iterator_state = s
  end
  return VSI.basis_collected
end
