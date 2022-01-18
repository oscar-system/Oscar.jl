export iterate_basis

################################################################################
#
#  All Monomials
#
################################################################################

# Return an iterator over all monomials of R of degree d
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

# Return the dimension of the graded component of degree d.
# If we cannot compute the Molien series (so far in the modular case), we return
# -1.
function dimension_via_molien_series(::Type{T}, R::InvRing, d::Int) where T <: IntegerUnion
  if !ismolien_series_implemented(R)
    return -1
  end

  Qt, t = PowerSeriesRing(QQ, d + 1, "t")
  F = molien_series(R)
  k = coeff(numerator(F)(t)*inv(denominator(F)(t)), d)
  @assert isintegral(k)
  return T(numerator(k))::T
end

@doc Markdown.doc"""
     iterate_basis(IR::InvRing, d::Int, algo::Symbol = :default)

Given an invariant ring `IR` and an integer `d`, return an iterator over a basis
for the invariants in degree `d`.
The used algorithm can be specified using the optional argument `algo`. Possible
values are `:reynolds` which uses the reynolds operator to construct the basis
(only available in the non-modular case) and `:linear_algebra` which uses plain
linear algebra. With the default value `:default` the heuristically best algorithm
is selected.

When using the reynolds operator the basis is constructed element-by-element.
With linear algebra this is not possible and the whole basis will be constructed
directly when calling the function.

See also [`basis`](@ref).

# Examples
```
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> B = iterate_basis(IR, 6)
Iterator over a basis of the component of degree 6 of
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> collect(B)
4-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]^2*x[2]^2*x[3]^2
 x[1]^4*x[2]*x[3] + x[1]*x[2]^4*x[3] + x[1]*x[2]*x[3]^4
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3 + x[2]^3*x[3]^3
 x[1]^6 + x[2]^6 + x[3]^6

julia> M = matrix(GF(3), [0 1 0; -1 0 0; 0 0 -1])
[0   1   0]
[2   0   0]
[0   0   2]

julia> G = MatrixGroup(3, GF(3), [M])
Matrix group of degree 3 over Galois field with characteristic 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Galois field with characteristic 3
with generators
gfp_mat[[0 1 0; 2 0 0; 0 0 2]]

julia> B = iterate_basis(IR, 2)
Iterator over a basis of the component of degree 2 of
Invariant ring of
Matrix group of degree 3 over Galois field with characteristic 3
with generators
gfp_mat[[0 1 0; 2 0 0; 0 0 2]]

julia> collect(B)
2-element Vector{MPolyElem_dec{gfp_elem, gfp_mpoly}}:
 x[1]^2 + x[2]^2
 x[3]^2
```
"""
function iterate_basis(R::InvRing, d::Int, algo::Symbol = :default)
  @assert d >= 0 "Degree must be non-negativ"

  if algo == :default
    if ismodular(R)
      algo = :linear_algebra
    else
      # Use the estimate in KS99, Section 17.2
      # We use the "worst case" estimate, so 2d|G|/s instead of sqrt(2d|G|/s)
      # for the reynolds operator since we have to assume that the user really
      # wants to iterate the whole basis.
      # "Experience" showed that one should drop the 2 in the bound.
      # It probably also depends on the type of the field etc., so one could
      # fine-tune here...
      # But this should be a good heuristic anyways and the users can
      # always choose for themselves :)
      s = length(action(R))
      g = order(Int, group(R))
      n = degree(group(R))
      k = binomial(n + d - 1, n - 1)
      if k > d*g/s
        algo = :reynolds
      else
        algo = :linear_algebra
      end
    end
  end

  if algo == :reynolds
    return iterate_basis_reynolds(R, d)
  elseif algo == :linear_algebra
    return iterate_basis_linear_algebra(R, d)
  else
    error("Unsupported argument :$(algo) for algo.")
  end
end

function iterate_basis_reynolds(R::InvRing, d::Int)
  @assert !ismodular(R)
  @assert d >= 0 "Degree must be non-negativ"

  monomials = all_monomials(polynomial_ring(R), d)

  k = dimension_via_molien_series(Int, R, d)
  @assert k != -1

  N = zero_matrix(base_ring(polynomial_ring(R)), 0, 0)
  return InvRingBasisIterator{typeof(R), typeof(monomials), eltype(monomials), typeof(N)}(R, d, k, true, monomials, Vector{eltype(monomials)}(), N)
end

# Sadly, we can't really do much iteratively here.
function iterate_basis_linear_algebra(IR::InvRing, d::Int)
  @assert d >= 0 "Degree must be non-negativ"

  R = polynomial_ring(IR)

  k = dimension_via_molien_series(Int, IR, d)
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
    M = zero_matrix(K, BI.dim, length(g))
    pivot_rows = zeros(Int, length(g))
    g = inv(leading_coefficient(g))*g
    i = 1
    for (c, m) in zip(coefficients(g), monomials(g))
      monomial_to_basis[m] = i
      M[1, i] = c
      i += 1
    end
    pivot_rows[1] = 1
    @assert M[1, 1] == 1
    v = zeros(K, ncols(M))
    return g, (M, monomial_to_basis, state, pivot_rows, 1, v)
  end
end

function iterate_reynolds(BI::InvRingBasisIterator, state)
  @assert BI.reynolds

  M = state[1]
  r = state[5]
  if r == BI.dim
    return nothing
  end

  pivot_rows = state[4]
  v = state[6]

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

    v = zero!.(v)
    new_cols = 0
    for (c, m) in zip(coefficients(g), monomials(g))
      if !haskey(monomial_to_basis, m)
        monomial_to_basis[m] = length(v) + 1
        new_cols += 1
        push!(v, c)
      else
        col = monomial_to_basis[m]
        v[col] = c
      end
    end
    if !iszero(new_cols)
      append!(pivot_rows, zeros(Int, new_cols))
      M = hcat(M, zero_matrix(K, nrows(M), new_cols))
    end
    if Hecke._add_row_to_rref!(M, v, pivot_rows, r + 1)
      return inv(leading_coefficient(g))*g, (M, monomial_to_basis, monomial_state, pivot_rows, r + 1, v)
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
  # Have to (should...) divide by the leading coefficient again:
  # The matrix was in echelon form, but the columns were not necessarily sorted
  # w.r.t. the monomial ordering.
  return inv(leading_coefficient(f))*f, 2
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
  return inv(leading_coefficient(f))*f, state + 1
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

################################################################################
#
# MSetPartitions
#
################################################################################

iterate_partitions(M::MSet) = MSetPartitions(M)

Base.eltype(MSP::MSetPartitions{T}) where T = Vector{MSet{T}}

function Base.iterate(MSP::MSetPartitions)
  if isempty(MSP.M)
    return nothing
  end

  return [ MSP.M ], MSetPartitionsState(MSP)
end

# This is basically Knu11, p. 429, Algorithm 7.2.1.5M
# M2 - 6 in the  comments correspond to the steps in the pseudocode
function Base.iterate(MSP::MSetPartitions{T}, state::MSetPartitionsState) where T
  c = state.c
  u = state.u
  v = state.v
  f = state.f
  a = state.a
  b = state.b
  l = state.l

  m = length(MSP.num_to_key)
  n = length(MSP.M)

  # M5
  j = b - 1
  while iszero(v[j])
    j -= 1
  end
  while j == a && v[j] == 1
    # M6
    if l == 1
      return nothing
    end
    l -= 1
    b = a
    a = f[l]

    j = b - 1
    while iszero(v[j])
      j -= 1
    end
  end
  v[j] = v[j] - 1
  for k = j + 1:b - 1
    v[k] = u[k]
  end

  # M2
  range_increased = true
  while range_increased
    k = b
    range_increased = false
    v_changed = false
    for j = a:b - 1
      u[k] = u[j] - v[j]
      if iszero(u[k])
        v_changed = true
        continue
      end
      c[k] = c[j]
      if v_changed
        v[k] = u[k]
      else
        v[k] = min(v[j], u[k])
        v_changed = (u[k] < v[j])
      end
      k += 1
    end

    # M3
    if k > b
      range_increased = true
      a = b
      b = k
      l += 1
      f[l + 1] = b
    end
  end

  # M4
  part = Vector{typeof(MSP.M)}()
  for j = 1:l
    N = MSet{T}()
    for k = f[j]:f[j + 1] - 1
      if iszero(v[k])
        continue
      end
      N.dict[MSP.num_to_key[c[k]]] = v[k]
    end
    push!(part, N)
  end

  state.c = c
  state.u = u
  state.v = v
  state.f = f
  state.a = a
  state.b = b
  state.l = l

  return part, state
end
