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
  oneR = elem_type(base_ring(R))[ one(base_ring(R)) ]
  exps = Vector{Int}[ zeros(Int, nvars(R)) ]

  function new_elt(w::Vector{Int})
    exps[1] = w
    return R(oneR, exps)
  end
  return ( new_elt(w) for w in WeightedIntegerVectors(nvars(R), d))
end

# Iterate over the basis of the degree d component of an invariant ring.
# algo can be either :default, :reynolds or :linear_algebra.
function InvRingBasisIterator(R::InvRing{FldT, GrpT, PolyElemT, PolyRingT, ActionT, SingularActionT}, d::Int, algo::Symbol = :default) where {FldT, GrpT, PolyElemT, PolyRingT, ActionT, SingularActionT}
  @assert d >= 0 "Degree must be non-negativ"

  Qt, t = PowerSeriesRing(QQ, d + 1, "t")
  F = molien_series(R)
  k = coeff(numerator(F)(t)*inv(denominator(F)(t)), d) # Dimension of the d-th component
  @assert isintegral(k)
  k = Int(numerator(k))
  reynolds = false
  if algo == :reynolds
    reynolds = true
  elseif algo == :linear_algebra
    reynolds = false
  elseif algo == :default
    # TODO: Fine tune this: Depending on d and the group order it is better
    # to use "via_linear_algebra" also in the non-modular case.
    if ismodular(R)
      reynolds = false
    else
      reynolds = true
    end
  else
    error("Unsupported argument :$(algo) for algo.")
  end

  monomials = all_monomials(polynomial_ring(R), d)
  return InvRingBasisIterator{FldT, GrpT, PolyElemT, PolyRingT, ActionT, SingularActionT, typeof(monomials)}(R, d, k, reynolds, monomials)
end

Base.eltype(BI::InvRingBasisIterator{FldT, GrpT, PolyElemT}) where {FldT, GrpT, PolyElemT} = PolyElemT

Base.length(BI::InvRingBasisIterator) = BI.dim

function Base.iterate(BI::InvRingBasisIterator)
  if BI.dim == 0
    return nothing
  end

  @assert BI.reynolds "Not implemented (yet)"

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
    return g, InvRingBasisIteratorState{typeof(M), typeof(g), typeof(state)}(M, monomial_to_basis, state)
  end
end

function Base.iterate(BI::InvRingBasisIterator, state::InvRingBasisIteratorState)
  @assert BI.reynolds "Not implemented (yet)"

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
      return g, InvRingBasisIteratorState{typeof(N), typeof(g), typeof(monomial_state)}(N, state.monomial_to_basis, monomial_state)
    end
  end
end

iterate_basis(IR::InvRing, d::Int, algo::Symbol = :default) = InvRingBasisIterator(IR, d, algo)

