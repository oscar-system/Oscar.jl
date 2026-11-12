# ---------------------------------------------------------------------------
# Helper functions shared by the primary decomposition algorithms: independent
# sets, extension to the rational function field K(u) and contraction back,
# saturation by leading coefficients and ideal comparison.
# ---------------------------------------------------------------------------

# The terms of `f` as pairs `(coefficient, exponent vector)`. The iterators
# `coefficients` and `exponent_vectors` are not reliably aligned for all
# polynomial types, so the terms are accessed by index.
_pd_terms(f::MPolyRingElem) = ((coeff(f, i), exponent_vector(f, i)) for i in 1:length(f))

# Equality of ideals.
_pd_ideals_equal(I::MPolyIdeal, J::MPolyIdeal) = is_subset(I, J) && is_subset(J, I)

# Whether `I` is the unit ideal.
_pd_is_unit(I::MPolyIdeal) = is_one(I)

# Remove zero generators and make sure that the ideal has a Groebner basis cached.
function _pd_clean(I::MPolyIdeal)
  R = base_ring(I)
  g = filter(!is_zero, gens(I))
  isempty(g) && return ideal(R, elem_type(R)[])
  return ideal(R, g)
end

# The exponent vectors of the leading monomials of a Groebner basis of `I`
# with respect to the degree reverse lexicographical ordering.
function _pd_leading_exponents(I::MPolyIdeal)
  R = base_ring(I)
  o = degrevlex(R)
  G = groebner_basis(I; ordering=o)
  return [exponent_vector(leading_monomial(g; ordering=o), 1) for g in G]
end

# A maximal strongly independent set of `I` of size `dim(I)` as a sorted vector
# of variable indices. Variables are preferred from the back, as in Singular.
function _pd_independent_set(I::MPolyIdeal)
  R = base_ring(I)
  n = ngens(R)
  d = dim(I)
  d <= 0 && return Int[]
  leads = _pd_leading_exponents(I)
  supports = [Set(findall(!iszero, e)) for e in leads]
  for c in Iterators.reverse(collect(combinations(1:n, d)))
    cs = Set(c)
    any(s -> issubset(s, cs), supports) && continue
    return sort(c)
  end
  error("no independent set of size $d found")
end

# Product ordering with the variables in `X` before the variables in `U`.
function _pd_block_ordering(R::MPolyRing, U::Vector{Int})
  X = setdiff(1:ngens(R), U)
  isempty(U) && return degrevlex(R)
  isempty(X) && return degrevlex(R)
  return degrevlex(gens(R)[X]) * degrevlex(gens(R)[U])
end

# Leading coefficient of `g` as a polynomial in the variables `U` with respect to
# the leading monomial in the variables outside of `U`. The result lives in `R`.
function _pd_lead_coeff_in(g::MPolyRingElem, U::Vector{Int}, o)
  R = parent(g)
  n = ngens(R)
  X = setdiff(1:n, U)
  lm = exponent_vector(leading_monomial(g; ordering=o), 1)
  ctx = MPolyBuildCtx(R)
  for (c, e) in _pd_terms(g)
    all(e[i] == lm[i] for i in X) || continue
    e2 = copy(e)
    for i in X
      e2[i] = 0
    end
    push_term!(ctx, c, e2)
  end
  return finish(ctx)
end

# The distinct leading coefficients (in `K[U]`) that are not units of a Groebner
# basis of `I` with respect to the block ordering with the variables `U` last.
# Their lcm is the polynomial `h` with `I : h^oo` the contraction of the extension
# of `I` to `K(U)[x']`; it is kept as a list since factoring the coefficients
# separately is much cheaper than factoring their lcm.
function _pd_lead_coeffs(I::MPolyIdeal, U::Vector{Int})
  R = base_ring(I)
  isempty(U) && return elem_type(R)[]
  o = _pd_block_ordering(R, U)
  G = groebner_basis(I; ordering=o)
  cs = elem_type(R)[]
  for g in G
    c = _pd_lead_coeff_in(g, U, o)
    is_unit(c) || c in cs || push!(cs, c)
  end
  return cs
end

# The distinct irreducible factors of the polynomials `cs`, sorted by degree.
function _pd_irreducible_factors(cs::Vector)
  facs = eltype(cs)[]
  for c in cs, (f, _) in factor(c)
    f in facs || push!(facs, f)
  end
  return sort!(facs; by=total_degree)
end

# Saturation of `I` by the product of the polynomials `facs`, computed one factor
# at a time, which is much cheaper than saturating by the product.
function _pd_saturate(I::MPolyIdeal, facs::Vector)
  R = base_ring(I)
  S = I
  for g in facs
    S = saturation(S, ideal(R, [g]))
  end
  return S
end

# Saturation of `I` by the lcm `h` of the polynomials `cs` together with a
# polynomial `f` (a minimal divisor of the squarefree part of `h` with the same
# saturation) and the least `k` such that `I : f^k` is the saturation, so that
# `I = (I : h^oo) ∩ (I + <f^k>)`; see `minsat` in Singular's `primdec.lib`.
# Returns `(saturation, factors, f, k)` where `factors` are the irreducible
# factors of `f`. If `I` is already saturated, `k` is `0` and `factors` is empty.
# With `minimal=false` only the saturation is computed and `k` is `0`.
function _pd_split_by(I::MPolyIdeal, cs::Vector; minimal::Bool=true)
  R = base_ring(I)
  isempty(cs) && return I, elem_type(R)[], one(R), 0
  facs = _pd_irreducible_factors(cs)
  S = _pd_saturate(I, facs)
  if !minimal || is_subset(S, I)
    return S, minimal ? elem_type(R)[] : facs, prod(facs), 0
  end
  # drop factors that are not needed to obtain the saturation
  for g in reverse(facs)
    length(facs) == 1 && break
    rest = filter(!=(g), facs)
    is_subset(S, _pd_saturate(I, rest)) && (facs = rest)
  end
  f = prod(facs)
  # the least `k` with `I : f^k = S`
  k = 0
  Ik = I
  while !is_subset(S, Ik)
    Ik = quotient(Ik, ideal(R, [f]))
    k += 1
  end
  return S, facs, f, k
end

# ---------------------------------------------------------------------------
# Passing to K(u)[x] and back
# ---------------------------------------------------------------------------

# Data of the extension of the base field by the independent variables `U`.
struct _PDExtension
  R::MPolyRing        # the original ring K[x_1, ..., x_n]
  U::Vector{Int}      # the indices of the independent variables
  X::Vector{Int}      # the indices of the other variables
  Qu::MPolyRing       # K[u]
  F::Field            # K(u)
  Rx::MPolyRing       # K(u)[x']
end

function _pd_extension(R::MPolyRing, U::Vector{Int})
  K = coefficient_ring(R)
  X = setdiff(1:ngens(R), U)
  syms = symbols(R)
  Qu, _ = polynomial_ring(K, syms[U]; cached=false)
  F = fraction_field(Qu)
  Rx, _ = polynomial_ring(F, syms[X]; cached=false)
  return _PDExtension(R, U, X, Qu, F, Rx)
end

# Map a polynomial of `K[x]` to `K(u)[x']`.
function _pd_to_ext(E::_PDExtension, f::MPolyRingElem)
  isempty(E.U) && return f
  ug = gens(E.Qu)
  ctx = MPolyBuildCtx(E.Rx)
  for (c, e) in _pd_terms(f)
    m = one(E.Qu)
    for (j, i) in enumerate(E.U)
      e[i] > 0 && (m *= ug[j]^e[i])
    end
    push_term!(ctx, E.F(c * m), e[E.X])
  end
  return finish(ctx)
end

# Map a polynomial of `K(u)[x']` to a polynomial of `K[x]` by clearing
# denominators.
function _pd_from_ext(E::_PDExtension, g::MPolyRingElem)
  isempty(E.U) && return g
  R = E.R
  n = ngens(R)
  den = one(E.Qu)
  for c in coefficients(g)
    den = lcm(den, denominator(c))
  end
  ctx = MPolyBuildCtx(R)
  for (c, e) in _pd_terms(g)
    q = numerator(c) * divexact(den, denominator(c))
    for (cu, eu) in _pd_terms(q)
      ee = zeros(Int, n)
      ee[E.X] = e
      ee[E.U] = eu
      push_term!(ctx, cu, ee)
    end
  end
  return finish(ctx)
end

# The contraction `J K(u)[x'] ∩ K[x]` of an ideal `J` of `K(u)[x']` given by
# generators `gs` of `J`.
function _pd_contract(E::_PDExtension, gs::Vector)
  R = E.R
  L = ideal(R, filter(!is_zero, [_pd_from_ext(E, g) for g in gs]))
  isempty(E.U) && return L
  # when many saturations are needed, linear algebra gives a much better start
  boost = function (S)
    isempty(gs) && return nothing
    Jx = ideal(E.Rx, filter(!is_zero, gs))
    leads = _pd_leading_exponents(Jx)
    all(i -> any(l -> l[i] > 0 && count(!iszero, l) == 1, leads), 1:length(E.X)) || return nothing
    new = _pd_contraction_generators(E, _pd_zero_data(Jx), maximum(total_degree, gens(S)))
    return new === nothing ? nothing : S + ideal(R, new)
  end
  return _pd_saturate_leads(L, E.U, boost)
end

# The exponent vectors of the monomials in `n` variables of degree at most `D`,
# sorted by decreasing degree.
function _pd_monomials_upto(n::Int, D::Int)
  res = Vector{Vector{Int}}()
  cur = zeros(Int, n)
  function rec(i, left)
    if i > n
      push!(res, copy(cur))
      return
    end
    for k in 0:left
      cur[i] = k
      rec(i + 1, left - k)
    end
    cur[i] = 0
  end
  rec(1, D)
  return sort!(res; by=e -> (-sum(e), e))
end

# Whether the polynomial `f` of `K[x]` lies in the ideal of `K(u)[x']` with the
# zero-dimensional data `Z`, that is whether its normal form vanishes.
function _pd_in_extension(E::_PDExtension, Z, f::MPolyRingElem)
  ug = gens(E.Qu)
  acc = Any[zero(E.F) for _ in Z.basis]
  for (c, e) in _pd_terms(f)
    m = c * prod(ug[j]^e[i] for (j, i) in enumerate(E.U); init=one(E.Qu))
    r = _pd_monomial_row(Z, e[E.X])
    for i in eachindex(acc)
      iszero(r[i]) || (acc[i] += E.F(m) * r[i])
    end
  end
  return all(iszero, acc)
end

_pd_first_nonzero(A, i) = findfirst(j -> !iszero(A[i, j]), 1:ncols(A))

# Generators of the contraction `J K(u)[x'] ∩ K[x]` of degree at most `D`, for
# the zero-dimensional ideal `J` with data `Z`, found by linear algebra: `f` is
# in `J` iff the combination of the normal forms of its monomials vanishes, which
# is tested at random points of `K^U`. Generators are selected degree by degree,
# modulo the multiples of the previous ones, and verified exactly. Returns
# `nothing` if the field is too small for random points.
function _pd_contraction_generators(E::_PDExtension, Z, D::Int)
  R = E.R
  K = coefficient_ring(R)
  d = length(Z.basis)
  (K isa FinField && order(K) < 1000) && return nothing
  cols = _pd_monomials_upto(ngens(R), D)
  nc = length(cols)
  rows = Dict{Vector{Int},Vector{Any}}()
  for e in cols
    a = e[E.X]
    haskey(rows, a) || (rows[a] = _pd_monomial_row(Z, a))
  end
  npts = cld(nc, d) + 10
  M = zero_matrix(K, npts * d, nc)
  cnt = 0
  while cnt < npts
    p = [K(rand(1:1000)) for _ in E.U]
    vals = Dict{Vector{Int},Vector{Any}}()
    ok = true
    for (a, r) in rows
      v = Any[]
      for c in r
        den = evaluate(denominator(c), p)
        iszero(den) && (ok = false; break)
        push!(v, evaluate(numerator(c), p) // den)
      end
      ok || break
      vals[a] = v
    end
    ok || continue
    cnt += 1
    for (j, e) in enumerate(cols)
      w = prod(p[k]^e[i] for (k, i) in enumerate(E.U); init=one(K))
      v = vals[e[E.X]]
      for i in 1:d
        M[(cnt - 1) * d + i, j] = v[i] * w
      end
    end
  end
  gens_found = elem_type(R)[]
  xs = gens(R)
  for deg in 1:D
    sub = findall(e -> sum(e) <= deg, cols)
    nl, N = nullspace(M[:, sub])
    nl == 0 && continue
    subidx = Dict(cols[j] => k for (k, j) in enumerate(sub))
    # multiples of the generators found so far, as rows in the monomials of degree <= deg
    mult = Vector{Vector{elem_type(K)}}()
    for g in gens_found
      dg = total_degree(g)
      dg > deg - 1 && continue
      for m in _pd_monomials_upto(ngens(R), deg - dg)
        row = zeros(K, length(sub))
        for (c, e) in _pd_terms(g)
          row[subidx[e + m]] = c
        end
        push!(mult, row)
      end
    end
    nm = length(mult)
    A = zero_matrix(K, nm + nl, length(sub))
    for (i, row) in enumerate(mult), j in eachindex(row)
      A[i, j] = row[j]
    end
    for i in 1:nl, j in eachindex(sub)
      A[nm + i, j] = N[j, i]
    end
    r1, B1 = rref(A[1:nm, :])
    piv1 = Set{Int}(_pd_first_nonzero(B1, i) for i in 1:r1)
    r2, A2 = rref(A)
    for i in 1:r2
      j = _pd_first_nonzero(A2, i)
      j in piv1 && continue
      g = zero(R)
      for jj in 1:ncols(A2)
        c = A2[i, jj]
        iszero(c) && continue
        g += c * prod(xs[k]^cols[sub[jj]][k] for k in 1:ngens(R))
      end
      _pd_in_extension(E, Z, g) && push!(gens_found, g)
    end
  end
  return gens_found
end

# Seconds of saturation after which the contraction starts from generators found
# by linear algebra.
const _PD_BOOST_AFTER = 5.0

# The contraction `I·K(U)[x'] ∩ K[x]`, that is the saturation of `I` by the leading
# coefficients of a Groebner basis. Any ideal between `I` and the contraction has
# the same contraction, so after each saturation that changes the ideal the
# leading coefficients are computed again, now from a Groebner basis that is
# usually much smaller. A factor that has been saturated stays saturated.
function _pd_saturate_leads(I::MPolyIdeal, U::Vector{Int}, boost=nothing)
  R = base_ring(I)
  S = I
  done = elem_type(R)[]
  t0 = time()
  while true
    facs = filter(g -> !(g in done), _pd_irreducible_factors(_pd_lead_coeffs(S, U)))
    isempty(facs) && return S
    for g in facs
      push!(done, g)
      if boost !== nothing && time() - t0 > _PD_BOOST_AFTER
        # `boost` returns an ideal between `S` and the contraction, or nothing
        B = boost(S)
        boost = nothing
        if B !== nothing
          S = B
          break
        end
      end
      S2 = saturation(S, ideal(R, [g]))
      if !is_subset(S2, S)
        S = S2
        break
      end
    end
  end
end

# Extend an ideal of `K[x]` to `K(u)[x']`.
function _pd_extend_ideal(E::_PDExtension, I::MPolyIdeal)
  isempty(E.U) && return I
  return ideal(E.Rx, filter(!is_zero, [_pd_to_ext(E, g) for g in gens(I)]))
end

# Remove from `L` the ideals that contain a different ideal of `L`, as well as
# duplicates; used for minimal primes.
function _pd_minimal_elements(L::Vector{<:MPolyIdeal})
  res = eltype(L)[]
  for (i, P) in enumerate(L)
    keep = true
    for (j, Q) in enumerate(L)
      i == j && continue
      if is_subset(Q, P)
        if !is_subset(P, Q) || j < i
          keep = false
          break
        end
      end
    end
    keep && push!(res, P)
  end
  return res
end

# The product of the distinct irreducible factors of `h`.
_pd_squarefree_multivariate(h::MPolyRingElem) = is_unit(h) ? h : prod(f for (f, _) in factor(h))
