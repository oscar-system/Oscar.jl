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

# The lcm of the leading coefficients (in `K[U]`) of a Groebner basis of `I` with
# respect to the block ordering with the variables `U` last.
function _pd_lead_lcm(I::MPolyIdeal, U::Vector{Int})
  R = base_ring(I)
  isempty(U) && return one(R)
  o = _pd_block_ordering(R, U)
  G = groebner_basis(I; ordering=o)
  h = one(R)
  for g in G
    c = _pd_lead_coeff_in(g, U, o)
    is_unit(c) && continue
    h = lcm(h, c)
  end
  return h
end

# Saturation of `I` by `h` together with a polynomial `f` (a minimal divisor of
# the squarefree part of `h` with the same saturation) and the least `k` such that
# `I : f^k` is the saturation, so that `I = (I : h^oo) ∩ (I + <f^k>)`; see
# `minsat` in Singular's `primdec.lib`. Returns `(saturation, factors, f, k)`
# where `factors` are the irreducible factors of `f`. With `minimal=false` the
# search for a minimal divisor is skipped, so `f` is the squarefree part of `h`.
function _pd_split_by(I::MPolyIdeal, h::MPolyRingElem; minimal::Bool=true)
  R = base_ring(I)
  is_unit(h) && return I, typeof(h)[], one(R), 0
  facs = sort!([f for (f, _) in factor(h)]; by=total_degree)
  f = prod(facs)
  S, k = saturation_with_index(I, ideal(R, [f]))
  # drop factors that are not needed to obtain the saturation
  if minimal && length(facs) > 1
    for g in reverse(facs)
      length(facs) == 1 && break
      rest = filter(!=(g), facs)
      S2 = saturation(I, ideal(R, [prod(rest)]))
      if is_subset(S, S2)
        facs = rest
      end
    end
    f = prod(facs)
    S, k = saturation_with_index(I, ideal(R, [f]))
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
  h = _pd_lead_lcm(L, E.U)
  S, _, _, _ = _pd_split_by(L, h; minimal=false)
  return S
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
