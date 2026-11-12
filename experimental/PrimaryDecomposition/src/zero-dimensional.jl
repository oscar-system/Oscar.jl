# ---------------------------------------------------------------------------
# Decomposition of zero-dimensional ideals over a field `F` (which is `K` or a
# rational function field `K(u)`). The primes are found by factoring the
# characteristic polynomial of a generic linear form acting on `R/I`, see the
# approach of Monico used by `zerodec` in Singular's `primdec.lib`.
# ---------------------------------------------------------------------------

struct _PDZeroData{T<:MPolyRingElem}
  I::MPolyIdeal{T}
  o::MonomialOrdering
  basis::Vector{Vector{Int}}           # exponent vectors of the standard monomials
  index::Dict{Vector{Int},Int}
  mats::Vector{Any}                    # lazily computed multiplication matrices
  reducers::Vector{Tuple{Vector{Int},Vector{Tuple{Any,Vector{Int}}}}}  # leading exponent, -tail/lc of each GB element
  rows::Dict{Vector{Int},Vector{Any}}  # normal forms of monomials as coordinate vectors
end

# The standard monomials of a zero-dimensional ideal given by the exponent
# vectors `leads` of the leading monomials of a Groebner basis.
function _pd_standard_monomials(leads::Vector{Vector{Int}}, n::Int)
  res = Vector{Int}[]
  isdiv(e) = any(l -> all(l[i] <= e[i] for i in 1:n), leads)
  stack = [zeros(Int, n)]
  seen = Set{Vector{Int}}(stack)
  while !isempty(stack)
    e = pop!(stack)
    isdiv(e) && continue
    push!(res, e)
    for i in 1:n
      e2 = copy(e)
      e2[i] += 1
      if !(e2 in seen)
        push!(seen, e2)
        push!(stack, e2)
      end
    end
  end
  return res
end

function _pd_zero_data(I::MPolyIdeal)
  R = base_ring(I)
  o = degrevlex(R)
  G = groebner_basis(I; ordering=o)
  leads = [exponent_vector(leading_monomial(g; ordering=o), 1) for g in G]
  basis = _pd_standard_monomials(leads, ngens(R))
  sort!(basis; by=e -> (sum(e), e))
  index = Dict(e => i for (i, e) in enumerate(basis))
  reducers = map(zip(G, leads)) do (g, l)
    lc = leading_coefficient(g; ordering=o)
    tail = Tuple{Any,Vector{Int}}[(-c // lc, e) for (c, e) in _pd_terms(g) if e != l]
    return (l, tail)
  end
  return _PDZeroData(I, o, basis, index, Any[nothing for _ in 1:ngens(R)], reducers,
                     Dict{Vector{Int},Vector{Any}}())
end

# The normal form of the monomial with exponent vector `e` as a coordinate vector
# in the standard monomials, obtained by reducing with the Groebner basis inside
# the quotient. This keeps the coefficients far smaller than a normal form of
# polynomials over a rational function field.
function _pd_monomial_row(Z::_PDZeroData, e::Vector{Int})
  haskey(Z.rows, e) && return Z.rows[e]
  F = coefficient_ring(base_ring(Z.I))
  res = Any[zero(F) for _ in Z.basis]
  k = get(Z.index, e, 0)
  if k != 0
    res[k] = one(F)
  else
    for (l, tail) in Z.reducers
      all(l[i] <= e[i] for i in eachindex(e)) || continue
      t = e - l
      for (c, s) in tail
        r = _pd_monomial_row(Z, s + t)
        for i in eachindex(res)
          iszero(r[i]) || (res[i] += c * r[i])
        end
      end
      break
    end
  end
  Z.rows[e] = res
  return res
end

# Multiplication matrix of the variable `j` on `R/I` in the monomial basis
# (acting on row vectors).
function _pd_mult_matrix(Z::_PDZeroData, j::Int)
  Z.mats[j] === nothing || return Z.mats[j]
  R = base_ring(Z.I)
  F = coefficient_ring(R)
  d = length(Z.basis)
  M = zero_matrix(F, d, d)
  for (k, b) in enumerate(Z.basis)
    e = copy(b)
    e[j] += 1
    r = _pd_monomial_row(Z, e)
    for i in 1:d
      iszero(r[i]) || (M[k, i] = r[i])
    end
  end
  Z.mats[j] = M
  return M
end

_pd_monomial(R::MPolyRing, e::Vector{Int}) = monomial(R, e)

# Whether the squarefree-part computation is meaningful: in positive
# characteristic over a rational function field irreducible polynomials may be
# inseparable, which this implementation does not handle.
function _pd_check_separable(f)
  F = coefficient_ring(parent(f))
  characteristic(F) == 0 && return
  F isa Union{FinField, FracField{<:MPolyRingElem}} || return
  F isa FinField && return
  is_zero(derivative(f)) && throw(NotImplementedError(:primary_decomposition,
    "inseparable polynomial over a rational function field in positive characteristic"))
end

# Product of the distinct irreducible factors of the univariate polynomial `f`
# (its squarefree part) over `F`.
function _pd_squarefree_part(f)
  fa = factor(f)
  res = one(parent(f))
  for (g, _) in fa
    _pd_check_separable(g)
    res *= g
  end
  return divexact(res, leading_coefficient(res))
end

# Minimal polynomial of the matrix `M` (acting on row vectors) with respect to the
# cyclic vector `e_1`, which is the coordinate vector of 1 in the monomial basis.
function _pd_min_poly(T, M::MatElem)
  F = base_ring(M)
  d = nrows(M)
  v = [i == 1 ? one(F) : zero(F) for i in 1:d]
  rows = Vector{Vector{elem_type(F)}}()   # reduced vectors
  combs = Vector{Vector{elem_type(F)}}()  # the combinations producing them
  piv = Int[]
  k = 0
  cur = v
  while true
    w = copy(cur)
    c = [zero(F) for _ in 1:d + 1]
    c[k+1] = one(F)
    for (r, cr, p) in zip(rows, combs, piv)
      if !is_zero(w[p])
        s = w[p]
        for i in 1:d
          w[i] -= s * r[i]
        end
        for i in eachindex(c)
          c[i] -= s * cr[i]
        end
      end
    end
    p = findfirst(!is_zero, w)
    if p === nothing
      # c is the dependency, monic in degree k
      return sum(c[i] * gen(T)^(i - 1) for i in 1:k+1)
    end
    s = inv(w[p])
    w .*= s
    c .*= s
    push!(rows, w); push!(combs, c); push!(piv, p)
    # next Krylov vector: cur * M
    nxt = [zero(F) for _ in 1:d]
    for i in 1:d
      is_zero(cur[i]) && continue
      for j in 1:d
        nxt[j] += cur[i] * M[i, j]
      end
    end
    cur = nxt
    k += 1
  end
end

# The pairs `(j, g)` such that `g(x_j)` for the squarefree part `g` of the minimal
# polynomial of the variable `x_j` on `R/I`; together with `I` they generate the
# radical (Seidenberg's lemma, for a field of characteristic zero or a finite
# field).
function _pd_radical_generators(Z::_PDZeroData, T)
  R = base_ring(Z.I)
  res = Tuple{Int,Any}[]
  for j in 1:ngens(R)
    m = _pd_min_poly(T, _pd_mult_matrix(Z, j))
    degree(m) <= 1 && continue
    g = _pd_squarefree_part(m)
    degree(g) == degree(m) && continue
    push!(res, (j, g))
  end
  return res
end

function _pd_radical_ideal(Z::_PDZeroData, rg)
  R = base_ring(Z.I)
  return ideal(R, vcat(copy(gens(Z.I)), [evaluate(g, gens(R)[j]) for (j, g) in rg]))
end

# Radical of a zero-dimensional ideal over a field of characteristic zero or a
# finite field (Seidenberg's lemma).
function _pd_zero_radical(I::MPolyIdeal)
  Z = _pd_zero_data(I)
  isempty(Z.basis) && return I
  T, _ = polynomial_ring(coefficient_ring(base_ring(I)), :T; cached=false)
  # the small generators from the minimal polynomials make the contraction cheaper
  return _pd_radical_ideal(Z, _pd_radical_generators(Z, T))
end

# The matrix of the polynomial `g` evaluated at the square matrix `M`.
function _pd_matrix_poly(g, M::MatElem)
  F = base_ring(M)
  d = nrows(M)
  cs = collect(coefficients(g))
  H = cs[end] * identity_matrix(F, d)
  for k in length(cs)-1:-1:1
    H = H * M
    for i in 1:d
      H[i, i] += cs[k]
    end
  end
  return H
end

# The nilradical of `R/I` as a subspace of the coordinate space, in reduced row
# echelon form, or `nothing` if `I` is radical. It is spanned by the products of
# the basis with the radical generators, so no Groebner basis of the radical is
# needed.
function _pd_nilradical(Z::_PDZeroData, rg)
  isempty(rg) && return nothing
  F = coefficient_ring(base_ring(Z.I))
  d = length(Z.basis)
  S = zero_matrix(F, d * length(rg), d)
  for (t, (j, g)) in enumerate(rg)
    G = _pd_matrix_poly(g, _pd_mult_matrix(Z, j))
    for a in 1:d, b in 1:d
      S[(t-1)*d+a, b] = G[a, b]
    end
  end
  r, Rm = rref(S)
  return (r, Rm)
end

# The nilradical via the trace form, in the format of `_pd_nilradical` (`nothing`
# if `I` is radical): the kernel of `tr(x^(a+b))` for the standard monomials `a`,
# `b`. This is valid in characteristic zero and in characteristic larger than the
# dimension of `R/I`; otherwise `missing` is returned. All traces are sums of the
# cached normal forms of monomials, so no matrix products are needed.
function _pd_nilradical_trace(Z::_PDZeroData)
  F = coefficient_ring(base_ring(Z.I))
  d = length(Z.basis)
  c = characteristic(F)
  (c != 0 && c <= d) && return missing
  T = zero_matrix(F, d, d)
  for a in 1:d, b in a:d
    s = zero(F)
    for i in 1:d
      s += _pd_monomial_row(Z, Z.basis[a] + Z.basis[b] + Z.basis[i])[i]
    end
    T[a, b] = s
    T[b, a] = s
  end
  nl, K = nullspace(T)
  nl == 0 && return nothing
  r, Rm = rref(transpose(K))
  return (r, Rm)
end

# The radical `J` of `I` (given by its nilradical or generated by the radical
# generators) together with the nilradical of `R/I`.
function _pd_radical_data(Z::_PDZeroData, T)
  N = _pd_nilradical_trace(Z)
  if N === missing
    rg = _pd_radical_generators(Z, T)
    return _pd_radical_ideal(Z, rg), _pd_nilradical(Z, rg)
  end
  N === nothing && return Z.I, nothing
  return _pd_ideal_with_subspace(Z, N[2][1:N[1], :]), N
end

# The matrix of the endomorphism `M` of `R/I` on the quotient by the nilradical
# `N`, in the basis given by the non-pivot coordinates.
function _pd_free_columns(N, d::Int)
  N === nothing && return collect(1:d), Int[]
  r, Rm = N
  piv = [findfirst(k -> !is_zero(Rm[i, k]), 1:d) for i in 1:r]
  return setdiff(1:d, piv), piv
end

# The rows of `G`, a matrix on the quotient by the nilradical, as vectors of
# `R/I` with zero coordinates at the pivot columns of the nilradical.
function _pd_embed_free(G::MatElem, N, d::Int)
  free, _ = _pd_free_columns(N, d)
  S = zero_matrix(base_ring(G), nrows(G), d)
  for a in 1:nrows(G), (b, c) in enumerate(free)
    S[a, c] = G[a, b]
  end
  return S
end

function _pd_quotient_matrix(M::MatElem, N)
  N === nothing && return M
  F = base_ring(M)
  r, Rm = N
  d = nrows(M)
  free, piv = _pd_free_columns(N, d)
  Q = zero_matrix(F, length(free), length(free))
  for (a, c) in enumerate(free)
    v = [M[c, k] for k in 1:d]
    for i in 1:r
      s = v[piv[i]]
      is_zero(s) && continue
      for k in 1:d
        v[k] -= s * Rm[i, k]
      end
    end
    for (b, c2) in enumerate(free)
      Q[a, b] = v[c2]
    end
  end
  return Q
end

# A random linear form in the variables with coefficients in `F`.
function _pd_random_linear_form(R::MPolyRing, tries::Int)
  F = coefficient_ring(R)
  range = 5 * tries + 5
  p = zero(R)
  for x in gens(R)
    c = F(rand(-range:range))
    p += c * x
  end
  if tries > 2 && characteristic(F) != 0 && F isa FracField
    # coefficients that are constants may not suffice over small fields; the
    # extra term makes the matrices larger, so it is only used after a failure
    u = gen(base_ring(F), 1)
    p += F(u^(tries % 3 + 1)) * rand(gens(R))
  end
  return p
end

# The multiplication matrix of `p` on `R/I`.
function _pd_poly_mult_matrix(Z::_PDZeroData, p::MPolyRingElem)
  R = base_ring(Z.I)
  F = coefficient_ring(R)
  d = length(Z.basis)
  M = zero_matrix(F, d, d)
  for (k, b) in enumerate(Z.basis), (c, e) in _pd_terms(p)
    r = _pd_monomial_row(Z, b + e)
    for i in 1:d
      iszero(r[i]) || (M[k, i] += c * r[i])
    end
  end
  return M
end

# Generators of `I` together with the elements of the subspace of `R/I` spanned by
# the rows of `rows` (in the monomial basis), as the combinations of standard
# monomials given by its reduced row echelon form with the largest monomials as
# pivots. This is a much smaller description than a normal form of a single
# polynomial and generates `I` plus the ideal of `R/I` that `rows` span.
function _pd_ideal_with_subspace(Z::_PDZeroData, rows::MatElem)
  R = base_ring(Z.I)
  d = length(Z.basis)
  S = zero_matrix(coefficient_ring(R), nrows(rows), d)
  for a in 1:nrows(rows), b in 1:d
    S[a, b] = rows[a, d + 1 - b]
  end
  r, Rm = rref(S)
  mons = [_pd_monomial(R, Z.basis[d + 1 - c]) for c in 1:d]
  polys = [sum(Rm[i, c] * mons[c] for c in 1:d if !is_zero(Rm[i, c])) for i in 1:r]
  return ideal(R, vcat(copy(gens(Z.I)), polys))
end

# The rows of the nilradical `N` (see `_pd_nilradical`) on top of the rows of `M`.
function _pd_stack_nilradical(N, M::MatElem)
  N === nothing && return M
  r, Rm = N
  d = ncols(M)
  S = zero_matrix(base_ring(M), r + nrows(M), d)
  for a in 1:r, b in 1:d
    S[a, b] = Rm[a, b]
  end
  for a in 1:nrows(M), b in 1:d
    S[r + a, b] = M[a, b]
  end
  return S
end

@doc raw"""
    _pd_zero_dimensional_decomposition(I::MPolyIdeal) -> Vector{Tuple{MPolyIdeal, MPolyIdeal}}

For a zero-dimensional ideal `I` return the list of pairs `(Q, P)` of the
primary components `Q` of `I` and their associated primes `P`. The base ring
of `I` is a polynomial ring over a field in which polynomials in one variable
can be factored.
"""
function _pd_zero_dimensional_decomposition(I::MPolyIdeal; primary::Bool=true)
  R = base_ring(I)
  F = coefficient_ring(R)
  I = _pd_clean(I)
  is_one(I) && return Tuple{typeof(I),typeof(I)}[]
  Z = _pd_zero_data(I)
  d = length(Z.basis)
  T, _ = polynomial_ring(F, :T; cached=false)
  J, N = _pd_radical_data(Z, T)
  tries = 0
  while true
    tries += 1
    tries > 60 && error("could not find a separating linear form (is the field too small?)")
    p = _pd_random_linear_form(R, tries)
    Mp = _pd_poly_mult_matrix(Z, p)
    chJ = charpoly(T, _pd_quotient_matrix(Mp, N))
    fa = collect(factor(chJ))
    any(((_, e),) -> e != 1, fa) && continue
    if length(fa) == 1
      return [(primary ? I : J, J)]
    end
    Qm = _pd_quotient_matrix(Mp, N)
    if primary
      chI = charpoly(T, Mp)
      faI = [(divexact(g, leading_coefficient(g)), e) for (g, e) in factor(chI)]
    end
    res = Tuple{typeof(I),typeof(I)}[]
    ok = true
    for (g0, _) in fa
      g = divexact(g0, leading_coefficient(g0))
      # the prime is the nilradical together with the image of `g(p)` on the quotient
      G = _pd_embed_free(_pd_matrix_poly(g, Qm), N, d)
      P = _pd_ideal_with_subspace(Z, _pd_stack_nilradical(N, G))
      if primary
        idx = findfirst(((h, _),) -> h == g, faI)
        idx === nothing && (ok = false; break)
        Q = _pd_ideal_with_subspace(Z, _pd_matrix_poly(g^faI[idx][2], Mp))
      else
        Q = P
      end
      push!(res, (Q, P))
    end
    ok || continue
    return res
  end
end
