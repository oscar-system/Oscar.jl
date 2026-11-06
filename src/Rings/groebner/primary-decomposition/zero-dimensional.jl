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
  return _PDZeroData(I, o, basis, index, Any[nothing for _ in 1:ngens(R)])
end

# Multiplication matrix of the variable `j` on `R/I` in the monomial basis
# (acting on row vectors).
function _pd_mult_matrix(Z::_PDZeroData, j::Int)
  Z.mats[j] === nothing || return Z.mats[j]
  R = base_ring(Z.I)
  F = coefficient_ring(R)
  d = length(Z.basis)
  x = gens(R)[j]
  mons = [_pd_monomial(R, e) for e in Z.basis]
  nfs = normal_form([x * m for m in mons], Z.I; ordering=Z.o)
  M = zero_matrix(F, d, d)
  for (k, f) in enumerate(nfs)
    for (c, e) in _pd_terms(f)
      M[k, Z.index[e]] = c
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

# Evaluate the univariate polynomial `f` at the polynomial `p` modulo `I`.
function _pd_evaluate_mod(Z::_PDZeroData, f, p::MPolyRingElem)
  R = parent(p)
  cs = collect(coefficients(f))
  h = zero(R)
  for k in length(cs):-1:1
    h = normal_form(h * p + R(cs[k]), Z.I; ordering=Z.o)
  end
  return h
end

# Radical of a zero-dimensional ideal over a field of characteristic zero or a
# finite field (Seidenberg's lemma).
function _pd_zero_radical(I::MPolyIdeal)
  R = base_ring(I)
  F = coefficient_ring(R)
  Z = _pd_zero_data(I)
  isempty(Z.basis) && return I
  T, _ = polynomial_ring(F, :T; cached=false)
  gs = copy(gens(I))
  for j in 1:ngens(R)
    m = _pd_min_poly(T, _pd_mult_matrix(Z, j))
    degree(m) <= 1 && continue
    g = _pd_squarefree_part(m)
    degree(g) == degree(m) && continue
    push!(gs, evaluate(g, gens(R)[j]))
  end
  return ideal(R, gs)
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
  if characteristic(F) != 0 && F isa FracField
    # coefficients that are constants may not suffice over small fields
    u = gen(base_ring(F), 1)
    p += F(u^(tries % 3 + 1)) * rand(gens(R))
  end
  return p
end

# Characteristic polynomial of the multiplication by `p` on `R/I`.
function _pd_charpoly(Z::_PDZeroData, T, p::MPolyRingElem)
  R = base_ring(Z.I)
  F = coefficient_ring(R)
  d = length(Z.basis)
  # build the multiplication matrix of `p` from the matrices of the variables;
  # constants and general monomials are handled by reducing
  mons = [_pd_monomial(R, e) for e in Z.basis]
  nfs = normal_form([p * m for m in mons], Z.I; ordering=Z.o)
  M = zero_matrix(F, d, d)
  for (k, f) in enumerate(nfs)
    for (c, e) in _pd_terms(f)
      M[k, Z.index[e]] = c
    end
  end
  return charpoly(T, M)
end

@doc raw"""
    _pd_zero_dimensional_decomposition(I::MPolyIdeal) -> Vector{Tuple{MPolyIdeal, MPolyIdeal}}

For a zero-dimensional ideal `I` return the list of pairs `(Q, P)` of the
primary components `Q` of `I` and their associated primes `P`. The base ring
of `I` is a polynomial ring over a field in which polynomials in one variable
can be factored.
"""
function _pd_zero_dimensional_decomposition(I::MPolyIdeal)
  R = base_ring(I)
  F = coefficient_ring(R)
  I = _pd_clean(I)
  is_one(I) && return Tuple{typeof(I),typeof(I)}[]
  J = _pd_zero_radical(I)
  ZJ = _pd_zero_data(J)
  ZI = _pd_zero_data(I)
  T, _ = polynomial_ring(F, :T; cached=false)
  tries = 0
  while true
    tries += 1
    tries > 60 && error("could not find a separating linear form (is the field too small?)")
    p = _pd_random_linear_form(R, tries)
    chJ = _pd_charpoly(ZJ, T, p)
    fa = collect(factor(chJ))
    any(((_, e),) -> e != 1, fa) && continue
    if length(fa) == 1
      return [(I, J)]
    end
    chI = _pd_charpoly(ZI, T, p)
    faI = [(divexact(g, leading_coefficient(g)), e) for (g, e) in factor(chI)]
    res = Tuple{typeof(I),typeof(I)}[]
    ok = true
    for (g0, _) in fa
      g = divexact(g0, leading_coefficient(g0))
      idx = findfirst(((h, _),) -> h == g, faI)
      idx === nothing && (ok = false; break)
      P = ideal(R, vcat(gens(J), [_pd_evaluate_mod(ZJ, g, p)]))
      Q = ideal(R, vcat(gens(I), [_pd_evaluate_mod(ZI, g^faI[idx][2], p)]))
      push!(res, (Q, P))
    end
    ok || continue
    return res
  end
end
