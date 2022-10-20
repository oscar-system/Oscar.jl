module Orderings

using Oscar, Markdown
import Oscar: Ring, MPolyRing, MPolyElem, weights, IntegerUnion, base_ring,
       support, matrix
export anti_diagonal, lex, degrevlex, deglex, revlex, negdeglex,
       neglex, negrevlex, negdegrevlex, wdeglex, wdegrevlex,
       negwdeglex, negwdegrevlex, matrix_ordering, monomial_ordering,
       isweighted, is_global, is_local, is_mixed,
       permutation_of_terms, weight_ordering, canonical_matrix,
       MonomialOrdering, ModuleOrdering, singular, opposite_ordering,
       is_elimination_ordering

abstract type AbsOrdering end

abstract type AbsGenOrdering <: AbsOrdering end

abstract type AbsModOrdering <: AbsOrdering end

# lex, deglex, ...
struct SymbOrdering{S} <: AbsGenOrdering
  vars::Vector{Int}
  function SymbOrdering(S::Symbol, v)
    S in (:lex, :deglex, :degrevlex, :revlex,
          :neglex, :negdeglex, :negdegrevlex, :negrevlex) ||
        throw(ArgumentError("unsupported ordering $S"))
    return new{S}(v)
  end
end

# wdeglex, wdegrevlex, ...
struct WSymbOrdering{S} <: AbsGenOrdering
  vars::Vector{Int}
  weights::Vector{Int}
  function WSymbOrdering(S::Symbol, v, w::Vector{Int})
    S in (:wdeglex, :wdegrevlex, :negwdeglex, :negwdegrevlex) ||
        throw(ArgumentError("unsupported ordering $S"))
    length(v) == length(w) ||
        throw(ArgumentError("number of variables should match the number of weights"))
    all(>(0), v) || throw(ArgumentError("all weights should be positive"))
    return new{S}(v, w)
  end
end

# in general denotes a partial ordering on vars
# if fullrank is true, then the matrix should be invertible
struct MatrixOrdering <: AbsGenOrdering
  vars::Vector{Int}
  matrix::fmpz_mat
  fullrank::Bool
  function MatrixOrdering(v, m::fmpz_mat, fullrank::Bool)
    length(v) == ncols(m) ||
        throw(ArgumentError("number of variables should match the number of columns"))
    if fullrank
      if nrows(m) > ncols(m)
        m = _canonical_matrix(m)
      end
      if nrows(m) < ncols(m)
        throw(ArgumentError("weight matrix is rank deficient"))
      else
        @assert nrows(m) == ncols(m)
        !iszero(det(m)) || throw(ArgumentError("weight matrix is not invertible"))
      end
    end
    return new(v, m, fullrank)
  end
end

function _canonical_matrix(w)
  ww = matrix(ZZ, 0, ncols(w), [])
  for i in 1:nrows(w)
    if is_zero_row(w, i)
      continue
    end
    nw = w[i, :]
    c = content(nw)
    if !isone(c)
      nw = divexact(nw, c)
    end
    for j in 1:nrows(ww)
      h = findfirst(x->ww[j, x] != 0, 1:ncols(w))
      if !iszero(nw[1, h])
        nw = abs(ww[j, h])*nw - sign(ww[j, h])*nw[1, h]*ww[j, :]
      end
    end
    if !iszero(nw)
      c = content(nw)
      if !isone(c)
        nw = divexact(nw, c)
      end
      ww = vcat(ww, nw)
    end
  end
  return ww
end


# convert (y,x,z) => (2,1,3) and check uniqueness
function _unique_var_indices(a::AbstractVector{<:MPolyElem})
  !isempty(a) || error("need at least one variable")
  z = Int[var_index(i) for i in a]
  allunique(z) || error("variables must be unique")
  return z
end

function isweighted(ord::Symbol)
   return ord == :wdeglex || ord == :wdegrevlex ||
          ord == :negwdeglex || ord == :negwdegrevlex
end

"""
The product of `a` and `b` (`vcat` of the the matrices)
"""
mutable struct ProdOrdering <: AbsGenOrdering
  a::AbsGenOrdering
  b::AbsGenOrdering
end

Base.:*(a::AbsGenOrdering, b::AbsGenOrdering) = ProdOrdering(a, b)

@doc Markdown.doc"""
    anti_diagonal(R::Ring, n::Int)

A square matrix with `1` on the anti-diagonal.
"""
function anti_diagonal(R::Ring, n::Int)
  a = zero_matrix(R, n, n)
  for i=1:n
    a[i, n-i+1] = one(R)
  end
  return a
end

function _support_indices(o::ProdOrdering)
  i = _support_indices(o.a)
  j = _support_indices(o.b)
  if i == j
    return i
  else
    return union(i, j)
  end
end

function _matrix(nvars::Int, o::ProdOrdering)
  return vcat(_matrix(nvars, o.a), _matrix(nvars, o.b))
end

"""
Orderings actually applied to polynomial rings (as opposed to variable indices)
"""
mutable struct MonomialOrdering{S}
  R::S
  o::AbsGenOrdering
end

base_ring(a::MonomialOrdering) = a.R

@doc Markdown.doc"""
    support(o::MonomialOrdering)

Return the vector of variables on which `o` is defined.
"""
function support(o::MonomialOrdering)
  return [gen(base_ring(o), i) for i in _support_indices(o.o)]
end

@doc Markdown.doc"""
    *(M::MonomialOrdering, N::MonomialOrdering)

The product ordering `M*N` tries to order by `M` first, and in the case of a
tie uses `N`. Corresponds to a vertical concatenation of the weight matrices.
"""
function Base.:*(M::MonomialOrdering, N::MonomialOrdering)
  M.R == N.R || error("wrong rings")
  return MonomialOrdering(M.R, M.o*N.o)
end

######## non-weighted orderings ########

function _support_indices(o::SymbOrdering)
  return o.vars
end

@doc Markdown.doc"""
    monomial_ordering(v::AbstractVector{<:MPolyElem}, s::Symbol)

Defines an ordering to be applied to the variables in `v`. The symbol `s`
should be one of `:lex`, `:deglex`, `:degrevlex`, `:revlex`, `:neglex`,
`:negdeglex`, `:negdegrevlex`, `:negrevlex`.
"""
function monomial_ordering(v::AbstractVector{<:MPolyElem}, s::Symbol)
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(s, i))
end

function monomial_ordering(R::MPolyRing, s::Symbol)
  return MonomialOrdering(R, SymbOrdering(s, 1:nvars(R)))
end

#### lex ####

@doc Markdown.doc"""
    lex(V::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Given a vector `V` of variables, return the lexicographical ordering on the set of monomials in these variables.

    lex(R::MPolyRing) -> MonomialOrdering

Return the lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = lex([w, x])
lex([w, x])

julia> o2 = lex(gens(R)[3:4])
lex([y, z])

julia> o3 = lex(R)
lex([w, x, y, z])
```
"""
function lex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:lex, i))
end

function lex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:lex, 1:nvars(R)))
end

function _matrix(nvars::Int, o::SymbOrdering{:lex})
  m = zero_matrix(ZZ, length(o.vars), nvars)
  i = 1
  for j in o.vars
    m[i, j] = 1
    i += 1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::SymbOrdering{:lex})
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### deglex ####

@doc Markdown.doc"""
    deglex(V::AbstractVector{<:MPolyElem}) -> MonomialOrdering
    
Given a vector `V` of variables, return the degree lexicographical ordering on the set of monomials in these variables.

    deglex(R::MPolyRing) -> MonomialOrdering

Return the degree lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = deglex([w, x])
deglex([w, x])

julia> o2 = deglex(gens(R)[3:4])
deglex([y, z])

julia> o3 = deglex(R)
deglex([w, x, y, z])
```
"""
function deglex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:deglex, i))
end

function deglex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:deglex, 1:nvars(R)))
end

function _matrix(nvars::Int, o::SymbOrdering{:deglex})
  m = zero_matrix(ZZ, 1 + length(o.vars), nvars)
  i = 2
  for j in o.vars
    m[1, j] = 1
    m[i, j] = 1
    i += 1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::SymbOrdering{:deglex})
  tdk = sum(exponent(f, k, i) for i in o.vars)
  tdl = sum(exponent(f, l, i) for i in o.vars)
  if tdk > tdl
    return 1
  elseif tdk < tdl
    return -1
  end
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### degrevlex ####

@doc Markdown.doc"""
    degrevlex(V::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Given a vector `V` of variables, return the degree reverse lexicographical ordering on the set of monomials in these variables
    degrevlex(R::MPolyRing) -> MonomialOrdering

Return the degree reverse lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = degrevlex([w, x])
degrevlex([w, x])

julia> o2 = degrevlex(gens(R)[3:4])
degrevlex([y, z])

julia> o3 = degrevlex(R)
degrevlex([w, x, y, z])
```
"""
function degrevlex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:degrevlex, i))
end

function degrevlex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:degrevlex, 1:nvars(R)))
end

function _matrix(nvars::Int, o::SymbOrdering{:degrevlex})
  m = zero_matrix(ZZ, 1 + length(o.vars), nvars)
  i = 1 + length(o.vars)
  for j in o.vars
    m[1, j] = 1
    m[i, j] = -1
    i -= 1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::SymbOrdering{:degrevlex})
  tdk = sum(exponent(f, k, i) for i in o.vars)
  tdl = sum(exponent(f, l, i) for i in o.vars)
  if tdk > tdl
    return 1
  elseif tdk < tdl
    return -1
  end
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### revlex ####

@doc Markdown.doc"""
    revlex(V::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Given a vector `V` of variables, return the reverse lexicographical ordering on the set of monomials in these variables.

    revlex(R::MPolyRing) -> MonomialOrdering

Return the reverse lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = revlex([w, x])
revlex([w, x])

julia> o2 = revlex(gens(R)[3:4])
revlex([y, z])

julia> o3 = revlex(R)
revlex([w, x, y, z])
```
"""
function revlex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:revlex, i))
end

function revlex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:revlex, 1:nvars(R)))
end

function _matrix(nvars::Int, o::SymbOrdering{:revlex})
  m = zero_matrix(ZZ, length(o.vars), nvars)
  i = length(o.vars)
  for j in o.vars
    m[i, j] = 1
    i -= 1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::SymbOrdering{:revlex})
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### neglex ####

@doc Markdown.doc"""
    neglex(V::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Given a vector `V` of variables, return the negative lexicographical ordering on the set of monomials in these variables.

    neglex(R::MPolyRing) -> MonomialOrdering

Return the negative lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = neglex([w, x])
neglex([w, x])

julia> o2 = neglex(gens(R)[3:4])
neglex([y, z])

julia> o3 = neglex(R)
neglex([w, x, y, z])
```
"""
function neglex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:neglex, i))
end

function neglex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:neglex, 1:nvars(R)))
end

function _matrix(nvars::Int, o::SymbOrdering{:neglex})
  m = zero_matrix(ZZ, length(o.vars), nvars)
  i = 1
  for j in o.vars
    m[i, j] = -1
    i += 1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::SymbOrdering{:neglex})
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### negrevlex ####

@doc Markdown.doc"""
    negrevlex(V::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Given a vector `V` of variables, return the negative reverse lexicographical ordering on the set of monomials in these variables.

    negrevlex(R::MPolyRing) -> MonomialOrdering

Return the negative reverse lexicographical ordering  on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = negrevlex([w, x])
negrevlex([w, x])

julia> o2 = negrevlex(gens(R)[3:4])
negrevlex([y, z])

julia> o3 = negrevlex(R)
negrevlex([w, x, y, z])
```
"""
function negrevlex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:negrevlex, i))
end

function negrevlex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:negrevlex, 1:nvars(R)))
end

function _matrix(nvars::Int, o::SymbOrdering{:negrevlex})
  m = zero_matrix(ZZ, length(o.vars), nvars)
  i = length(o.vars)
  for j in o.vars
    m[i, j] = -1
    i -= 1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::SymbOrdering{:negrevlex})
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### negdegrevlex ####

@doc Markdown.doc"""
    negdegrevlex(V::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Given a vector `V` of variables, return the negative degree reverse lexicographical ordering on the set of monomials in these variables.

    negdegrevlex(R::MPolyRing) -> MonomialOrdering

Return the negative degree reverse lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = negdegrevlex([w, x])
negdegrevlex([w, x])

julia> o2 = negdegrevlex(gens(R)[3:4])
negdegrevlex([y, z])

julia> o3 = negdegrevlex(R)
negdegrevlex([w, x, y, z])
```
"""
function negdegrevlex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:negdegrevlex, i))
end

function negdegrevlex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:negdegrevlex, 1:nvars(R)))
end

function _matrix(nvars::Int, o::SymbOrdering{:negdegrevlex})
  m = zero_matrix(ZZ, 1 + length(o.vars), nvars)
  i = 1 + length(o.vars)
  for j in o.vars
    m[1, j] = -1
    m[i, j] = -1
    i -= 1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::SymbOrdering{:negdegrevlex})
  tdk = sum(exponent(f, k, i) for i in o.vars)
  tdl = sum(exponent(f, l, i) for i in o.vars)
  if tdk > tdl
    return -1
  elseif tdk < tdl
    return 1
  end
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### negdeglex ####

@doc Markdown.doc"""
    negdeglex(V::AbstractVector{<:MPolyElem}) -> MonomialOrdering    

Given a vector `V` of variables, return the negative degree lexicographical ordering on the set of monomials in these variables.

    negdeglex(R::MPolyRing) -> MonomialOrdering

Return the negative degree lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = negdeglex([w, x])
negdeglex([w, x])

julia> o2 = negdeglex(gens(R)[3:4])
negdeglex([y, z])

julia> o3 = negdeglex(R)
negdeglex([w, x, y, z])
```
"""
function negdeglex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:negdeglex, i))
end

function negdeglex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:negdeglex, 1:nvars(R)))
end

function _matrix(nvars::Int, o::SymbOrdering{:negdeglex})
  m = zero_matrix(ZZ, 1 + length(o.vars), nvars)
  i = 2
  for j in o.vars
    m[1, j] = -1
    m[i, j] = 1
    i += 1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::SymbOrdering{:negdeglex})
  tdk = sum(exponent(f, k, i) for i in o.vars)
  tdl = sum(exponent(f, l, i) for i in o.vars)
  if tdk > tdl
    return -1
  elseif tdk < tdl
    return 1
  end
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

######## weighted orderings ########

function _support_indices(o::WSymbOrdering)
  return o.vars
end

@doc Markdown.doc"""
    monomial_ordering(v::AbstractVector{<:MPolyElem}, s::Symbol, w::Vector{Int})
    monomial_ordering(R::MPolyRing, s::Symbol, w::Vector{Int}) -> MonomialOrdering

Defines a weighted ordering to be applied to the variables in `v`. The weight
vector `w` should be the same length as `v`, and the symbol `s` should be one
of `:wdeglex`, `:wdegrevlex`, `:negwdeglex`, `:negwdegrevlex`.
"""
function monomial_ordering(v::AbstractVector{<:MPolyElem}, s::Symbol, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(s, i, w))
end

function monomial_ordering(R::MPolyRing, s::Symbol, w::Vector{Int})
  return MonomialOrdering(R, WSymbOrdering(s, 1:nvars(R), w))
end

function _cmp_weighted_degree(f::MPolyElem, k::Int, l::Int, vars::Vector{Int}, w::Vector{Int})
  n = length(vars)
  @assert n == length(w)
  t = 0
  for j in 1:n
    v = vars[j]
    t += w[j]*(exponent(f, k, v) - exponent(f, l, v))
  end
  return t > 0 ? 1 : t < 0 ? -1 : 0;
end

#### wdeglex, Wp ####

@doc Markdown.doc"""
    wdeglex(V::AbstractVector{<:MPolyElem}, W::Vector{Int}) -> MonomialOrdering
    
Given a vector `V` of variables and a vector `W` of positive integers, return the corresponding weighted 
lexicographical ordering on the set of monomials in the given variables.

    wdeglex(R, W::Vector{Int}) -> MonomialOrdering

If `W` is a vector of positive integers, return the corresponding weighted 
lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = wdeglex([w, x], [1, 2])
wdeglex([w, x], [1, 2])

julia> o2 = wdeglex(gens(R)[3:4], [3, 4])
wdeglex([y, z], [3, 4])

julia> o3 = wdeglex(R, [1, 2, 3, 4])
wdeglex([w, x, y, z], [1, 2, 3, 4])
```
"""
function wdeglex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:wdeglex, i, w))
end

function wdeglex(R::MPolyRing, w::Vector{Int})
  return MonomialOrdering(R, WSymbOrdering(:wdeglex, 1:nvars(R), w))
end

function _matrix(nvars::Int, o::WSymbOrdering{:wdeglex})
  n = length(o.vars)
  m = zero_matrix(ZZ, 1 + n, nvars)
  for i in 1:n
    j = o.vars[i]
    m[1, j] = o.weights[i]
    m[1 + i, j] = 1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::WSymbOrdering{:wdeglex})
  c = _cmp_weighted_degree(f, k, l, o.vars, o.weights)
  c == 0 || return c
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### wdegrevlex, wp ####

@doc Markdown.doc"""
    wdegrevlex(V::AbstractVector{<:MPolyElem}, W::Vector{Int}) -> MonomialOrdering

Given a vector `V` of variables and a vector `W` of positive integers, return the corresponding weighted reverse 
lexicographical ordering on the set of monomials in the given variables.

    wdegrevlex(R::MPolyRing, W::Vector{Int}) -> MonomialOrdering

If `W` is a vector of positive integers, return the corresponding weighted reverse 
lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = wdegrevlex([w, x], [1, 2])
wdegrevlex([w, x], [1, 2])

julia> o2 = wdegrevlex(gens(R)[3:4], [3, 4])
wdegrevlex([y, z], [3, 4])

julia> o3 = wdegrevlex(R, [1, 2, 3, 4])
wdegrevlex([w, x, y, z], [1, 2, 3, 4])
```
"""
function wdegrevlex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:wdegrevlex, i, w))
end

function wdegrevlex(R::MPolyRing, w::Vector{Int})
  return MonomialOrdering(R, WSymbOrdering(:wdegrevlex, 1:nvars(R), w))
end

function _matrix(nvars::Int, o::WSymbOrdering{:wdegrevlex})
  n = length(o.vars)
  m = zero_matrix(ZZ, 1 + n, nvars)
  for i in 1:n
    j = o.vars[i]
    m[1, j] = o.weights[i]
    m[2 + n - i, j] = -1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::WSymbOrdering{:wdegrevlex})
  c = _cmp_weighted_degree(f, k, l, o.vars, o.weights)
  c == 0 || return c
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### negwdeglex, Ws ####

@doc Markdown.doc"""
    negwdeglex(V::AbstractVector{<:MPolyElem}, W::Vector{Int}) -> MonomialOrdering
    
Given a vector `V` of variables and a vector `W` of positive integers, return the corresponding
negative weighted lexicographical ordering on the set of monomials in the given variables.

    negwdeglex(R::MPolyRing, W::Vector{Int}) -> MonomialOrdering

If `W` is a vector of positive integers, return the corresponding negative weighted 
lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = negwdeglex([w, x], [1, 2])
negwdeglex([w, x], [1, 2])

julia> o2 = negwdeglex(gens(R)[3:4], [3, 4])
negwdeglex([y, z], [3, 4])

julia> o3 = negwdeglex(R, [1, 2, 3, 4])
negwdeglex([w, x, y, z], [1, 2, 3, 4])
```
"""
function negwdeglex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:negwdeglex, i, w))
end

function negwdeglex(R::MPolyRing, w::Vector{Int})
  return MonomialOrdering(R, WSymbOrdering(:negwdeglex, 1:nvars(R), w))
end

function _matrix(nvars::Int, o::WSymbOrdering{:negwdeglex})
  n = length(o.vars)
  m = zero_matrix(ZZ, 1 + n, nvars)
  for i in 1:n
    j = o.vars[i]
    m[1, j] = -o.weights[i]
    m[1 + i, j] = 1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::WSymbOrdering{:negwdeglex})
  c = _cmp_weighted_degree(f, l, k, o.vars, o.weights)
  c == 0 || return c
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### negwdegrevlex, ws ####

@doc Markdown.doc"""
    negwdegrevlex(V::AbstractVector{<:MPolyElem}, W::Vector{Int}) -> MonomialOrdering
    
Given a vector `V` of variables and a vector `W` of positive integers, return the corresponding negative
weighted reverse lexicographical ordering on the set of monomials in the given variables.

    negwdegrevlex(R::MPolyRing, W::Vector{Int}) -> MonomialOrdering

If `W` is a vector of positive integers, return the corresponding negative weighted
reverse lexicographical ordering on the set of monomials in the variables of `R`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
(Multivariate Polynomial Ring in w, x, y, z over Rational Field, fmpq_mpoly[w, x, y, z])

julia> o1 = negwdegrevlex([w, x], [1, 2])
negwdegrevlex([w, x], [1, 2])

julia> o2 = negwdegrevlex(gens(R)[3:4], [3, 4])
negwdegrevlex([y, z], [3, 4])

julia> o3 = negwdegrevlex(R, [1, 2, 3, 4])
negwdegrevlex([w, x, y, z], [1, 2, 3, 4])
```
"""
function negwdegrevlex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:negwdegrevlex, i, w))
end

function negwdegrevlex(R::MPolyRing, w::Vector{Int})
  return MonomialOrdering(R, WSymbOrdering(:negwdegrevlex, 1:nvars(R), w))
end

function _matrix(nvars::Int, o::WSymbOrdering{:negwdegrevlex})
  n = length(o.vars)
  m = zero_matrix(ZZ, 1 + n, nvars)
  for i in 1:n
    j = o.vars[i]
    m[1, j] = -o.weights[i]
    m[2 + n - i, j] = -1
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::WSymbOrdering{:negwdegrevlex})
  c = _cmp_weighted_degree(f, l, k, o.vars, o.weights)
  c == 0 || return c
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### matrix, M ####

function _support_indices(o::MatrixOrdering)
  return o.vars
end

@doc Markdown.doc"""
    matrix_ordering(V::AbstractVector{<:MPolyElem}, M::Union{Matrix{T}, MatElem{T}}; check = true) where T -> MonomialOrdering

Given a vector `V` of variables and an integer matrix `M` such that `length(V) = ncols(M) = rank(M)`, 
return the corresponding matrix ordering on the set of monomials in the given variables. 

    matrix_ordering(R::MPolyRing, M::Union{Matrix{T}, MatElem{T}}; check = true) -> MonomialOrdering

Given an integer matrix `M` such that `nvars(R) = ncols(M) = rank(M)`, 
return the matrix ordering on the set of variables of `R` which is defined by `M`.

!!! note
    The matrix `M` need not be square.

!!! note
    If `check = false` is supplied, the rank check is omitted, and the resulting
    ordering may only be partial.

# Examples
```jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"];

julia> M =[1 1 1; 0 0 -1; 0 -1 0]
3×3 Matrix{Int64}:
 1   1   1
 0   0  -1
 0  -1   0

julia> o = matrix_ordering(R, M)
matrix_ordering([x, y, z], [1 1 1; 0 0 -1; 0 -1 0])

julia> o == degrevlex(R)
true

julia> canonical_matrix(matrix_ordering([x, y, z], [1 2 3; 2 4 6]; check = false))
[1   2   3]
```
"""
function matrix_ordering(v::AbstractVector{<:MPolyElem}, M::Union{Matrix{T}, MatElem{T}}; check = true) where T
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), MatrixOrdering(i, fmpz_mat(M), check))
end

function matrix_ordering(R::MPolyRing, M::Union{Matrix{T}, MatElem{T}}; check = true) where T
  return MonomialOrdering(R, MatrixOrdering(1:nvars(R), fmpz_mat(M), check))
end

@doc Markdown.doc"""
    weight_ordering(w::Vector{Int}, ord::MonomialOrdering)

Returns the ordering on `support(ord)` obtained by first comparing the
`w`-weighted degrees and then using `ord` in the case of tie.
"""
function weight_ordering(w::Vector{Int}, o::MonomialOrdering)
  i = _support_indices(o.o)
  m = fmpz_mat(1, length(w), w)
  return MonomialOrdering(base_ring(o), MatrixOrdering(i, m, false))*o
end

function _matrix(nvars::Int, o::MatrixOrdering)
  r = nrows(o.matrix)
  c = length(o.vars)
  @assert c == ncols(o.matrix)
  m = zero_matrix(ZZ, r, nvars)
  for i in 1:r, j in 1:c
    m[i, o.vars[j]] = o.matrix[i, j]
  end
  return m
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::MatrixOrdering)
  M = o.matrix
  r = nrows(M)
  c = ncols(M)
  for i in 1:r
    v = o.vars[1]
    t = M[i,1]*(exponent(f, k, v) - exponent(f, l, v))
    for j in 2:c
      v = o.vars[j]
      t += M[i, j]*(exponent(f, k, v) - exponent(f, l, v))
    end
    if t > 0
      return 1
    elseif t < 0
      return -1
    end
  end
  return 0 
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::ProdOrdering)
  cmp = _cmp_monomials(f, k, l, o.a)
  cmp == 0 || return cmp
  return _cmp_monomials(f, k, l, o.b)
end

###################################

@doc Markdown.doc"""
    matrix(M::MonomialOrdering)

Return a corresponding weight matrix for the given ordering.
"""
function matrix(M::MonomialOrdering)
  return _matrix(nvars(base_ring(M)), M.o)
end

@doc Markdown.doc"""
    simplify(M::MonomialOrdering) -> MonomialOrdering

Returns a matrix ordering with a unique weight matrix.
"""
function Hecke.simplify(M::MonomialOrdering)
  w = canonical_matrix(M)
  return MonomialOrdering(M.R, MatrixOrdering(collect(1:ncols(w)), w, true))
end

function canonical_matrix(nvars::Int, M::AbsOrdering)
  return _canonical_matrix(_matrix(nvars, M))
end

@doc Markdown.doc"""
    canonical_matrix(M::MonomialOrdering)

Return the corresponding canonical weight matrix for the given ordering.
"""
function canonical_matrix(M::MonomialOrdering)
  return canonical_matrix(nvars(base_ring(M)), M.o)
end

import Base.==
function ==(M::MonomialOrdering, N::MonomialOrdering)
  return canonical_matrix(M) == canonical_matrix(N)
end

function Base.hash(M::MonomialOrdering, u::UInt)
  return hash(canonical_matrix(M), u)
end

function _cmp_var(M, j::Int)
  for i in 1:nrows(M)
    e = M[i,j]
    e > 0 && return +1
    e < 0 && return -1
  end
  return 0
end

@doc Markdown.doc"""
    is_global(ord::MonomialOrdering)

Return `true` if `ord` is global, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> o = matrix_ordering([x, y], [1 1; 0 -1])
matrix_ordering([x, y], [1 1; 0 -1])

julia> is_global(o)
true
```
"""
function is_global(ord::MonomialOrdering)
  M = matrix(ord)
  for i in _support_indices(ord.o)
    if _cmp_var(M, i) <= 0
      return false
    end
  end
  return true
end

@doc Markdown.doc"""
    is_local(ord::MonomialOrdering)

Return `true` if `v<1` for each `v` in `support(ord)`, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> o = matrix_ordering([x, y], [-1 -1; 0 -1])
matrix_ordering([x, y], [-1 -1; 0 -1])

julia> is_local(o)
true
```
"""
function is_local(ord::MonomialOrdering)
  M = matrix(ord)
  for i in _support_indices(ord.o)
    if _cmp_var(M, i) >= 0
      return false
    end
  end
  return true
end

@doc Markdown.doc"""
    is_mixed(ord::MonomialOrdering)

Return `!is_global(ord) && !is_local(ord)`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> o = matrix_ordering([x, y], [1 -1; 0 -1])
matrix_ordering([x, y], [1 -1; 0 -1])

julia> is_mixed(o)
true
```
"""
function is_mixed(ord::MonomialOrdering)
  M = matrix(ord)
  s = _support_indices(ord.o)
  for i in s
    if _cmp_var(M, i) <= 0
      @goto is_not_global
    end
  end
  return false
@label is_not_global
  for i in s
    if _cmp_var(M, i) >= 0
      return true
    end
  end
  return false
end


# ex: nvars = 7, sigmaC = {    2, 3,       6   }
#             =>  sigma = { 1,       4, 5,    7}
#                varmap = {-1,+1,+2,-2,-3,+3,-4}
function _elimination_data(n::Int, sigmaC::Vector)
  varmap = zeros(Int, n)
  for si in sigmaC
    varmap[si] = 1
  end
  sigma = Int[]
  sigmaC = Int[]
  for i in 1:n
    if varmap[i] == 0
      push!(sigma, i)
      varmap[i] = -length(sigma)
    else
      push!(sigmaC, i)
      varmap[i] = +length(sigmaC)
    end
  end
  return varmap, sigma, sigmaC
end

@doc Markdown.doc"""
    is_elimination_ordering(M::MonomialOrdering, v::Vector{Int})
    is_elimination_ordering(M::MonomialOrdering, v::Vector{<:MPolyElem})

Return `true` if the given ordering eliminates the variables in `v`.
"""
function is_elimination_ordering(o::MonomialOrdering, sigmaC::Vector{Int})
# This is probably correct: o is elimination ordering for sigma^C <=>
#   For each j in sigma^C, the first non-zero weight M[i,j] in column j is >= 0,
#   and the weights in columns corresponding to sigma in rows <= i are zero.
  M = canonical_matrix(o)
  m = nrows(M)
  n = ncols(M)
  varmap, sigma, sigmaC = _elimination_data(n, sigmaC)
  maxrow = 1
  for j in sigmaC
    found = false
    for i in 1:m
      Mij = M[i,j]
      if Mij > 0
        maxrow = max(maxrow, i)
        found = true
        break # continue j loop
      elseif Mij < 0
        return false
      end
    end
    found || return false
  end
  for i in 1:maxrow
    for k in sigma
      iszero(M[i,k]) || return false
    end
  end
  return true
end

function is_elimination_ordering(o::MonomialOrdering, sigmaC::Vector{<:MPolyElem})
  return is_elimination_ordering(o, Int[var_index(v) for v in sigmaC])
end

###################################################

# Module orderings (not module Orderings)

mutable struct ModOrdering{T} <: AbsModOrdering
   gens::T
   ord::Symbol
   function ModOrdering(u::T, s::Symbol) where {T <: AbstractVector{Int}}
     r = new{T}()
     r.gens = u
     r.ord = s
     return r
   end
end

mutable struct ModuleOrdering{S}
   M::S
   o::AbsOrdering # must allow gen*mon or mon*gen product ordering
end

base_ring(a::ModuleOrdering) = a.M

mutable struct ModProdOrdering <: AbsModOrdering
   a::AbsOrdering
   b::AbsOrdering
end

Base.:*(a::AbsGenOrdering, b::AbsModOrdering) = ModProdOrdering(a, b)

Base.:*(a::AbsModOrdering, b::AbsGenOrdering) = ModProdOrdering(a, b)

function module_ordering(a::AbstractVector{Int}, s::Symbol)
   i = minimum(a)
   I = maximum(a)
   if a == i:I #test if variables are consecutive or not.
     return ModOrdering(i:I, s)
   end
   return ModOrdering(collect(a), s)
 end

function ordering(a::AbstractVector{<:AbstractAlgebra.ModuleElem}, s...)
   R = parent(first(a))
   g = gens(R)
   aa = [findfirst(x -> x == y, g) for y = a]
   if nothing in aa
     error("only generators allowed")
   end
   return module_ordering(aa, s...)
 end

 function lex(v::AbstractVector{<:AbstractAlgebra.ModuleElem})
   return ModuleOrdering(parent(first(v)), ordering(v, :lex))
end

function revlex(v::AbstractVector{<:AbstractAlgebra.ModuleElem})
   return ModuleOrdering(parent(first(v)), ordering(v, :revlex))
end

function Base.:*(M::ModuleOrdering, N::MonomialOrdering)
   base_ring(M.M) == N.R || error("wrong rings")
   return ModuleOrdering(M.M, M.o*N.o)
end

function Base.:*(M::MonomialOrdering, N::ModuleOrdering)
   base_ring(N.M) == M.R || error("wrong rings")
   return ModuleOrdering(N.M, M.o*N.o)
end

@doc Markdown.doc"""
    induced_ring_ordering(M::ModuleOrdering)

Return the induced ring ordering.
"""
function induced_ring_ordering(M::ModuleOrdering)
  R = base_ring(M.M)
  if M.o.a isa AbsGenOrdering
    return MonomialOrdering(R, M.o.a)
  else
    return MonomialOrdering(R, M.o.b)
  end
end

@doc Markdown.doc"""
    is_global(M::ModuleOrdering)

Return `true` if the given ordering is global, i.e. if the
induced ring ordering is global.
"""
function is_global(M::ModuleOrdering)
  return is_global(induced_ring_ordering(M))
end


@doc Markdown.doc"""
    permutation_of_terms(f::MPolyElem, ord::MonomialOrdering)

Return the permutation that puts the terms of `f` in the order `ord`.
"""
function permutation_of_terms(f::MPolyElem, ord::MonomialOrdering)
  p = collect(1:length(f))
  sort!(p, lt = (k, l) -> (Orderings._cmp_monomials(f, k, l, ord.o) < 0), rev = true)
  return p
end

############ printing ############

function _expressify(o::SymbOrdering{S}, sym)  where S
  return Expr(:call, S, Expr(:vect, (sym[i] for i in o.vars)...))
end

function _expressify(o::WSymbOrdering{S}, sym)  where S
  return Expr(:call, S, Expr(:vect, (sym[i] for i in o.vars)...),
                        Expr(:vect, o.weights...))
end

function _expressify(o::MatrixOrdering, sym)
  return Expr(:call, :matrix_ordering, Expr(:vect, (sym[i] for i in o.vars)...),
                                          AbstractAlgebra.expressify(o.matrix))
end

function _expressify(o::Union{ProdOrdering, ModProdOrdering}, sym)
  return Expr(:call, :*, _expressify(o.a, sym), _expressify(o.b, sym))
end

function AbstractAlgebra.expressify(M::MonomialOrdering; context = nothing)
  return _expressify(M.o, symbols(base_ring(M)))
end

@enable_all_show_via_expressify MonomialOrdering

function _expressify(o::ModOrdering, sym)
  return Expr(:call, o.ord, Expr(:vect, (Expr(:call, :gen, i) for i in o.gens)...))
end

function AbstractAlgebra.expressify(M::ModuleOrdering; context = nothing)
  return _expressify(M.o, symbols(base_ring(base_ring(M))))
end

@enable_all_show_via_expressify ModuleOrdering

############ opposite ordering ############

# TODO return singular-friendly orderings if that is how they come in
function _opposite_ordering(nvars::Int, o::SymbOrdering{T}) where T
  return SymbOrdering(T, nvars+1 .- o.vars)
end

# TODO ditto
function _opposite_ordering(nvars::Int, o::WSymbOrdering{T}) where T
  return WSymbOrdering(T, nvars+1 .- o.vars, o.weights)
end

function _opposite_ordering(nvars::Int, o::MatrixOrdering)
  M = o.matrix
  M = reduce(hcat, [M[:,i] for i in ncols(M):-1:1])
  return MatrixOrdering(reverse(nvars+1 .- o.vars), M, false)
end

function _opposite_ordering(n::Int, o::ProdOrdering)
  return ProdOrdering(_opposite_ordering(n, o.a), _opposite_ordering(n, o.b))
end

@doc Markdown.doc"""
    opposite_ordering(R::MPolyRing, o::MonomialOrdering)

Return an ordering on `R` whose weight matrix is the column-wise reverse of the
weight matrix of `o`.
"""
function opposite_ordering(R::MPolyRing, o::MonomialOrdering)
  @assert nvars(R) == nvars(base_ring(o))
  return MonomialOrdering(R, _opposite_ordering(nvars(R), o.o))
end

############ Singular conversions ############

#### Oscar -> Singular ####

mutable struct order_conversion_ctx
  last_var::Int
  has_c_or_C::Bool
  def::Singular.sordering
end

function _is_consecutive_from(a::Vector{Int}, b::Int)
  n = length(a)
  n > 0 || return false
  for i in 1:n
    a[i] == b + i || return false
  end
  return true
end

function _try_singular_easy(Q::order_conversion_ctx, o::SymbOrdering{S}) where S
  _is_consecutive_from(o.vars, Q.last_var) || return (false, Q.def)
  n = length(o.vars)
  Q.last_var += n
  return S == :lex          ? (true, Singular.ordering_lp(n)) :
         S == :revlex       ? (true, Singular.ordering_rp(n)) :
         S == :deglex       ? (true, Singular.ordering_Dp(n)) :
         S == :degrevlex    ? (true, Singular.ordering_dp(n)) :
         S == :neglex       ? (true, Singular.ordering_ls(n)) :
         S == :negrevlex    ? (true, Singular.ordering_rs(n)) :
         S == :negdeglex    ? (true, Singular.ordering_Ds(n)) :
         S == :negdegrevlex ? (true, Singular.ordering_ds(n)) :
                              (false, Q.def)
end

function _try_singular_easy(Q::order_conversion_ctx, o::WSymbOrdering{S}) where S
  _is_consecutive_from(o.vars, Q.last_var) || return (false, Q.def)
  n = length(o.vars)
  Q.last_var += n
  return S == :wdeglex       ? (true, Singular.ordering_Wp(o.weights)) :
         S == :wdegrevlex    ? (true, Singular.ordering_wp(o.weights)) :
         S == :negwdeglex    ? (true, Singular.ordering_Ws(o.weights)) :
         S == :negwdegrevlex ? (true, Singular.ordering_ws(o.weights)) :
                               (false, Q.def)
end

function _try_singular_easy(Q::order_conversion_ctx, o::MatrixOrdering)
  _is_consecutive_from(o.vars, Q.last_var) || return (false, Q.def)
  n = length(o.vars)
  if nrows(o.matrix) == n && (o.fullrank || !iszero(det(o.matrix)))
    Q.last_var += n
    return (true, Singular.ordering_M(o.matrix; check=false))
  elseif nrows(o.matrix) == 1
    return (true, Singular.ordering_a([Int(o.matrix[1,i]) for i in 1:n]))
  end
  return (false, Q.def)
end

function _try_singular_easy(Q::order_conversion_ctx, o::Orderings.ModOrdering)
  Q.has_c_or_C && return (false, Q.def)
  Q.has_c_or_C = true
  o.gens == 1:length(o.gens) || return (false, Q.def)
  return o.ord == :lex    ? (true, Singular.ordering_C(length(o.gens))) :
         o.ord == :revlex ? (true, Singular.ordering_c(length(o.gens))) :
                            (false, Q.def)
end

function _try_singular_easy(Q::order_conversion_ctx, o::Union{ProdOrdering, ModProdOrdering})
  (ok, a) = _try_singular_easy(Q, o.a)
  ok || return (false, Q.def)
  (ok, b) = _try_singular_easy(Q, o.b)
  ok || return (false, Q.def)
  return (true, a*b)
end

function singular(ord::MonomialOrdering)
  Q = order_conversion_ctx(0, false, Singular.ordering_lp())
  (ok, sord) = _try_singular_easy(Q, ord.o)
  if ok && Q.last_var == nvars(base_ring(ord))
    return sord
  else
    # Singular.jl checks det != 0 (and is square)
    return Singular.ordering_M(canonical_matrix(ord))
  end
end

function singular(ord::ModuleOrdering)
  Q = order_conversion_ctx(0, false, Singular.ordering_lp())
  (ok, sord) = _try_singular_easy(Q, ord.o)
  n = nvars(base_ring(base_ring(ord)))
  if ok && Q.last_var == n
    @assert Q.has_c_or_C
    return sord
  else
    error("failed to convert module ordering")
  end
end

#### Singular -> Oscar, uses unofficial parts of Singular.jl ####

function _convert_sblock(nvars::Int, o::Singular.sorder_block, lastvar::Int)
  newlastvar = lastvar+o.size
  i = collect(lastvar+1:newlastvar)
  if o.order == Singular.ringorder_lp
    return SymbOrdering(:lex, i), newlastvar
  elseif o.order == Singular.ringorder_rp
    return SymbOrdering(:revlex, i), newlastvar
  elseif o.order == Singular.ringorder_Dp
    return SymbOrdering(:deglex, i), newlastvar
  elseif o.order == Singular.ringorder_dp
    return SymbOrdering(:degrevlex, i), newlastvar
  elseif o.order == Singular.ringorder_ls
    return SymbOrdering(:neglex, i), newlastvar
  elseif o.order == Singular.ringorder_rs
    return SymbOrdering(:negrevlex, i), newlastvar
  elseif o.order == Singular.ringorder_Ds
    return SymbOrdering(:negdeglex, i), newlastvar
  elseif o.order == Singular.ringorder_ds
    return SymbOrdering(:negdegrevlex, i), newlastvar

  elseif o.order == Singular.ringorder_Wp
    return WSymbOrdering(:wdeglex, i, o.weights), newlastvar
  elseif o.order == Singular.ringorder_wp
    return WSymbOrdering(:wdegrevlex, i, o.weights), newlastvar
  elseif o.order == Singular.ringorder_Ws
    return WSymbOrdering(:negwdeglex, i, o.weights), newlastvar
  elseif o.order == Singular.ringorder_ws
    return WSymbOrdering(:negwdegrevlex, i, o.weights), newlastvar

  elseif o.order == Singular.ringorder_a
    # just adds a row to the matrix without increasing `lastvar`
    newlastvar = lastvar+length(o.weights)
    newlastvar <= nvars || error("too many weights in Singular.ordering_a")
    i = collect(lastvar+1:newlastvar)
    m = fmpz_mat(1, length(o.weights), o.weights)
    return MatrixOrdering(i, m, false), lastvar
  elseif o.order == Singular.ringorder_M
    m = fmpz_mat(o.size, o.size, o.weights)
    return MatrixOrdering(i, m, true), newlastvar

  elseif o.order == Singular.ringorder_C
    return nothing, lastvar
  elseif o.order == Singular.ringorder_c
    return nothing, lastvar

  else
    error("cannot convert singular block $(o.order)")
  end
end

@doc Markdown.doc"""
    monomial_ordering(R::MPolyRing, ord::Singular.sordering)

Return an ordering on `R` equivalent to the Singular.jl ordering `ord`.
"""
function monomial_ordering(R::MPolyRing, ord::Singular.sordering)
  n = nvars(R)
  z, lastvar = nothing, 0
  for i in 1:length(ord.data)
    x, lastvar = _convert_sblock(n, ord.data[i], lastvar)
    isnothing(x) && continue
    z = isnothing(z) ? x : ProdOrdering(z, x)
  end
  lastvar == n || error("number of variables in ordering does not match")
  return MonomialOrdering(R, z::AbsGenOrdering)
end

end  # module Orderings
