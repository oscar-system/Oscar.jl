module Orderings

using Oscar, Markdown
import Oscar: Ring, MPolyRing, MPolyElem, weights, IntegerUnion, base_ring
export anti_diagonal, lex, degrevlex, deglex, revlex, negdeglex,
       neglex, negrevlex, negdegrevlex, wdeglex, wdegrevlex,
       negwdeglex, negwdegrevlex, matrix_ordering, monomial_ordering,
       weight_matrix, isweighted, is_global, is_local, is_mixed,
       permutation_of_terms, weighted_ordering, canonical_weight_matrix,
       MonomialOrdering, ModuleOrdering, singular

abstract type AbsOrdering end

abstract type AbsGenOrdering <: AbsOrdering end

abstract type AbsModOrdering <: AbsOrdering end

# lex, deglex, ...
struct SymbOrdering{S} <: AbsGenOrdering
  vars::Vector{Int}
  function SymbOrdering(S::Symbol, v::Vector{Int})
    return new{S}(v)
  end
end

# wdeglex, wdegrevlex, ...
struct WSymbOrdering{S} <: AbsGenOrdering
  vars::Vector{Int}
  weights::Vector{Int}
  function WSymbOrdering(S::Symbol, v::Vector{Int}, w::Vector{Int})
    @assert length(v) == length(w)
    return new{S}(v, w)
  end
end

struct MatrixOrdering <: AbsGenOrdering
  vars::Vector{Int}
  matrix::fmpz_mat
  function MatrixOrdering(v::Vector{Int}, m::fmpz_mat)
    @assert length(v) == ncols(m)
    return new(v, m)
  end
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

function flat(a::SymbOrdering)
  return [a]
end

function flat(a::WSymbOrdering)
  return [a]
end

function flat(a::MatrixOrdering)
  return [a]
end

function flat(a::ProdOrdering)
   return vcat(flat(a.a), flat(a.b))
end  

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



function _weight_matrix(nvars::Int, o::ProdOrdering)
  return vcat(_weight_matrix(nvars, o.a), _weight_matrix(nvars, o.b))
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
    :*(M::MonomialOrdering, N::MonomialOrdering)

The product ordering `M*N` tries to order by `M` first, and in the case of a
tie uses `N`. Corresponds to a vertical concatenation of the weight matrices.
"""
function Base.:*(M::MonomialOrdering, N::MonomialOrdering)
  M.R == N.R || error("wrong rings")
  return MonomialOrdering(M.R, M.o*N.o)
end

######## non-weighted orderings ########

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

#### lex ####

@doc Markdown.doc"""
    lex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `lex` (lexicographic) ordering on the variables given.
"""
function lex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:lex, i))
end

function _weight_matrix(nvars::Int, o::SymbOrdering{:lex})
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
    deglex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `deglex` ordering on the variables given.
"""
function deglex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:deglex, i))
end

function _weight_matrix(nvars::Int, o::SymbOrdering{:deglex})
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
    degrevlex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `degrevlex` ordering on the variables given.
"""
function degrevlex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:degrevlex, i))
end

function _weight_matrix(nvars::Int, o::SymbOrdering{:degrevlex})
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
    revlex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `revlex` ordering on the variables given.
"""
function revlex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:revlex, i))
end

function _weight_matrix(nvars::Int, o::SymbOrdering{:revlex})
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
    neglex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `neglex` ordering on the variables given.
"""
function neglex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:neglex, i))
end

function _weight_matrix(nvars::Int, o::SymbOrdering{:neglex})
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
    negrevlex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `negrevlex` ordering on the variables given.
"""
function negrevlex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:negrevlex, i))
end

function _weight_matrix(nvars::Int, o::SymbOrdering{:negrevlex})
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
    negdegrevlex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `negdegrevlex` ordering on the variables given.
"""
function negdegrevlex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:negdegrevlex, i))
end

function _weight_matrix(nvars::Int, o::SymbOrdering{:negdegrevlex})
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
    negdeglex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `negdeglex` ordering on the variables given.
"""
function negdeglex(v::AbstractVector{<:MPolyElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:negdeglex, i))
end

function _weight_matrix(nvars::Int, o::SymbOrdering{:negdeglex})
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

@doc Markdown.doc"""
    monomial_ordering(v::AbstractVector{<:MPolyElem}, s::Symbol, w::Vector{Int})

Defines a weighted ordering to be applied to the variables in `v`. The weight
vector `w` should be the same length as `v`, and the symbol `s` should be one
of `:wdeglex`, `:wdegrevlex`, `:negwdeglex`, `:negwdegrevlex`.
"""
function monomial_ordering(v::AbstractVector{<:MPolyElem}, s::Symbol, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(s, i, w))
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
    wdeglex(v::AbstractVector{<:MPolyElem}, w::Vector{Int}) -> MonomialOrdering

Defines the `wdeglex` ordering on the variables given with the weights `w`.
"""
function wdeglex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:wdeglex, i, w))
end

function _weight_matrix(nvars::Int, o::WSymbOrdering{:wdeglex})
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
    wdegrevlex(v::AbstractVector{<:MPolyElem}, w::Vector{Int}) -> MonomialOrdering

Defines the `wdegrevlex` ordering on the variables given with the weights `w`.
"""
function wdegrevlex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:wdegrevlex, i, w))
end

function _weight_matrix(nvars::Int, o::WSymbOrdering{:wdegrevlex})
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
    negwdeglex(v::AbstractVector{<:MPolyElem}, w::Vector{Int}) -> MonomialOrdering

Defines the `negwdeglex` ordering on the variables given with the weights `w`.
"""
function negwdeglex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:negwdeglex, i, w))
end

function _weight_matrix(nvars::Int, o::WSymbOrdering{:negwdeglex})
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
    negwdegrevlex(v::AbstractVector{<:MPolyElem}, w::Vector{Int}) -> MonomialOrdering

Defines the `negwdegrevlex` ordering on the variables given with the weights `w`.
"""
function negwdegrevlex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:negwdegrevlex, i, w))
end

function _weight_matrix(nvars::Int, o::WSymbOrdering{:negwdegrevlex})
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

@doc Markdown.doc"""
    matrix_ordering(v::AbstractVector{<:MPolyElem}, M::Union{Matrix{T}, MatElem{T}})

Defines the matrix ordering on the variables given with the matrix `M`. The
matrix need not be square nor have full row rank, thus the resulting ordering
may be only a partial ordering on the given variables.
"""
function matrix_ordering(v::AbstractVector{<:MPolyElem}, M::Union{Matrix{T}, MatElem{T}}) where T
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), MatrixOrdering(i, fmpz_mat(M)))
end

@doc Markdown.doc"""
    weighted_ordering(v::AbstractVector{<:MPolyElem}, w::Vector{Int})

Defines the (partial) weighted ordering on the variables given with the weight
vector `w`. This is equivalent to a matrix ordering with just one row.
"""
function weighted_ordering(v::AbstractVector{<:MPolyElem}, w::Vector{Int}) where T
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), MatrixOrdering(i, fmpz_mat(1, length(w), w)))
end

function _weight_matrix(nvars::Int, o::MatrixOrdering)
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
    weight_matrix(M::MonomialOrdering)

Return a corresponding weight matrix for the given ordering.
"""
function weight_matrix(M::MonomialOrdering)
  return _weight_matrix(nvars(base_ring(M)), M.o)
end

@doc Markdown.doc"""
    simplify(M::MonomialOrdering) -> MonomialOrdering

Returns a matrix ordering with a unique weight matrix.
"""
function Hecke.simplify(M::MonomialOrdering)
  w = canonical_weight_matrix(M)
  return MonomialOrdering(M.R, MatrixOrdering(1:ncols(w), w))
end

function canonical_weight_matrix(nvars::Int, M::AbsOrdering)
  w = _weight_matrix(nvars, M)
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

@doc Markdown.doc"""
    canonical_weight_matrix(M::MonomialOrdering)

Return the corresponding canonical weight matrix for the given ordering.
"""
function canonical_weight_matrix(M::MonomialOrdering)
  return canonical_weight_matrix(nvars(base_ring(M)), M.o)
end

import Base.==
function ==(M::MonomialOrdering, N::MonomialOrdering)
  return canonical_weight_matrix(M) == canonical_weight_matrix(N)
end

function Base.hash(M::MonomialOrdering, u::UInt)
  return hash(canonical_weight_matrix(M), u)
end

# Return two arrays of indices containing the variables which are > 1 and < 1
# respectively.
function _global_and_local_vars(M::MonomialOrdering)
  globals = Int[]
  locals = Int[]
  n = nvars(base_ring(M))
  d = zeros(Int, n)
  w = weight_matrix(M)
  m = nrows(w)
  @assert n == ncols(w)
  for i in 1:m, j in 1:n
    d[j] == 0 || continue
    if w[i, j] > 0
      d[j] = 1
      push!(globals, j)
    elseif w[i, j] < 0
      d[j] = -1
      push!(locals, j)
    end
  end
  return globals, locals
end

@doc Markdown.doc"""
    is_global(M::MonomialOrdering)

Return `true` if the given ordering is global, i.e. if $1 < x$ for
each variable $x$ in the ring `M.R` for which `M` is defined.
"""
function is_global(M::MonomialOrdering)
  globals, locals = _global_and_local_vars(M)
  return length(globals) == nvars(base_ring(M))
end

@doc Markdown.doc"""
    is_local(M::MonomialOrdering)

Return `true` if the given ordering is local, i.e. if $1 > x$ for
each variable $x$ in the ring `M.R` for which `M` is defined.
"""
function is_local(M::MonomialOrdering)
  globals, locals = _global_and_local_vars(M)
  return length(locals) == nvars(base_ring(M))
end

@doc Markdown.doc"""
    is_mixed(M::MonomialOrdering)

Return `true` if the given ordering is mixed, i.e. if $1 < x_i$ for
a variable $x_i$ and $1 > x_j$ for another variable $x_j$ in the ring `M.R`
for which `M` is defined.
"""
function is_mixed(M::MonomialOrdering)
  globals, locals = _global_and_local_vars(M)
  return !isempty(globals) && !isempty(locals)
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
   if I-i+1 == length(a) #test if variables are consecutive or not.
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

function flat(a::ModOrdering)
   return [a]
end
function flat(a::ModProdOrdering)
   return vcat(flat(a.a), flat(a.b))
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

function _expressify(o::MatrixOrdering, sym) where S
  return Expr(:call, nrows(o.matrix) == 1 ? :weighted_ordering : :matrix_ordering,
                     Expr(:vect, (sym[i] for i in o.vars)...),
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

############ Singular conversions ############

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
  if nrows(o.matrix) == n && !iszero(det(o.matrix))
    Q.last_var += n
    return (true, Singular.ordering_M(o.matrix, false))
  elseif nrows(o.matrix) == 1
    return (true, Singular.ordering_a([Int(o.matrix[1,i]) for i in 1:n]))
  end
  return (false, Q.def)
end

function _try_singular_easy(Q::order_conversion_ctx, o::Orderings.ModOrdering)
  Q.has_c_or_C && return (false, Q.def)
  Q.has_c_or_C = true
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
    return Singular.ordering_M(canonical_weight_matrix(ord))
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
    return sord
  end
end

end  # module Orderings
