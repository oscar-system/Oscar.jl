module Orderings

using Oscar
import Oscar: Ring, MPolyRing, MPolyRingElem, weights, IntegerUnion, base_ring,
       support, matrix

export ModuleOrdering
export MonomialOrdering
export anti_diagonal
export canonical_matrix
export deglex
export degrevlex
export index_of_leading_term
export induce
export induced_ring_ordering
export deginvlex
export invlex
export is_elimination_ordering
export is_global
export is_local
export is_mixed
export is_total
export _is_weighted
export lex
export matrix_ordering
export monomial_ordering
export negdeglex
export negdegrevlex
export neginvlex
export neglex
export negwdeglex
export negwdegrevlex
export opposite_ordering
export permutation_of_terms
export singular
export wdeglex
export wdegrevlex
export weight_ordering

abstract type AbsOrdering end

abstract type AbsGenOrdering <: AbsOrdering end

abstract type AbsModOrdering <: AbsOrdering end

# lex, deglex, ...
struct SymbOrdering{S} <: AbsGenOrdering
  vars::Vector{Int}
  function SymbOrdering(S::Symbol, v)
    S in (:lex, :deglex, :degrevlex, :invlex, :revlex, :deginvlex,
          :neglex, :negdeglex, :negdegrevlex, :neginvlex, :negrevlex) ||
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
    @req length(v) == length(w) "number of variables should match the number of weights"
    @req all(>(0), v) "all weights should be positive"
    return new{S}(v, w)
  end
end

# in general denotes a partial ordering on vars
# if fullrank is true, then the matrix should be invertible
struct MatrixOrdering <: AbsGenOrdering
  vars::Vector{Int}
  matrix::ZZMatrix
  fullrank::Bool
  function MatrixOrdering(v, m::ZZMatrix, fullrank::Bool)
    @req length(v) == ncols(m) "number of variables should match the number of columns"
    if fullrank
      if nrows(m) > ncols(m)
        m = _canonical_matrix(m)
      end
      @req nrows(m) >= ncols(m) "weight matrix is rank deficient"
      @assert nrows(m) == ncols(m)
      @req !iszero(det(m)) "weight matrix is not invertible"
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
    nw = w[i:i, :]
    c = content(nw)
    if !isone(c)
      nw = divexact(nw, c)
    end
    for j in 1:nrows(ww)
      h = findfirst(x->ww[j, x] != 0, 1:ncols(w))
      if !iszero(nw[1, h])
        nw = abs(ww[j, h])*nw - sign(ww[j, h])*nw[1, h]*ww[j:j, :]
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
  @assert nrows(ww) <= ncols(ww)
  return ww
end


# convert (y,x,z) => (2,1,3) and check uniqueness
function _unique_var_indices(a::AbstractVector{<:MPolyRingElem})
  !isempty(a) || error("need at least one variable")
  z = Int[var_index(i) for i in a]
  allunique(z) || error("variables must be unique")
  return z
end

function _is_weighted(ord::Symbol)
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

@doc raw"""
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
  is_total::Bool
  is_total_is_known::Bool
  canonical_matrix::ZZMatrix

  function MonomialOrdering(R::S, o::AbsGenOrdering) where {S}
    return new{S}(R, o, false, false)
  end

  function MonomialOrdering(R::S, o::AbsGenOrdering, is_total::Bool) where {S}
    return new{S}(R, o, is_total, true)
  end
end

base_ring(a::MonomialOrdering) = a.R

base_ring_type(::Type{MonomialOrdering{S}}) where {S} = S

@doc raw"""
    support(o::MonomialOrdering)

Return the vector of variables on which `o` is defined.
"""
function support(o::MonomialOrdering)
  return [gen(base_ring(o), i) for i in _support_indices(o.o)]
end

@doc raw"""
    *(o1::MonomialOrdering, o2::MonomialOrdering)

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

@doc raw"""
    monomial_ordering(v::AbstractVector{<:MPolyRingElem}, s::Symbol)

Defines an ordering to be applied to the variables in `v`. The symbol `s`
should be one of `:lex`, `:deglex`, `:degrevlex`, `:invlex`, `:deginvlex`, `:neglex`,
`:negdeglex`, `:negdegrevlex`, `:neginvlex`.
"""
function monomial_ordering(v::AbstractVector{<:MPolyRingElem}, s::Symbol)
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(s, i))
end

function monomial_ordering(R::MPolyRing, s::Symbol)
  return MonomialOrdering(R, SymbOrdering(s, 1:nvars(R)))
end

#### lex ####

@doc raw"""
    lex(R::MPolyRing) -> MonomialOrdering

Return the lexicographical ordering on the set of monomials in the variables of `R`.

    lex(V::AbstractVector{<:MPolyRingElem}) -> MonomialOrdering

Given a vector `V` of variables, return the lexicographical ordering on the set of monomials in these variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = lex(R)
lex([w, x, y, z])

julia> canonical_matrix(o1)
[1   0   0   0]
[0   1   0   0]
[0   0   1   0]
[0   0   0   1]

julia> o2 = lex([w, x])
lex([w, x])

julia> o3 = lex(gens(R)[3:4])
lex([y, z])
```
"""
function lex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:lex, 1:nvars(R)))
end

function lex(v::AbstractVector{<:MPolyRingElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:lex, i))
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::SymbOrdering{:lex})
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### deglex ####

@doc raw"""
    deglex(R::MPolyRing) -> MonomialOrdering

Return the degree lexicographical ordering on the set of monomials in the variables of `R`.

    deglex(V::AbstractVector{<:MPolyRingElem}) -> MonomialOrdering
    
Given a vector `V` of variables, return the degree lexicographical ordering on the set of monomials in these variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = deglex(R)
deglex([w, x, y, z])

julia> canonical_matrix(o1)
[1    1    1    1]
[0   -1   -1   -1]
[0    0   -1   -1]
[0    0    0   -1]

julia> o2 = deglex([w, x])
deglex([w, x])

julia> o3 = deglex(gens(R)[3:4])
deglex([y, z])
```
"""
function deglex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:deglex, 1:nvars(R)))
end

function deglex(v::AbstractVector{<:MPolyRingElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:deglex, i))
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::SymbOrdering{:deglex})
  tdk = sum(exponent(f, k, i) for i in o.vars)
  tdl = sum(exponent(g, l, i) for i in o.vars)
  if tdk > tdl
    return 1
  elseif tdk < tdl
    return -1
  end
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### degrevlex ####

@doc raw"""
    degrevlex(R::MPolyRing) -> MonomialOrdering

Return the degree reverse lexicographical ordering on the set of monomials in the variables of `R`.

    degrevlex(V::AbstractVector{<:MPolyRingElem}) -> MonomialOrdering

Given a vector `V` of variables, return the degree reverse lexicographical ordering on the set of monomials in these variables

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = degrevlex(R)
degrevlex([w, x, y, z])

julia> canonical_matrix(o1)
[1    1    1    1]
[0    0    0   -1]
[0    0   -1    0]
[0   -1    0    0]

julia> o2 = degrevlex([w, x])
degrevlex([w, x])

julia> o3 = degrevlex(gens(R)[3:4])
degrevlex([y, z])
```
"""
function degrevlex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:degrevlex, 1:nvars(R)))
end

function degrevlex(v::AbstractVector{<:MPolyRingElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:degrevlex, i))
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::SymbOrdering{:degrevlex})
  tdk = sum(exponent(f, k, i) for i in o.vars)
  tdl = sum(exponent(g, l, i) for i in o.vars)
  if tdk > tdl
    return 1
  elseif tdk < tdl
    return -1
  end
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### invlex ####

@doc raw"""
    invlex(R::MPolyRing) -> MonomialOrdering

Return the inverse lexicographical ordering on the set of monomials in the variables of `R`.

    invlex(V::AbstractVector{<:MPolyRingElem}) -> MonomialOrdering

Given a vector `V` of variables, return the inverse lexicographical ordering on the set of monomials in these variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = invlex(R)
invlex([w, x, y, z])

julia> canonical_matrix(o1)
[0   0   0   1]
[0   0   1   0]
[0   1   0   0]
[1   0   0   0]

julia> o2 = invlex([w, x])
invlex([w, x])

julia> o3 = invlex(gens(R)[3:4])
invlex([y, z])
```
"""
function invlex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:invlex, 1:nvars(R)))
end

function invlex(v::AbstractVector{<:MPolyRingElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:invlex, i))
end


function _matrix(nvars::Int, o::SymbOrdering{:invlex})
  m = zero_matrix(ZZ, length(o.vars), nvars)
  i = length(o.vars)
  for j in o.vars
    m[i, j] = 1
    i -= 1
  end
  return m
end

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::SymbOrdering{:invlex})
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### deginvlex ####

@doc raw"""
    deginvlex(R::MPolyRing) -> MonomialOrdering

Return the degree inverse lexicographical ordering on the set of monomials in the variables of `R`.

    deginvlex(V::AbstractVector{<:MPolyRingElem}) -> MonomialOrdering

Given a vector `V` of variables, return the degree inverse lexicographical ordering on the set of monomials in these variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = deginvlex(R)
deginvlex([w, x, y, z])

julia> canonical_matrix(o1)
[1   1   1   1]
[0   0   0   1]
[0   0   1   0]
[0   1   0   0]

julia> o2 = deginvlex([w, x])
deginvlex([w, x])

julia> o3 = deginvlex(gens(R)[3:4])
deginvlex([y, z])
```
"""
function deginvlex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:deginvlex, 1:nvars(R)))
end

function deginvlex(v::AbstractVector{<:MPolyRingElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:deginvlex, i))
end

function _matrix(nvars::Int, o::SymbOrdering{:deginvlex})
  m = zero_matrix(ZZ, 1 + length(o.vars), nvars)
  for j in o.vars
    m[1, j] = 1
  end
  i = 1 + length(o.vars)
  for j in o.vars
    m[i, j] = 1
    i -= 1
  end
  return m
end

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::SymbOrdering{:deginvlex})
  tdk = sum(exponent(f, k, i) for i in o.vars)
  tdl = sum(exponent(g, l, i) for i in o.vars)
  if tdk > tdl
    return 1
  elseif tdk < tdl
    return -1
  end
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### neglex ####

@doc raw"""
    neglex(R::MPolyRing) -> MonomialOrdering

Return the negative lexicographical ordering on the set of monomials in the variables of `R`.

    neglex(V::AbstractVector{<:MPolyRingElem}) -> MonomialOrdering

Given a vector `V` of variables, return the negative lexicographical ordering on the set of monomials in these variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = neglex(R)
neglex([w, x, y, z])

julia> canonical_matrix(o1)
[-1    0    0    0]
[ 0   -1    0    0]
[ 0    0   -1    0]
[ 0    0    0   -1]

julia> o2 = neglex([w, x])
neglex([w, x])

julia> o3 = neglex(gens(R)[3:4])
neglex([y, z])
```
"""
function neglex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:neglex, 1:nvars(R)))
end

function neglex(v::AbstractVector{<:MPolyRingElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:neglex, i))
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::SymbOrdering{:neglex})
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### neginvlex ####

@doc raw"""
    neginvlex(R::MPolyRing) -> MonomialOrdering

Return the negative inverse lexicographical ordering on the set of monomials in the variables of `R`.

    neginvlex(V::AbstractVector{<:MPolyRingElem}) -> MonomialOrdering

Given a vector `V` of variables, return the negative inverse lexicographical ordering on the set of monomials in these variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = neginvlex(R)
neginvlex([w, x, y, z])

julia> canonical_matrix(o1)
[ 0    0    0   -1]
[ 0    0   -1    0]
[ 0   -1    0    0]
[-1    0    0    0]

julia> o2 = neginvlex([w, x])
neginvlex([w, x])

julia> o3 = neginvlex(gens(R)[3:4])
neginvlex([y, z])
```
"""
function neginvlex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:neginvlex, 1:nvars(R)))
end

function neginvlex(v::AbstractVector{<:MPolyRingElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:neginvlex, i))
end

function _matrix(nvars::Int, o::SymbOrdering{:neginvlex})
  m = zero_matrix(ZZ, length(o.vars), nvars)
  i = length(o.vars)
  for j in o.vars
    m[i, j] = -1
    i -= 1
  end
  return m
end

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::SymbOrdering{:neginvlex})
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### negdegrevlex ####

@doc raw"""
    negdegrevlex(R::MPolyRing) -> MonomialOrdering

Return the negative degree reverse lexicographical ordering on the set of monomials in the variables of `R`.

    negdegrevlex(V::AbstractVector{<:MPolyRingElem}) -> MonomialOrdering

Given a vector `V` of variables, return the negative degree reverse lexicographical ordering on the set of monomials in these variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = negdegrevlex(R)
negdegrevlex([w, x, y, z])

julia> canonical_matrix(o1)
[-1   -1   -1   -1]
[ 0    0    0   -1]
[ 0    0   -1    0]
[ 0   -1    0    0]

julia> o2 = negdegrevlex([w, x])
negdegrevlex([w, x])

julia> o3 = negdegrevlex(gens(R)[3:4])
negdegrevlex([y, z])
```
"""
function negdegrevlex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:negdegrevlex, 1:nvars(R)))
end

function negdegrevlex(v::AbstractVector{<:MPolyRingElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:negdegrevlex, i))
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::SymbOrdering{:negdegrevlex})
  tdk = sum(exponent(f, k, i) for i in o.vars)
  tdl = sum(exponent(g, l, i) for i in o.vars)
  if tdk > tdl
    return -1
  elseif tdk < tdl
    return 1
  end
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### negdeglex ####

@doc raw"""
    negdeglex(R::MPolyRing) -> MonomialOrdering

Return the negative degree lexicographical ordering on the set of monomials in the variables of `R`.

    negdeglex(V::AbstractVector{<:MPolyRingElem}) -> MonomialOrdering    

Given a vector `V` of variables, return the negative degree lexicographical ordering on the set of monomials in these variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = negdeglex(R)
negdeglex([w, x, y, z])

julia> canonical_matrix(o1)
[-1   -1   -1   -1]
[ 0   -1   -1   -1]
[ 0    0   -1   -1]
[ 0    0    0   -1]

julia> o2 = negdeglex([w, x])
negdeglex([w, x])

julia> o3 = negdeglex(gens(R)[3:4])
negdeglex([y, z])
```
"""
function negdeglex(R::MPolyRing)
  return MonomialOrdering(R, SymbOrdering(:negdeglex, 1:nvars(R)))
end

function negdeglex(v::AbstractVector{<:MPolyRingElem})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), SymbOrdering(:negdeglex, i))
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::SymbOrdering{:negdeglex})
  tdk = sum(exponent(f, k, i) for i in o.vars)
  tdl = sum(exponent(g, l, i) for i in o.vars)
  if tdk > tdl
    return -1
  elseif tdk < tdl
    return 1
  end
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
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

@doc raw"""
    monomial_ordering(v::AbstractVector{<:MPolyRingElem}, s::Symbol, w::Vector{Int})
    monomial_ordering(R::MPolyRing, s::Symbol, w::Vector{Int}) -> MonomialOrdering

Defines a weighted ordering to be applied to the variables in `v`. The weight
vector `w` should be the same length as `v`, and the symbol `s` should be one
of `:wdeglex`, `:wdegrevlex`, `:negwdeglex`, `:negwdegrevlex`.
"""
function monomial_ordering(v::AbstractVector{<:MPolyRingElem}, s::Symbol, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(s, i, w))
end

function monomial_ordering(R::MPolyRing, s::Symbol, w::Vector{Int})
  return MonomialOrdering(R, WSymbOrdering(s, 1:nvars(R), w))
end

function _cmp_weighted_degree(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, vars::Vector{Int}, w::Vector{Int})
  n = length(vars)
  @assert n == length(w)
  t = 0
  for j in 1:n
    v = vars[j]
    t += w[j]*(exponent(f, k, v) - exponent(g, l, v))
  end
  return t > 0 ? 1 : t < 0 ? -1 : 0;
end

#### wdeglex, Wp ####

@doc raw"""
    wdeglex(R::MPolyRing, W::Vector{Int}) -> MonomialOrdering

If `W` is a vector of positive integers, return the corresponding weighted 
lexicographical ordering on the set of monomials in the variables of `R`.

    wdeglex(V::AbstractVector{<:MPolyRingElem}, W::Vector{Int}) -> MonomialOrdering
    
Given a vector `V` of variables and a vector `W` of positive integers, return the corresponding weighted 
lexicographical ordering on the set of monomials in the given variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = wdeglex(R, [1, 2, 3, 4])
wdeglex([w, x, y, z], [1, 2, 3, 4])

julia> canonical_matrix(o1)
[1    2    3    4]
[0   -2   -3   -4]
[0    0   -3   -4]
[0    0    0   -1]

julia> o2 = wdeglex([w, x], [1, 2])
wdeglex([w, x], [1, 2])

julia> o3 = wdeglex(gens(R)[3:4], [3, 4])
wdeglex([y, z], [3, 4])
```
"""
function wdeglex(R::MPolyRing, w::Vector{Int})
  return MonomialOrdering(R, WSymbOrdering(:wdeglex, 1:nvars(R), w))
end

function wdeglex(v::AbstractVector{<:MPolyRingElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:wdeglex, i, w))
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::WSymbOrdering{:wdeglex})
  c = _cmp_weighted_degree(f, k, g, l, o.vars, o.weights)
  c == 0 || return c
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### wdegrevlex, wp ####

@doc raw"""
    wdegrevlex(R::MPolyRing, W::Vector{Int}) -> MonomialOrdering

If `W` is a vector of positive integers, return the corresponding weighted reverse 
lexicographical ordering on the set of monomials in the variables of `R`.

    wdegrevlex(V::AbstractVector{<:MPolyRingElem}, W::Vector{Int}) -> MonomialOrdering

Given a vector `V` of variables and a vector `W` of positive integers, return the corresponding weighted reverse 
lexicographical ordering on the set of monomials in the given variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = wdegrevlex(R, [1, 2, 3, 4])
wdegrevlex([w, x, y, z], [1, 2, 3, 4])

julia> canonical_matrix(o1)
[1    2    3    4]
[0    0    0   -1]
[0    0   -1    0]
[0   -1    0    0]

julia> o2 = wdegrevlex([w, x], [1, 2])
wdegrevlex([w, x], [1, 2])

julia> o3 = wdegrevlex(gens(R)[3:4], [3, 4])
wdegrevlex([y, z], [3, 4])
```
"""
function wdegrevlex(R::MPolyRing, w::Vector{Int})
  return MonomialOrdering(R, WSymbOrdering(:wdegrevlex, 1:nvars(R), w))
end

function wdegrevlex(v::AbstractVector{<:MPolyRingElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:wdegrevlex, i, w))
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::WSymbOrdering{:wdegrevlex})
  c = _cmp_weighted_degree(f, k, g, l, o.vars, o.weights)
  c == 0 || return c
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return -1
    elseif ek < el
      return 1
    end
  end
  return 0
end

#### negwdeglex, Ws ####

@doc raw"""
    negwdeglex(R::MPolyRing, W::Vector{Int}) -> MonomialOrdering

If `W` is a vector of positive integers, return the corresponding negative weighted 
lexicographical ordering on the set of monomials in the variables of `R`.

    negwdeglex(V::AbstractVector{<:MPolyRingElem}, W::Vector{Int}) -> MonomialOrdering
    
Given a vector `V` of variables and a vector `W` of positive integers, return the corresponding
negative weighted lexicographical ordering on the set of monomials in the given variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = negwdeglex(R, [1, 2, 3, 4])
negwdeglex([w, x, y, z], [1, 2, 3, 4])

julia> canonical_matrix(o1)
[-1   -2   -3   -4]
[ 0   -2   -3   -4]
[ 0    0   -3   -4]
[ 0    0    0   -1]

julia> o2 = negwdeglex([w, x], [1, 2])
negwdeglex([w, x], [1, 2])

julia> o3 = negwdeglex(gens(R)[3:4], [3, 4])
negwdeglex([y, z], [3, 4])
```
"""
function negwdeglex(R::MPolyRing, w::Vector{Int})
  return MonomialOrdering(R, WSymbOrdering(:negwdeglex, 1:nvars(R), w))
end

function negwdeglex(v::AbstractVector{<:MPolyRingElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:negwdeglex, i, w))
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::WSymbOrdering{:negwdeglex})
  c = _cmp_weighted_degree(g, l, f, k, o.vars, o.weights)
  c == 0 || return c
  for i in o.vars
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
    if ek > el
      return 1
    elseif ek < el
      return -1
    end
  end
  return 0
end

#### negwdegrevlex, ws ####

@doc raw"""
    negwdegrevlex(R::MPolyRing, W::Vector{Int}) -> MonomialOrdering

If `W` is a vector of positive integers, return the corresponding negative weighted
reverse lexicographical ordering on the set of monomials in the variables of `R`.

    negwdegrevlex(V::AbstractVector{<:MPolyRingElem}, W::Vector{Int}) -> MonomialOrdering
    
Given a vector `V` of variables and a vector `W` of positive integers, return the corresponding negative
weighted reverse lexicographical ordering on the set of monomials in the given variables.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[w, x, y, z])

julia> o1 = negwdegrevlex(R, [1, 2, 3, 4])
negwdegrevlex([w, x, y, z], [1, 2, 3, 4])

julia> canonical_matrix(o1)
[-1   -2   -3   -4]
[ 0    0    0   -1]
[ 0    0   -1    0]
[ 0   -1    0    0]

julia> o2 = negwdegrevlex([w, x], [1, 2])
negwdegrevlex([w, x], [1, 2])

julia> o3 = negwdegrevlex(gens(R)[3:4], [3, 4])
negwdegrevlex([y, z], [3, 4])
```
"""
function negwdegrevlex(R::MPolyRing, w::Vector{Int})
  return MonomialOrdering(R, WSymbOrdering(:negwdegrevlex, 1:nvars(R), w))
end

function negwdegrevlex(v::AbstractVector{<:MPolyRingElem}, w::Vector{Int})
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), WSymbOrdering(:negwdegrevlex, i, w))
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::WSymbOrdering{:negwdegrevlex})
  c = _cmp_weighted_degree(g, l, f, k, o.vars, o.weights)
  c == 0 || return c
  for i in reverse(o.vars)
    ek = exponent(f, k, i)
    el = exponent(g, l, i)
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

@doc raw"""
    matrix_ordering(R::MPolyRing, M::Union{Matrix{T}, MatElem{T}}; check::Bool = true) where T -> MonomialOrdering

Given an integer matrix `M` such that `nvars(R) = ncols(M) = rank(M)`, 
return the matrix ordering on the set of variables of `R` which is defined by `M`.

    matrix_ordering(V::AbstractVector{<:MPolyRingElem}, M::Union{Matrix{T}, MatElem{T}}; check::Bool = true) where T -> MonomialOrdering

Given a vector `V` of variables and an integer matrix `M` such that `length(V) = ncols(M) = rank(M)`, 
return the matrix ordering on the set of monomials in the given variables which is defined by `M`.

!!! note
    The matrix `M` need not be square.

!!! note
    If `check = false` is supplied, the rank check is omitted, and the resulting
    ordering may only be partial.

# Examples
```jldoctest
julia> R, (w, x, y, z) = QQ["w", "x", "y", "z"];

julia> M =[1 1 1 1; 0 -1 -1 -1; 0 0 -1 -1; 0 0 0 -1]
4×4 Matrix{Int64}:
 1   1   1   1
 0  -1  -1  -1
 0   0  -1  -1
 0   0   0  -1

julia> o1 = matrix_ordering(R, M)
matrix_ordering([w, x, y, z], [1 1 1 1; 0 -1 -1 -1; 0 0 -1 -1; 0 0 0 -1])

julia> N =[1 1; 0 -1]
2×2 Matrix{Int64}:
 1   1
 0  -1

julia> o2 = matrix_ordering([w, x], N)
matrix_ordering([w, x], [1 1; 0 -1])

julia> canonical_matrix(o2)
[1    1   0   0]
[0   -1   0   0]

julia> o3 = matrix_ordering(gens(R)[3:4], N)
matrix_ordering([y, z], [1 1; 0 -1])

julia> o3 = matrix_ordering(gens(R)[3:4], N)
matrix_ordering([y, z], [1 1; 0 -1])

julia> canonical_matrix(o3)
[0   0   1    1]
[0   0   0   -1]
```
"""
function matrix_ordering(R::MPolyRing, M::Union{Matrix{T}, MatElem{T}}; check::Bool = true) where T
  return MonomialOrdering(R, MatrixOrdering(1:nvars(R), ZZMatrix(M), check))
end

function matrix_ordering(v::AbstractVector{<:MPolyRingElem}, M::Union{Matrix{T}, MatElem{T}}; check::Bool = true) where T
  i = _unique_var_indices(v)
  return MonomialOrdering(parent(first(v)), MatrixOrdering(i, ZZMatrix(M), check))
end

@doc raw"""
    weight_ordering(W::Vector{Int}, ord::MonomialOrdering) -> MonomialOrdering

Given an integer vector `W` and a monomial ordering `ord` on a set of monomials in
`length(W)` variables, return the monomial ordering `ord_W` on this set of monomials
which is obtained by first comparing the `W`-weighted degrees and then using `ord` 
in the case of a tie.

!!! note 
    The ordering `ord_W` is   
    - global if all entries of `W` are positive, or if they are all non-negative and `ord` is global,
    - an elimination ordering for the set of variables which correspond to positive entries of `W`.   

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> W = [1, 0, -1];

julia> o = lex(R)
lex([x, y, z])

julia> matrix(o)
[1   0   0]
[0   1   0]
[0   0   1]

julia> oW = weight_ordering(W, o)
matrix_ordering([x, y, z], [1 0 -1])*lex([x, y, z])

julia> matrix(oW)
[1   0   -1]
[1   0    0]
[0   1    0]
[0   0    1]

julia> canonical_matrix(oW)
[1   0   -1]
[0   0    1]
[0   1    0]

julia> o2 = weight_ordering([1, -1], lex([x, z]))
matrix_ordering([x, z], [1 -1])*lex([x, z])

julia> canonical_matrix(o2)
[1   0   -1]
[0   0    1]
```
"""
function weight_ordering(w::Vector{Int}, o::MonomialOrdering)
  i = _support_indices(o.o)
  m = ZZMatrix(1, length(w), w)
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

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::MatrixOrdering)
  M = o.matrix
  r = nrows(M)
  c = ncols(M)
  for i in 1:r
    v = o.vars[1]
    t = M[i,1]*(exponent(f, k, v) - exponent(g, l, v))
    for j in 2:c
      v = o.vars[j]
      t += M[i, j]*(exponent(f, k, v) - exponent(g, l, v))
    end
    if t > 0
      return 1
    elseif t < 0
      return -1
    end
  end
  return 0 
end

function _cmp_monomials(f::MPolyRingElem, k::Int, g::MPolyRingElem, l::Int, o::ProdOrdering)
  c = _cmp_monomials(f, k, g, l, o.a)
  c == 0 || return c
  return _cmp_monomials(f, k, g, l, o.b)
end

###################################

@doc raw"""
    matrix(ord::MonomialOrdering)

Return a matrix defining `ord` as a matrix ordering.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> o1 = degrevlex(R)
degrevlex([x, y, z])

julia> matrix(o1)
[ 1    1    1]
[ 0    0   -1]
[ 0   -1    0]
[-1    0    0]

julia> o2 = degrevlex([x, y])
degrevlex([x, y])

julia> matrix(o2)
[ 1    1   0]
[ 0   -1   0]
[-1    0   0]
```
"""
function matrix(M::MonomialOrdering)
  return _matrix(nvars(base_ring(M)), M.o)
end

@doc raw"""
    simplify(M::MonomialOrdering) -> MonomialOrdering

Return a matrix ordering with a unique weight matrix.
"""
function Hecke.simplify(M::MonomialOrdering)
  w = canonical_matrix(M)
  return MonomialOrdering(M.R, MatrixOrdering(collect(1:ncols(w)), w, true))
end

function canonical_matrix(nvars::Int, M::AbsOrdering)
  return _canonical_matrix(_matrix(nvars, M))
end

@doc raw"""
    canonical_matrix(ord::MonomialOrdering)

Return the canonical matrix defining `ord` as a matrix ordering.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> o1 = degrevlex(R)
degrevlex([x, y, z])

julia> canonical_matrix(o1)
[1    1    1]
[0    0   -1]
[0   -1    0]

julia> o2 = degrevlex([x, y])
degrevlex([x, y])

julia> canonical_matrix(o2)
[1    1   0]
[0   -1   0]
```
"""
function canonical_matrix(M::MonomialOrdering)
  if !isdefined(M, :canonical_matrix)
    M.canonical_matrix = canonical_matrix(nvars(base_ring(M)), M.o)
  end
  return M.canonical_matrix
end

# is_obviously_equal checks whether the two orderings are defined in the exact
# same way, so `true` means that the orderings are equal and `false` means
# that they might be equal but not 'obviously' so.

# Generic fall back
function is_obviously_equal(M::AbsOrdering, N::AbsOrdering)
  return false
end

function is_obviously_equal(M::SymbOrdering{S}, N::SymbOrdering{S}) where {S}
  return M.vars == N.vars
end

function is_obviously_equal(M::WSymbOrdering{S}, N::WSymbOrdering{S}) where {S}
  return M.vars == N.vars && M.weights == N.weights
end

function is_obviously_equal(M::MatrixOrdering, N::MatrixOrdering)
  return M.vars == N.vars && M.matrix == N.matrix
end

function is_obviously_equal(M::ProdOrdering, N::ProdOrdering)
  return is_obviously_equal(M.a, N.a) && is_obviously_equal(M.b, N.b)
end

import Base.==
function ==(M::MonomialOrdering, N::MonomialOrdering)
  if M === N || is_obviously_equal(M.o, N.o)
    return true
  end
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

@doc raw"""
    is_global(ord::MonomialOrdering)

Return `true` if `ord` is global, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> o = matrix_ordering(R, [1 1; 0 -1])
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

@doc raw"""
    is_local(ord::MonomialOrdering)

Return `true` if `ord` is local, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> o = matrix_ordering(R, [-1 -1; 0 -1])
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

@doc raw"""
    is_mixed(ord::MonomialOrdering)

Return `true` if `ord` is mixed, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> o = matrix_ordering(R, [1 -1; 0 -1])
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

@doc raw"""
    is_total(ord::MonomialOrdering)

Return `true` if `ord` is total ordering, `false` otherwise.
"""
function is_total(ord::MonomialOrdering)
  if !ord.is_total_is_known
    m = canonical_matrix(ord)
    ord.is_total = nrows(m) == ncols(m)
    ord.is_total_is_known = true
  end
  return ord.is_total
end

@doc raw"""
    cmp(ord::MonomialOrdering, a::MPolyRingElem, b::MPolyRingElem)

    cmp(ord::ModuleOrdering, a::FreeModElem{T}, b::FreeModElem{T}) where T <: MPolyRingElem

Compare monomials `a` and `b` with regard to the ordering `ord`: Return `-1` for `a < b`
and `1` for `a > b` and `0` for `a == b`. An error is thrown if `ord` is
a partial ordering that does not distinguish `a` from `b`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> cmp(lex([x,y]), x, one(R))
1

julia> try cmp(lex([x,y]), z, one(R)); catch e; e; end
ErrorException("z and 1 are incomparable with respect to lex([x, y])")

julia> cmp(lex([x,y,z]), z, one(R))
1

julia> F = free_module(R, 2)
Free module of rank 2 over R

julia> cmp(lex(R)*invlex(F), F[1], F[2])
-1
```
"""
function Base.cmp(ord::MonomialOrdering, a::MPolyRingElem, b::MPolyRingElem)
  @assert base_ring(ord) === parent(a) === parent(b)
  @assert is_monomial(a)
  @assert is_monomial(b)
  c = _cmp_monomials(a, 1, b, 1, ord.o)
  if c == 0 && a != b
    error("$a and $b are incomparable with respect to $ord")
  end
  return c
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

@doc raw"""
    is_elimination_ordering(ord::MonomialOrdering, V::Vector{<:MPolyRingElem})

Given a vector `V` of polynomials which are variables, return `true` if `ord` is an elimination ordering for the variables in `V`.
Return `false`, otherwise.

    is_elimination_ordering(ord::MonomialOrdering, V:Vector{Int})
    
Given a vector `V` of indices which specify variables, return `true` if `ord` is an elimination ordering for the specified variables.
Return `false`, otherwise.

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> o1 = lex(R)
lex([w, x, y, z])

julia> is_elimination_ordering(o1, [w, x])
true

julia> o2 = weight_ordering([1, 1, 0, 0], degrevlex(R))
matrix_ordering([w, x, y, z], [1 1 0 0])*degrevlex([w, x, y, z])

julia> is_elimination_ordering(o2, [w, x])
true

julia> o3 = weight_ordering([1, -1, 0, 0], degrevlex(R))
matrix_ordering([w, x, y, z], [1 -1 0 0])*degrevlex([w, x, y, z])

julia> is_elimination_ordering(o3, [w, x])
false

julia> o4 = degrevlex([w, x])*degrevlex([y, z])
degrevlex([w, x])*degrevlex([y, z])

julia> is_elimination_ordering(o4, [w, x])
true

julia> o5 = degrevlex([w, x])*negdegrevlex([y, z])
degrevlex([w, x])*negdegrevlex([y, z])

julia> is_elimination_ordering(o5, [w, x])
true

julia> o6 = negdegrevlex([w, x])*negdegrevlex([y, z])
negdegrevlex([w, x])*negdegrevlex([y, z])

julia> is_elimination_ordering(o6, [w, x])
false
```
"""
function is_elimination_ordering(o::MonomialOrdering, sigmaC::Vector{<:MPolyRingElem})
  return is_elimination_ordering(o, Int[var_index(v) for v in sigmaC])
end

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

function induce(vars::Vector{Int}, o::SymbOrdering{S}) where S
  return SymbOrdering(S, vars)
end

function induce(vars::Vector{Int}, o::WSymbOrdering{S}) where S
  @req length(vars) == length(o.vars) "Number of variables must coincide"
  return WSymbOrdering(S, vars, o.weights)
end

function induce(vars::Vector{Int}, o::MatrixOrdering)
  @req length(vars) == length(o.vars) "Number of variables must coincide"
  return MatrixOrdering(vars, o.matrix, o.fullrank)
end

function induce(vars::Vector{Int}, o::ProdOrdering)
  m = length(_support_indices(o.a))
  n = length(_support_indices(o.b))
  @req m + n == length(vars) "Number of variables must coincide"
  return induce(vars[1:m], o.a)*induce(vars[m + 1:end], o.b)
end

@doc raw"""
    induce(vars::AbstractVector{<:MPolyRingElem}, ord::ModuleOrdering)

Return the monomial ordering on the variables `vars` induced by transferring
the ordering `ord`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> S, (a, b, c) = polynomial_ring(GF(5), ["a", "b", "c"]);

julia> ord = degrevlex([x, y])*neglex([z]);

julia> induce([a, b, c], ord)
degrevlex([a, b])*neglex([c])
```
"""
function induce(vars::AbstractVector{<: MPolyRingElem}, o::MonomialOrdering)
  @assert !isempty(vars)
  v = _unique_var_indices(vars)
  o2 = induce(v, o.o)
  return MonomialOrdering(parent(vars[1]), o2)
end

###################################################

# Module orderings (not module Orderings)

# For ordering the generators (I think)
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
  o::AbsModOrdering # must allow gen*mon or mon*gen product ordering

  canonical_matrix::ZZMatrix

  function ModuleOrdering(M::S, o::AbsModOrdering) where {S}
    return new{S}(M, o)
  end
end

base_ring(a::ModuleOrdering) = a.M

base_ring_type(::Type{ModuleOrdering{S}}) where {S} = S

mutable struct ModProdOrdering <: AbsModOrdering
  a::AbsOrdering
  b::AbsOrdering
end

Base.:*(a::AbsGenOrdering, b::AbsModOrdering) = ModProdOrdering(a, b)

Base.:*(a::AbsModOrdering, b::AbsGenOrdering) = ModProdOrdering(a, b)

# For equality checking and hashing
# we produce a matrix representation by embedding the underlying free module (and its ordering) into a polynomial ring
# then we build the matrix in this ring

function _embedded_ring_ordering(o::ModuleOrdering)
  return _embedded_ring_ordering(o.o)
end

function _embedded_ring_ordering(o::ModOrdering)
  return SymbOrdering(o.ord, o.gens)
end

function _embedded_ring_ordering(o::AbsGenOrdering)
  return deepcopy(o)
end

function _embedded_ring_ordering(o::ModProdOrdering)
  ea = _embedded_ring_ordering(o.a)
  shift = maximum(ea.vars)
  eb = _embedded_ring_ordering(o.b)
  eb.vars .+= shift 
  return ea*eb
end

function _canonical_matrix_intern(o::ModuleOrdering)
  nvrs = ngens(o.M) + ngens(base_ring(o.M))
  eo = _embedded_ring_ordering(o)
  return canonical_matrix(nvrs, eo)
end

function canonical_matrix(o::ModuleOrdering)
  if !isdefined(o, :canonical_matrix)
    if isone(ngens(o.M))
      o.canonical_matrix = canonical_matrix(induced_ring_ordering(o))
    else
      o.canonical_matrix = _canonical_matrix_intern(o)
    end
  end
  return o.canonical_matrix
end

# is_obviously_equal checks whether the two orderings are defined in the exact
# same way, so `true` means that the orderings are equal and `false` means
# that they might be equal but not 'obviously' so.
# The generic fall back is_obviously_equal(::AbsOrdering, ::AbsOrdering) = false
# is defined above.

function is_obviously_equal(o1::ModOrdering{T}, o2::ModOrdering{T}) where {T}
  return o1.gens == o2.gens && o1.ord == o2.ord
end

function is_obviously_equal(o1::ModProdOrdering, o2::ModProdOrdering)
  return is_obviously_equal(o1.a, o2.a) && is_obviously_equal(o1.b, o2.b)
end

function Base.:(==)(o1::ModuleOrdering, o2::ModuleOrdering)
  if o1 === o2 || is_obviously_equal(o1.o, o2.o)
    return true
  end
  return canonical_matrix(o1) == canonical_matrix(o2)
end

function Base.hash(o::ModuleOrdering, h::UInt)
  return hash(canonical_matrix(o), h)
end

#### _cmp_vector_monomials: cmp f[k]*gen(m) with g[l]*gen(n)

function _cmp_vector_monomials(
  m::Int, f::MPolyRingElem, k::Int,
  n::Int, g::MPolyRingElem, l::Int,
  o::ModOrdering)

  if o.ord == :lex
    return m < n ? 1 : m > n ? -1 : 0
  elseif o.ord == :invlex
    return m > n ? 1 : m < n ? -1 : 0
  else
    error("oops")
  end
end

function _cmp_vector_monomials(
  m::Int, f::MPolyRingElem, k::Int,
  n::Int, g::MPolyRingElem, l::Int,
  o::AbsGenOrdering)

  return _cmp_monomials(f, k, g, l, o)
end

function _cmp_vector_monomials(
  m::Int, f::MPolyRingElem, k::Int,
  n::Int, g::MPolyRingElem, l::Int,
  o::ModProdOrdering)

  c = _cmp_vector_monomials(m, f, k, n, g, l, o.a)
  c == 0 || return c
  return _cmp_vector_monomials(m, f, k, n, g, l, o.b)
end


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

function invlex(v::AbstractVector{<:AbstractAlgebra.ModuleElem})
   return ModuleOrdering(parent(first(v)), ordering(v, :invlex))
end

function Base.:*(M::ModuleOrdering, N::MonomialOrdering)
   base_ring(M.M) == N.R || error("wrong rings")
   return ModuleOrdering(M.M, M.o*N.o)
end

function Base.:*(M::MonomialOrdering, N::ModuleOrdering)
   base_ring(N.M) == M.R || error("wrong rings")
   return ModuleOrdering(N.M, M.o*N.o)
end

@doc raw"""
    induced_ring_ordering(ord::ModuleOrdering)

Return the ring ordering induced by `ord`.  

# Examples
```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> F = free_module(R, 3)
Free module of rank 3 over R

julia> o = invlex(gens(F))*degrevlex(R)
invlex([gen(1), gen(2), gen(3)])*degrevlex([w, x, y, z])

julia> induced_ring_ordering(o)
degrevlex([w, x, y, z])
```
"""
function induced_ring_ordering(M::ModuleOrdering)
  return MonomialOrdering(base_ring(M.M), induced_ring_ordering(M.o))
end

function induced_ring_ordering(o::ModOrdering)
  error("no induced ring ordering exists")
end

function induced_ring_ordering(o::AbsGenOrdering)
  return o
end

function induced_ring_ordering(o::ModProdOrdering)
  if o.a isa ModOrdering
    return induced_ring_ordering(o.b)
  elseif o.b isa ModOrdering
    return induced_ring_ordering(o.a)
  else
    return induced_ring_ordering(o.a)*induced_ring_ordering(o.b)
  end
end

@doc raw"""
    is_global(M::ModuleOrdering)

Return `true` if the given ordering is global, i.e. if the
induced ring ordering is global.
"""
function is_global(M::ModuleOrdering)
  return is_global(induced_ring_ordering(M))
end


@doc raw"""
    permutation_of_terms(f::MPolyRingElem, ord::MonomialOrdering)

Return the permutation that puts the terms of `f` in the order `ord`.
"""
function permutation_of_terms(f::MPolyRingElem, ord::MonomialOrdering)
  p = collect(1:length(f))
  sort!(p, lt = (k, l) -> (Orderings._cmp_monomials(f, k, f, l, ord.o) < 0), rev = true)
  return p
end

@doc raw"""
    index_of_leading_term(f::MPolyRingElem, ord::MonomialOrdering)

Return the index of the leading term of `f` with repsect to the order `ord`. The
returned index is itself relative to the ordering in `AbstractAlgebra.terms(f)`.
"""
function index_of_leading_term(f::MPolyRingElem, ord::MonomialOrdering)
  n = length(f)
  @req n > 0 "zero polynomial does not have a leading term"
  res = 1
  for i in 2:n
    if Orderings._cmp_monomials(f, res, f, i, ord.o) < 0
      res = i
    end
  end
  return res
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

# TODO return more singular-friendly orderings if that is how they come in
function _opposite_ordering(nvars::Int, o::SymbOrdering{T}) where T
  return SymbOrdering(T, nvars+1 .- o.vars)
end

function _opposite_ordering(nvars::Int, o::SymbOrdering{:deglex})
  n = length(o.vars)
  newvars = reverse(nvars+1 .- o.vars)
  return MatrixOrdering(newvars, ZZMatrix(1, n, ones(Int, n)), false)*
         SymbOrdering(:invlex, newvars)
end

# TODO ditto
function _opposite_ordering(nvars::Int, o::WSymbOrdering{T}) where T
  return WSymbOrdering(T, nvars+1 .- o.vars, o.weights)
end

function _opposite_ordering(nvars::Int, o::MatrixOrdering)
  M = o.matrix
  M = reduce(hcat, [M[:,i:i] for i in ncols(M):-1:1])
  return MatrixOrdering(reverse(nvars+1 .- o.vars), M, false)
end

function _opposite_ordering(n::Int, o::ProdOrdering)
  return ProdOrdering(_opposite_ordering(n, o.a), _opposite_ordering(n, o.b))
end

@doc raw"""
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
         S == :invlex       ? (true, Singular.ordering_ip(n)) :
         S == :deginvlex    ? (true, Singular.ordering_Ip(n)) :
         S == :deglex       ? (true, Singular.ordering_Dp(n)) :
         S == :degrevlex    ? (true, Singular.ordering_dp(n)) :
         S == :deginvlex    ? (true, Singular.ordering_Ip(n)) :
         S == :neglex       ? (true, Singular.ordering_ls(n)) :
         S == :neginvlex    ? (true, Singular.ordering_is(n)) :
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
         o.ord == :invlex ? (true, Singular.ordering_c(length(o.gens))) :
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
  elseif o.order == Singular.ringorder_ip
    return SymbOrdering(:invlex, i), newlastvar
  elseif o.order == Singular.ringorder_Dp
    return SymbOrdering(:deglex, i), newlastvar
  elseif o.order == Singular.ringorder_Ip
    return SymbOrdering(:deginvlex, i), newlastvar
  elseif o.order == Singular.ringorder_dp
    return SymbOrdering(:degrevlex, i), newlastvar
  elseif o.order == Singular.ringorder_ls
    return SymbOrdering(:neglex, i), newlastvar
  elseif o.order == Singular.ringorder_is
    return SymbOrdering(:neginvlex, i), newlastvar
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
    m = ZZMatrix(1, length(o.weights), o.weights)
    return MatrixOrdering(i, m, false), lastvar
  elseif o.order == Singular.ringorder_M
    m = ZZMatrix(o.size, o.size, o.weights)
    return MatrixOrdering(i, m, true), newlastvar

  elseif o.order == Singular.ringorder_C
    return nothing, lastvar
  elseif o.order == Singular.ringorder_c
    return nothing, lastvar

  else
    error("cannot convert singular block $(o.order)")
  end
end

@doc raw"""
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

import Oscar.Orderings: induce # needed at least for group characters

