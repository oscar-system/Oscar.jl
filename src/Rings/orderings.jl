module Orderings

using Oscar, Markdown
import Oscar: Ring, MPolyRing, MPolyElem, weights, IntegerUnion
export anti_diagonal, lex, degrevlex, deglex, revlex, negdeglex,
       neglex, negrevlex, negdegrevlex, wdeglex, wdegrevlex,
       negwdeglex, negwdegrevlex, matrix_ordering, weights, isweighted,
       MonomialOrdering, ModuleOrdering, singular,
       isglobal, islocal, ismixed

abstract type AbsOrdering end

abstract type AbsGenOrdering <: AbsOrdering end

abstract type AbsModOrdering <: AbsOrdering end

"""
Ring-free monomial ordering: just the indices of the variables are given.
`T` can be a `UnitRange` to make Singular happy or any `Array` if the
variables are not consecutive
"""
mutable struct GenOrdering{T} <: AbsGenOrdering
  vars::T
  ord::Symbol
  wgt::fmpz_mat
  w::Vector{Int}
  function GenOrdering(u::T, s::Symbol) where {T <: AbstractVector{Int}}
    r = new{typeof(u)}()
    r.vars = u
    r.ord = s
    return r
  end
  function GenOrdering(u::T, m::fmpz_mat; ord::Symbol = :weight) where {T <: AbstractVector{Int}}
    r = new{typeof(u)}()
    @assert ncols(m) == length(u)
    r.vars = u
    r.ord = ord
    r.wgt = m
    return r
  end
  function GenOrdering(u::T, w::Vector{Int}; ord::Symbol) where {T <: AbstractVector{Int}}
   r = new{typeof(u)}()
   r.vars = u
   r.ord = ord
   r.w = w
   return r
 end
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

#not really user facing
function ordering(a::AbstractVector{Int}, s::Union{Symbol, fmpz_mat})
  i = minimum(a)
  I = maximum(a)
  if I-i+1 == length(a) #test if variables are consecutive or not.
    return GenOrdering(i:I, s)
  end
  return GenOrdering(collect(a), s)
end

#not really user facing
function ordering(a::AbstractVector{Int}, s::Symbol, w::fmpz_mat)
  i = minimum(a)
  I = maximum(a)
  if I-i+1 == length(a)
    return GenOrdering(i:I, w, ord = s)
  end
  return GenOrdering(collect(a), w, ord = s)
end

#not really user facing
function ordering(a::AbstractVector{Int}, s::Symbol, w::Vector{Int})
   i = minimum(a)
   I = maximum(a)
   length(a) == length(w) || error("number of weights doesn't match number of variables")
   if I-i+1 == length(a)
     return GenOrdering(i:I, w, ord = s)
   end
   return GenOrdering(collect(a), w, ord = s)
 end
 

#not really user facing, flattens a product of product orderings into an array 
function flat(a::GenOrdering)
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

#not user facing
function weights(a::GenOrdering)
  if a.ord == :lex || a.ord == Symbol("Singular(lp)")
    return identity_matrix(ZZ, length(a.vars))
  end
  if a.ord == :deglex || a.ord == Symbol("Singular(Dp)")
    return [matrix(ZZ, 1, length(a.vars), ones(fmpz, length(a.vars)));
            identity_matrix(ZZ, length(a.vars)-1) zero_matrix(ZZ, length(a.vars)-1, 1)]
  end
  if a.ord == :degrevlex || a.ord == Symbol("Singular(dp)")
    return [matrix(ZZ, 1, length(a.vars), ones(fmpz, length(a.vars))) ;
            zero_matrix(ZZ, length(a.vars)-1, 1) -anti_diagonal(ZZ, length(a.vars)-1)]
  end
  if a.ord == :revlex || a.ord == Symbol("Singular(rp)")
    return anti_diagonal(ZZ, length(a.vars))
  end
  if a.ord == :wdeglex || a.ord == Symbol("Singular(Wp)")
    return [matrix(ZZ, 1, length(a.vars), a.w);
            identity_matrix(ZZ, length(a.vars)-1) zero_matrix(ZZ, length(a.vars)-1, 1)]
  end
  if a.ord == :wdegrevlex || a.ord == Symbol("Singular(wp)")
    return [matrix(ZZ, 1, length(a.vars), a.w);
            zero_matrix(ZZ, length(a.vars)-1, 1) anti_diagonal(ZZ, length(a.vars)-1)]
  end
  if a.ord == :neglex || a.ord == Symbol("Singular(ls)")
    return -identity_matrix(ZZ, length(a.vars))
  end
  if a.ord == :negdeglex || a.ord == Symbol("Singular(Ds)")
    return [-matrix(ZZ, 1, length(a.vars), ones(fmpz, length(a.vars)));
            -identity_matrix(ZZ, length(a.vars)-1) zero_matrix(ZZ, length(a.vars)-1, 1)]
  end
  if a.ord == :negdegrevlex || a.ord == Symbol("Singular(ds)")
    return [-matrix(ZZ, 1, length(a.vars), ones(fmpz, length(a.vars))) ;
            zero_matrix(ZZ, length(a.vars)-1, 1) -anti_diagonal(ZZ, length(a.vars)-1)]
  end
  if a.ord == :negrevlex || a.ord == Symbol("Singular(rs)")
    return -anti_diagonal(ZZ, length(a.vars))
  end
  if a.ord == :negwdeglex || a.ord == Symbol("Singular(Ws)")
    return [-matrix(ZZ, 1, length(a.vars), a.w);
            -identity_matrix(ZZ, length(a.vars)-1) zero_matrix(ZZ, length(a.vars)-1, 1)]
  end
  if a.ord == :negwdegrevlex || a.ord == Symbol("Singular(ws)")
    return [-matrix(ZZ, 1, length(a.vars), a.w);
            zero_matrix(ZZ, length(a.vars)-1, 1) -anti_diagonal(ZZ, length(a.vars)-1)]
  end
  if a.ord == :weight || a.ord == Symbol("Singular(a)") || a.ord == Symbol("Singular(M)")
    return a.wgt
  end
  error("Ordering not recognized")
end

#not user facing
function weights(a::AbsOrdering)
  aa = flat(a)
  m = matrix(ZZ, 0, 0, [])
  for o = aa
    if o isa GenOrdering
      w = weights(o)
      if maximum(o.vars) > ncols(m)
        m = hcat(m, zero_matrix(ZZ, nrows(m), maximum(o.vars) - ncols(m)))
      end
      mm = zero_matrix(ZZ, nrows(w), ncols(m))
      for r = 1:nrows(w)
        for c = 1:length(o.vars)
          mm[r, o.vars[c]] = w[r, c]
        end
      end
      m = vcat(m, mm)
    end
  end
  return m
end

"""
Orderings actually applied to polynomial rings (as opposed to variable indices)
"""
mutable struct MonomialOrdering{S}
  R::S
  o::AbsGenOrdering
end

#not really user facing, not exported
@doc Markdown.doc"""
    ordering(a::Vector{MPolyElem}, s::Symbol)
    ordering(a::Vector{MPolyElem}, m::fmpz_mat)
    ordering(a::Vector{MPolyElem}, s::Symbol, m::fmpz_mat)
    ordering(a::Vector{MPolyElem}, s::Symbol, w::Vector{Int})

Defines an ordering to be applied to the variables in `a`.
In the first form the symbol `s` has to be one of `:lex`, `:deglex`,
`:degrevlex`, `:revlex`, `:neglex`, `:negdeglex`, `:negdegrevlex`,
`:negrevlex`, `:wdeglex`, `:wdegrevlex`, `:negwdeglex`, `:negwdegrevlex`.
In the second form, a weight ordering using the given matrix is used.
In the last version, the symbol if of the form `Singular(..)`.
"""
function ordering(a::AbstractVector{<:MPolyElem}, s...)
  R = parent(first(a))
  g = gens(R)
  aa = [findfirst(x -> x == y, g) for y = a]
  if nothing in aa
    error("only variables allowed")
  end
  return ordering(aa, s...)
end

@doc Markdown.doc"""
    :*(M::MonomialOrdering, N::MonomialOrdering)

For orderings on the same ring, the product ordering obtained by concatenation
of the weight matrix.
"""
function Base.:*(M::MonomialOrdering, N::MonomialOrdering)
  M.R == N.R || error("wrong rings")
  return MonomialOrdering(M.R, M.o*N.o)
end

function Base.show(io::IO, M::MonomialOrdering)
  a = flat(M.o)
  if length(a) > 1
    print(io, "Product ordering: ")
    for i=1:length(a)-1
      show(io, M.R, a[i])
      print(io, " * ")
    end
  end
  show(io, M.R, a[end])
end

function Base.show(io::IO, R::MPolyRing, o::AbsOrdering)
  if o.ord == :weight
    print(io, "weight($(gens(R)[o.vars]) via $(o.wgt))")
  elseif isweighted(o.ord)
   print(io, "weight($(gens(R)[o.vars]) via $(o.w))")
  else
    print(io, "$(String(o.ord))($(gens(R)[o.vars]))")
  end
end

@doc Markdown.doc"""
    lex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `lex` (lexicographic) ordering on the variables given.
"""
function lex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :lex))
end
@doc Markdown.doc"""
    deglex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `deglex` ordering on the variables given.
"""
function deglex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :deglex))
end
@doc Markdown.doc"""
    degrevlex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `degrevlex` ordering on the variables given.
"""
function degrevlex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :degrevlex))
end
@doc Markdown.doc"""
    revlex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `revlex` ordering on the variables given.
"""
function revlex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :revlex))
end
@doc Markdown.doc"""
    neglex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `neglex` ordering on the variables given.
"""
function neglex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :neglex))
end
@doc Markdown.doc"""
    negrevlex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `negrevlex` ordering on the variables given.
"""
function negrevlex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :negrevlex))
end
@doc Markdown.doc"""
    negdegrevlex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `negdegrevlex` ordering on the variables given.
"""
function negdegrevlex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :negdegrevlex))
end
@doc Markdown.doc"""
    negdeglex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `negdeglex` ordering on the variables given.
"""
function negdeglex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :negdeglex))
end
@doc Markdown.doc"""
    wdeglex(v::AbstractVector{<:MPolyElem}, w::Vector{Int}) -> MonomialOrdering

Defines the `wdeglex` ordering on the variables given with the weights `w`.
"""
function wdeglex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  return MonomialOrdering(parent(first(v)), ordering(v, :wdeglex, w))
end
@doc Markdown.doc"""
    wdegrevlex(v::AbstractVector{<:MPolyElem}, w::Vector{Int}) -> MonomialOrdering

Defines the `wdegrevlex` ordering on the variables given with the weights `w`.
"""
function wdegrevlex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  return MonomialOrdering(parent(first(v)), ordering(v, :wdegrevlex, w))
end
@doc Markdown.doc"""
    negwdeglex(v::AbstractVector{<:MPolyElem}, w::Vector{Int}) -> MonomialOrdering

Defines the `negwdeglex` ordering on the variables given with the weights `w`.
"""
function negwdeglex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  return MonomialOrdering(parent(first(v)), ordering(v, :negwdeglex, w))
end
@doc Markdown.doc"""
    negwdegrevlex(v::AbstractVector{<:MPolyElem}, w::Vector{Int}) -> MonomialOrdering

Defines the `negwdegrevlex` ordering on the variables given with the weights `w`.
"""
function negwdegrevlex(v::AbstractVector{<:MPolyElem}, w::Vector{Int})
  return MonomialOrdering(parent(first(v)), ordering(v, :negwdegrevlex, w))
end
@doc Markdown.doc"""
    matrix_ordering(v::AbstractVector{<:MPolyElem}, M::fmpz_mat) -> MonomialOrdering

Defines the matrix ordering on the variables given with the matrix `M`.
"""
function matrix_ordering(v::AbstractVector{<:MPolyElem}, M::fmpz_mat)
  return MonomialOrdering(parent(first(v)), ordering(v, M))
end

@doc Markdown.doc"""
    singular(ord::Symbol, v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines an ordering given in terms of Singular primitives on the variables given.
`ord` can be one of `:lp`, `:rs`, `:ls`, `:dp`, `:rp`, `:ds`, `:Ds`, `:Dp`.
"""
function singular(ord::Symbol, v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, Symbol("Singular($(string(ord)))")))
end

@doc Markdown.doc"""
    singular(ord::Symbol, v::AbstractVector{<:MPolyElem}, w::AbstractMatrix{<:IntegerUnion}) -> MonomialOrdering

Defines an ordering given in terms of Singular weight ordering (`M`) with the
matrix given. `ord` has to be `:M` here.
"""
function singular(ord::Symbol, v::AbstractVector{<:MPolyElem}, w::AbstractMatrix{<:IntegerUnion})
  @assert ord == :M
  W = matrix(ZZ, size(w, 1), size(w, 2), w)
  return MonomialOrdering(parent(first(v)), ordering(v, Symbol("Singular($(string(ord)))"), W))
end

@doc Markdown.doc"""
    singular(ord::Symbol, v::AbstractVector{<:MPolyElem}, w::AbstractVector{<:IntegerUnion}) -> MonomialOrdering

Defines an ordering given in terms of Singular weight ordering (`a`) with the
weights given. `ord` has to be `:a` here. The weights will be supplemented by
`0`.
"""
function singular(ord::Symbol, v::AbstractVector{<:MPolyElem}, w::AbstractVector{<:IntegerUnion})
  @assert ord == :a
  W = map(fmpz, w)
  while length(v) > length(W)
    push!(W, 0)
  end
  return MonomialOrdering(parent(first(v)), ordering(v, Symbol("Singular($(string(ord)))"), matrix(ZZ, 1, length(W), W)))
end

@doc Markdown.doc"""
    weights(M::MonomialOrdering)
 
Compute a corresponding weight matrix for the given ordering.
"""
function weights(M::MonomialOrdering)
  return weights(M.o)
end

@doc Markdown.doc"""
    simplify(M::MonomialOrdering) -> MonomialOrdering

Compute a weight ordering with a unique weight matrix.    
"""
function Hecke.simplify(M::MonomialOrdering)
  ww = simplify_weight_matrix(M.o)
  return MonomialOrdering(M.R, ordering(1:ncols(ww), ww))
end

function simplify_weight_matrix(M::AbsOrdering)
  w = weights(M)
  ww = matrix(ZZ, 0, ncols(w), [])
  for i=1:nrows(w)
    if is_zero_row(w, i)
      continue
    end
    nw = w[i, :]
    c = content(nw)
    if c != 1
      nw = divexact(nw, c)
    end
    for j=1:nrows(ww)
      h = findfirst(x->ww[j, x] != 0, 1:ncols(w))
      if nw[1, h] != 0
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

import Base.==
function ==(M::MonomialOrdering, N::MonomialOrdering)
  return Hecke.simplify(M).o.wgt == Hecke.simplify(N).o.wgt
end

function Base.hash(M::MonomialOrdering, u::UInt)
  return hash(Hecke.simplify(M).o.wgt, u)
end

@doc Markdown.doc"""
    isglobal(M::MonomialOrdering)

Return `true` if the given ordering is global, i.e. if $1 < x$ for
each variable $x$ in the ring `M.R` for which `M` is defined.
"""
function isglobal(M::MonomialOrdering)
  fl = flat(M.o)
  R = M.R
  for o in fl
    for i in o.vars
      if isone(leading_monomial(gens(R)[i] + one(R), M))
        return false
      end
    end
  end
  return true
end

@doc Markdown.doc"""
    islocal(M::MonomialOrdering)

Return `true` if the given ordering is local, i.e. if $1 > x$ for
each variable $x$ in the ring `M.R` for which `M` is defined.
"""
function islocal(M::MonomialOrdering)
  fl = flat(M.o)
  R = M.R
  for o in fl
    for i in o.vars
      if !isone(leading_monomial(gens(R)[i] + one(R), M))
        return false
      end
    end
  end
  return true
end

@doc Markdown.doc"""
    ismixed(M::MonomialOrdering)

Return `true` if the given ordering is mixed, i.e. if $1 < x_i$ for
a variable $x_i$ and $1 > x_j$ for another variable $x_j$ in the ring `M.R`
for which `M` is defined.
"""
ismixed(M::MonomialOrdering) = !isglobal(M) && !islocal(M)

# Return two arrays of indices containing the variables which are > 1 and < 1
# respectively.
function _global_and_local_vars(M::MonomialOrdering)
  globals = Int[]
  locals = Int[]
  for i in 1:ngens(M.R)
    if isone(leading_monomial(gens(M.R)[i] + 1, M))
      push!(locals, i)
    else
      push!(globals, i)
    end
  end
  return globals, locals
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

mutable struct ModProdOrdering <: AbsModOrdering
   a::AbsOrdering
   b::AbsOrdering
end

function Base.show(io::IO, M::ModuleOrdering)
   a = flat(M.o)
   if length(a) > 1
     print(io, "Product ordering: ")
     for i=1:length(a)-1
       if isa(a[i], GenOrdering)
         show(io, base_ring(M.M), a[i])
       else
         show(io, M.M, a[i])
       end
       print(io, " * ")
     end
   end
   if isa(a[end], GenOrdering)
      show(io, base_ring(M.M), a[end])
   else
      show(io, M.M, a[end])
   end
end

function Base.show(io::IO, R::AbstractAlgebra.Module, o::AbsOrdering)
   if o.ord == :lex
      print(io, "Module ordering (lex)")
   elseif o.ord == :revlex
      print(io, "Module ordering (revlex)")
   else
      error("module ordering unknown")
   end
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

function max_used_variable(st::Int, o::GenOrdering)
   g = o.vars
   mi = minimum(g)
   ma = maximum(g)
   if mi == st && length(g) + st == ma+1
      return ma+1
   else
      return 0
   end
end

function max_used_variable(st::Int, o::ModOrdering)
   return -1
end

end  # module Orderings

###################################################

# Some isless functions for orderings given by a symbol:
# _isless_:ord(f, k, l) returns true if the k-th term is lower than the l-th
# term of f in the ordering :ord.

 # lp: weight matrix [1 0 0; 0 1 0; 0 0 1]
 function _isless_lex(f::MPolyElem, k::Int, l::Int)
   n = nvars(parent(f))
   for i = 1:n
     ek = exponent(f, k, i)
     el = exponent(f, l, i)
     if ek != el
       return ek < el
     end
   end
   return false
 end
 
 # ls: weight matrix [-1 0 0; 0 -1 0; 0 0 -1]   !!! confusing !!!
 function _isless_neglex(f::MPolyElem, k::Int, l::Int)
   n = nvars(parent(f))
   for i = 1:n
     ek = exponent(f, k, i)
     el = exponent(f, l, i)
     if ek != el
       return ek > el
     end
   end
   return false
 end
 
 # rp: weight matrix [0 0 1; 0 1 0; 1 0 0]   !!! confusing !!!
 function _isless_revlex(f::MPolyElem, k::Int, l::Int)
   n = nvars(parent(f))
   for i = n:-1:1
     ek = exponent(f, k, i)
     el = exponent(f, l, i)
     if ek != el
        return ek < el
     end
   end
   return false
 end
 
 # rs: weight matrix [0 0 -1; 0 -1 0; -1 0 0]   !!! confusing !!!
 function _isless_negrevlex(f::MPolyElem, k::Int, l::Int)
   n = nvars(parent(f))
   for i = n:-1:1
     ek = exponent(f, k, i)
     el = exponent(f, l, i)
     if ek != el
       return ek > el
     end
   end
   return false
 end
 
 # Dp: weight matrix [1 1 1; 1 0 0; 0 1 0; 0 0 1]
 function _isless_deglex(f::MPolyElem, k::Int, l::Int)
   tdk = total_degree(term(f, k))
   tdl = total_degree(term(f, l))
   if tdk < tdl
     return true
   elseif tdk > tdl
     return false
   end
   return _isless_lex(f, k, l)
 end
 
 # dp: weight matrix [1 1 1; 0 0 -1; 0 -1 0; -1 0 0]
 function _isless_degrevlex(f::MPolyElem, k::Int, l::Int)
   tdk = total_degree(term(f, k))
   tdl = total_degree(term(f, l))
   if tdk != tdl
     return tdk < tdl
   end
   return _isless_negrevlex(f, k, l)
 end
 
 # Ds: weight matrix [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1]
 function _isless_negdeglex(f::MPolyElem, k::Int, l::Int)
   tdk = total_degree(term(f, k))
   tdl = total_degree(term(f, l))
   if tdk != tdl
     return tdk > tdl
   end
   return _isless_lex(f, k, l)
 end
 
 # ds: weight matrix [-1 -1 -1; 0 0 -1; 0 -1 0; -1 0 0]
 function _isless_negdegrevlex(f::MPolyElem, k::Int, l::Int)
   tdk = total_degree(term(f, k))
   tdl = total_degree(term(f, l))
   if tdk != tdl
     return tdk > tdl
   end
   return _isless_negrevlex(f, k, l)
 end
 
 # Returns the degree of the k-th term of f weighted by w,
 # that is deg(x^a) = w_1a_1 + ... + w_na_n.
 # No sanity checks are performed!
 function weighted_degree(f::MPolyElem, k::Int, w::Vector{Int})
   ek = exponent_vector(f, k)
   return dot(ek, w)
 end

 function weighted_degree(f::MPolyElem, k::Int, r::UnitRange{Int}, w::Vector{Int})
   ek = exponent_vector(f, k)[r]
   return dot(ek, w)
 end

 # Wp: the matrix is [w; 1 0 0; 0 1 0; 0 0 1]
 function _isless_weightdeglex(f::MPolyElem, k::Int, l::Int, w::Vector{Int})
   dk = weighted_degree(f, k, w)
   dl = weighted_degree(f, l, w)
   if dk != dl
     return dk < dl
   end
   return _isless_lex(f, k, l)
 end

 # wp: the matrix is [w; 0 0 -1; 0 -1 0; -1 0 0]
 function _isless_weightdegrevlex(f::MPolyElem, k::Int, l::Int, w::Vector{Int})
   dk = weighted_degree(f, k, w)
   dl = weighted_degree(f, l, w)
   if dk != dl
     return dk < dl
   end
   return _isless_negrevlex(f, k, l)
 end
 
 # Ws: the matrix is [-w; 1 0 0; 0 1 0; 0 0 1]
 function _isless_weightnegdeglex(f::MPolyElem, k::Int, l::Int, w::Vector{Int})
   dk = weighted_degree(f, k, w)
   dl = weighted_degree(f, l, w)
   if dk != dl
     return dk > dl
   end
   return _isless_lex(f, k, l)
 end

 # ws: the matrix is [-w; 0 0 -1; 0 -1 0; -1 0 0]
 function _isless_weightnegdegrevlex(f::MPolyElem, k::Int, l::Int, w::Vector{Int})
   dk = weighted_degree(f, k, w)
   dl = weighted_degree(f, l, w)
   if dk != dl
     return dk > dl
   end
   return _isless_negrevlex(f, k, l)
 end
 
 function _isless_matrix(f::MPolyElem, k::Int, l::Int, M::Union{ Matrix{T}, MatElem{T} }) where T
   ek = exponent_vector(f, k)
   el = exponent_vector(f, l)
   n = nvars(parent(f))
   for i = 1:size(M, 1)
     eki = sum( M[i, j]*ek[j] for j = 1:n )
     eli = sum( M[i, j]*el[j] for j = 1:n )
     if eki == eli
       continue
     elseif eki > eli
       return false
     else
       return true
     end
   end
   return false
 end
 
 function _perm_of_terms(f::MPolyElem, ord_lt::Function)
   p = collect(1:length(f))
   sort!(p, lt = (k, l) -> ord_lt(f, k, l), rev = true)
   return p
 end
 
 # Requiring R for consistence with the other lt_from_ordering functions
 function lt_from_ordering(::MPolyRing, ord::Symbol)
   if ord == :lex || ord == :lp
     return _isless_lex
   elseif ord == :revlex || ord == :rp
     return _isless_revlex
   elseif ord == :deglex || ord == :Dp
     return _isless_deglex
   elseif ord == :degrevlex || ord == :dp
     return _isless_degrevlex
   elseif ord == :neglex || ord == :ls
     return _isless_neglex
   elseif ord == :negrevlex || ord == :rs
     return _isless_negrevlex
   elseif ord == :negdeglex || ord == :Ds
     return _isless_negdeglex
   elseif ord == :negdegrevlex || ord == :ds
     return _isless_negdegrevlex
   else
     error("Ordering $ord not available")
   end
 end
 
 function lt_from_ordering(R::MPolyRing, ord::Symbol, w::Vector{Int})
   @assert length(w) == nvars(R) "Number of weights has to match number of variables"
 
   if ord == :wdeglex || ord == :Wp
     @assert all(x -> x > 0, w) "Weights have to be positive"
     return (f, k, l) -> _isless_weightdeglex(f, k, l, w)
   elseif ord == :wdegrevlex || ord == :wp
     @assert all(x -> x > 0, w) "Weights have to be positive"
     return (f, k, l) -> _isless_weightdegrevlex(f, k, l, w)
   elseif ord == :negwdeglex || ord == :Ws
     @assert !iszero(w[1]) "First weight must not be 0"
     return (f, k, l) -> _isless_weightnegdeglex(f, k, l, w)
   elseif ord == :negwdegrevlex || ord == :ws
     @assert !iszero(w[1]) "First weight must not be 0"
     return (f, k, l) -> _isless_weightnegdegrevlex(f, k, l, w)
   else
     error("Ordering $ord not available")
   end
 end
 
 function lt_from_ordering(R::MPolyRing, M::Union{ Matrix{T}, MatElem{T} }) where T
   @assert size(M, 2) == nvars(R) "Matrix dimensions have to match number of variables"
 
   return (f, k, l) -> _isless_matrix(f, k, l, M)
 end

###################################################

# Comparisons for orderings given by an ordering
# object

# Return -1, 0, +1 if term k of f is <, ==, > term
# l of f wrt ord
# not user facing
function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::Oscar.Orderings.GenOrdering{T}) where T <: UnitRange{Int}
   tdk = sum(exponent(f, k, i) for i in o.vars)
   tdl = sum(exponent(f, l, i) for i in o.vars)
   if o.ord == :lex
      for i in o.vars
         ek = exponent(f, k, i)
         el = exponent(f, l, i)
         if ek > el
           return 1
         elseif ek < el
           return -1
         end
       end
   elseif o.ord == :deglex
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
   elseif o.ord == :revlex
      for i in reverse(o.vars)
         ek = exponent(f, k, i)
         el = exponent(f, l, i)
         if ek > el
           return 1
         elseif ek < el
           return -1
         end
       end
   elseif o.ord == :degrevlex
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
   elseif o.ord == :neglex
      for i in o.vars
         ek = exponent(f, k, i)
         el = exponent(f, l, i)
         if ek > el
           return -1
         elseif ek < el
           return 1
         end
       end
   elseif o.ord == :negdeglex
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
   elseif o.ord == :negrevlex
      for i in reverse(o.vars)
         ek = exponent(f, k, i)
         el = exponent(f, l, i)
         if ek > el
           return -1
         elseif ek < el
           return 1
         end
       end
   elseif o.ord == :negdegrevlex
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
   elseif o.ord == :weight
     if _isless_matrix(f, k, l, o.wgt)
       return -1
     elseif _isless_matrix(f, l, k, o.wgt)
       return 1
     end
   else
      dk = weighted_degree(f, k, o.vars, o.w)
      dl = weighted_degree(f, l, o.vars, o.w)
      if o.ord == :wdeglex
         if dk < dl
            return -1
          elseif dk > dl
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
      elseif o.ord == :wdegrevlex
         if dk < dl
            return -1
          elseif dk > dl
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
      elseif o.ord == :negwdeglex
         if dk > dl
            return -1
          elseif dk < dl
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
      elseif o.ord == :negwdegrevlex
         if dk > dl
            return -1
          elseif dk < dl
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
      else
         error("ordering not implemented")
      end
   end
   return 0 # equal wrt this (partial) ordering
end

# fallback for cases where T is not a UnitRange{Int}
# variables may be out of order
# just compare using weight matrix
function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::Oscar.Orderings.GenOrdering{T}) where T
   if _isless_matrix(f, k, l, weights(o))
      return -1
   elseif _isless_matrix(f, l, k, weights(o))
      return 1
   end
   return 0 # equal wrt this (partial) ordering
end

function _cmp_monomials(f::MPolyElem, k::Int, l::Int, o::Oscar.Orderings.ProdOrdering)
   V = Oscar.Orderings.flat(o)
   for oi in V
      cmp = _cmp_monomials(f, k, l, oi)
      if cmp != 0
         return cmp
      end
   end
   return 0
end

function lt_from_ordering(R::MPolyRing, o::Oscar.Orderings.AbsGenOrdering)
   return (f, k, l) -> _cmp_monomials(f, k, l, o) < 0
end

