export WordOrdering, degree_left_lex, degree_right_lex,
       left_elimination_ordering, right_elimination_ordering,
       weighted_degree_left_lex, weighted_degree_right_lex

###############################################################################
#
#  Word orderings
#
###############################################################################

abstract type AbsWordOrdering end

struct SymbWordOrdering{S} <: AbsWordOrdering
end

struct WSymbWordOrdering{S} <: AbsWordOrdering
  weights::Vector{Int}
end

mutable struct WordOrdering{T}
  base_ring::T
  o::AbsWordOrdering
end


# left lex is not a real ordering
function _cmp_words(
  o::SymbWordOrdering{:_left_lex},
  a::Vector{Int}, b::Vector{Int}
)
  an = length(a)
  bn = length(b)
  for i in 1:min(an,bn)
    if a[i] != b[i]
      return a[i] < b[i] ? 1 : -1
    end
  end
  return an < bn ? 1 : -1
end

# right lex is not a real ordering
function _cmp_words(
  o::SymbWordOrdering{:_right_lex},
  a::Vector{Int}, b::Vector{Int}
)
  an = length(a)
  bn = length(b)
  for i in 1:min(an,bn)
    if a[1+an-i] != b[1+bn-i]
      return a[1+an-i] < b[1+bn-i] ? 1 : -1
    end
  end
  return an < bn ? 1 : -1
end

function _word_weight(a::Vector{Int}, w::Vector{Int})
  s = 0
  for ai in a
    s += w[ai]
  end
  return s
end

function _min_and_max(a::Vector{Int})
  amin = typemax(Int)
  amax = typemin(Int)
  for ai in a
    amin = min(amin, ai)
    amax = max(amax, ai)
  end
  return amin, amax
end

####

@doc Markdown.doc"""
    degree_left_lex(R::FreeAssAlgebra)

The ordering which compares the length of words first, then uses left
lexicographic as a tie breaker.
"""
function degree_left_lex(R::FreeAssAlgebra)
  return WordOrdering(R, SymbWordOrdering{:degree_left_lex}())
end

function _cmp_words(
  o::SymbWordOrdering{:degree_left_lex},
  a::Vector{Int}, b::Vector{Int}
)
  if length(a) != length(b)
    return length(a) > length(b) ? 1 : -1
  end
  return _cmp_words(SymbWordOrdering{:_left_lex}(), a, b)
end

####

@doc Markdown.doc"""
    degree_right_lex(R::FreeAssAlgebra)

The ordering which compares the length of words first, then uses right
lexicographic as a tie breaker.
"""
function degree_right_lex(R::FreeAssAlgebra)
  return WordOrdering(R, SymbWordOrdering{:degree_right_lex}())
end

function _cmp_words(
  o::SymbWordOrdering{:degree_right_lex},
  a::Vector{Int}, b::Vector{Int}
)
  if length(a) != length(b)
    return length(a) > length(b) ? 1 : -1
  end
  return _cmp_words(SymbWordOrdering{:_right_lex}(), a, b)
end

####

@doc Markdown.doc"""
    left_elimination_ordering(R::FreeAssAlgebra)

The left elimination ordering on, say, monomials in $K<x,y,z>$ first compares
the counts of $x$, then in case of a tie compares the counts of $y$, then in
case of a tie compares the counts of $z$. The final tie breaker is left
lexicographic with $x>y>z$.
"""
function left_elimination_ordering(R::FreeAssAlgebra)
  return WordOrdering(R, SymbWordOrdering{:left_elimination_ordering}())
end

function _cmp_words(
  o::SymbWordOrdering{:left_elimination_ordering},
  a::Vector{Int}, b::Vector{Int}
)
  an = length(a)
  bn = length(b)

  if an == 0
    return bn == 0 ? 0 : -1
  elseif bn == 0
    return 1
  end

  amin, amax = _min_and_max(a)
  bmin, bmax = _min_and_max(b)

  if amin != bmin
    return amin < bmin ? 1 : -1
  end

  for i in amin:min(amax, bmax)
    ca = count(==(i), a)
    cb = count(==(i), b)
    if ca != cb
      return ca > cb ? 1 : -1
    end
  end

  if amax != bmax
    return amax > bmax ? 1 : -1
  end

  @assert length(a) == length(b)
  return _cmp_words(SymbWordOrdering{:_left_lex}(), a, b)
end

####

@doc Markdown.doc"""
    right_elimination_ordering(R::FreeAssAlgebra)

The right elimination ordering on, say, monomials in $K<x,y,z>$ first compares
the counts of $z$, then in case of a tie compares the counts of $y$, then in
case of a tie compares the counts of $x$. The final tie breaker is right
lexicographic with $z>y>x$.
"""
function right_elimination_ordering(R::FreeAssAlgebra)
  return WordOrdering(R, SymbWordOrdering{:right_elimination_ordering}())
end

function _cmp_words(
  o::SymbWordOrdering{:right_elimination_ordering},
  a::Vector{Int}, b::Vector{Int}
)
  an = length(a)
  bn = length(b)

  if an == 0
    return bn == 0 ? 0 : -1
  elseif bn == 0
    return 1
  end

  amin, amax = _min_and_max(a)
  bmin, bmax = _min_and_max(b)

  if amax != bmax
    return amax > bmax ? 1 : -1
  end

  for i in amax:-1:max(amin, bmin)
    ca = count(==(i), a)
    cb = count(==(i), b)
    if ca != cb
      return ca > cb ? 1 : -1
    end
  end

  if amin != bmin
    return amin < bmin ? 1 : -1
  end

  @assert length(a) == length(b)  # so we can just use -cmp(right_lex)
  return -_cmp_words(SymbWordOrdering{:_right_lex}(), a, b)
end

####

@doc Markdown.doc"""
    weighted_degree_left_lex(R::FreeAssAlgebra, w::Vector{Int})

The ordering which compares the weighted degree of words first, then uses left
lexicographic as a tie breaker. The weights must be strictly positive.
"""
function weighted_degree_left_lex(R::FreeAssAlgebra, w::Vector{Int})
  return WordOrdering(R, WSymbWordOrdering{:weighted_degree_left_lex}(w))
end

function _cmp_words(
  o::WSymbWordOrdering{:weighted_degree_left_lex},
  a::Vector{Int}, b::Vector{Int}
)
  wa = _word_weight(a, o.weights)
  wb = _word_weight(b, o.weights)
  if wa != wb
    return wa > wb ? 1 : -1
  end
  return _cmp_words(SymbWordOrdering{:_left_lex}(), a, b)
end

####

@doc Markdown.doc"""
    weighted_degree_right_lex(R::FreeAssAlgebra, w::Vector{Int})

The ordering which compares the weighted degree of words first, then uses right
lexicographic as a tie breaker. The weights must be strictly positive.
"""
function weighted_degree_right_lex(R::FreeAssAlgebra, w::Vector{Int})
  return WordOrdering(R, WSymbWordOrdering{:weighted_degree_right_lex}(w))
end

function _cmp_words(
  o::WSymbWordOrdering{:weighted_degree_right_lex},
  a::Vector{Int}, b::Vector{Int}
)
  wa = _word_weight(a, o.weights)
  wb = _word_weight(b, o.weights)
  if wa != wb
    return wa > wb ? 1 : -1
  end
  return _cmp_words(SymbWordOrdering{:_right_lex}(), a, b)
end

###############################################################################
#
#  Terms, Monomials, ect. in a specific ordering
#
###############################################################################

# don't modify the return of this function
function _collected_exponent_words(a::FreeAssAlgElem)
  return collect(exponent_words(a))
end

# nor this one
function _collected_coefficients(a::FreeAssAlgElem)
  return collect(coefficients(a))
end

_collected_exponent_words(a::AbstractAlgebra.Generic.FreeAssAlgElem) = view(a.exps, 1:length(a))
_collected_coefficients(a::AbstractAlgebra.Generic.FreeAssAlgElem) = view(a.coeffs, 1:length(a))

function permutation_of_terms(a::FreeAssAlgElem, o::WordOrdering)
  return sortperm(_collected_exponent_words(a);
                  lt = (a, b) -> (_cmp_words(o.o, a, b) > 0))
end

function coefficients(a::FreeAssAlgElem, o::WordOrdering)
  p = permutation_of_terms(a, o)
  c = _collected_coefficients(a)
  return (c[p[i]] for i in 1:length(p))
end

function exponent_words(a::FreeAssAlgElem, o::WordOrdering)
  p = permutation_of_terms(a, o)
  e = _collected_exponent_words(a)
  return (e[p[i]] for i in 1:length(p))
end

function coefficients_and_exponent_words(a::FreeAssAlgElem, o::WordOrdering)
  p = permutation_of_terms(a, o)
  c = _collected_coefficients(a)
  e = _collected_exponent_words(a)
  return ((c[p[i]], e[p[i]]) for i in 1:length(p))
end

function monomials(a::FreeAssAlgElem, o::WordOrdering)
  p = permutation_of_terms(a, o)
  e = _collected_exponent_words(a)
  R = parent(a)
  return (R([one(coefficient_ring(R))], [e[p[i]]]) for i in 1:length(p))
end

function terms(a::FreeAssAlgElem, o::WordOrdering)
  p = permutation_of_terms(a, o)
  e = _collected_exponent_words(a)
  c = _collected_exponent_words(a)
  R = parent(a)
  return (R(c[i], [e[p[i]]]) for i in 1:length(p))
end
