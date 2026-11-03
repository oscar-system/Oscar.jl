# ---------------------------------------------------------------------------
# User-facing wrappers and verification helpers corresponding to `primary_decomposition_native`,
# `minAss`, `testPrimary`, `idealsEqual` and friends of Singular's `primdec.lib`.
# ---------------------------------------------------------------------------

@doc raw"""
    primary_decomposition_native(I::MPolyIdeal; algorithm::Symbol = :gtz) -> Vector{Tuple{MPolyIdeal, MPolyIdeal}}

Compute a primary decomposition of `I` as a vector of pairs `(Q, P)` of a
primary component and its radical, with the algorithm `:gtz` (Gianni, Trager,
Zacharias) or `:sy` (Shimoyama, Yokoyama).
"""
function primary_decomposition_native(I::MPolyIdeal; algorithm::Symbol=:gtz)
  algorithm == :gtz && return primary_decomposition_gtz(I)
  algorithm == :sy && return primary_decomposition_sy(I)
  throw(ArgumentError("unknown algorithm $algorithm, expected :gtz or :sy"))
end

@doc raw"""
    minimal_associated(I::MPolyIdeal; algorithm::Symbol = :gtz) -> Vector{MPolyIdeal}

Compute the minimal associated primes of `I` with the algorithm `:gtz`.
"""
function minimal_associated(I::MPolyIdeal; algorithm::Symbol=:gtz)
  algorithm == :gtz && return minimal_associated_gtz(I)
  throw(ArgumentError("unknown algorithm $algorithm, expected :gtz"))
end

@doc raw"""
    test_primary(L::Vector{Tuple{MPolyIdeal, MPolyIdeal}}, I::MPolyIdeal) -> Bool

Check that the pairs `(Q, P)` in `L` consist of primary ideals `Q` with prime
radical `P` whose intersection is `I`, see `testPrimary` of Singular.
"""
function test_primary(L::Vector{<:Tuple{<:MPolyIdeal,<:MPolyIdeal}}, I::MPolyIdeal)
  isempty(L) && return is_one(I)
  inter = length(L) == 1 ? L[1][1] : intersect([Q for (Q, _) in L]...)
  _pd_ideals_equal(inter, I) || return false
  for (Q, P) in L
    _pd_ideals_equal(radical(Q), P) || return false
    is_prime(P) || return false
    is_primary(Q) || return false
  end
  return true
end
