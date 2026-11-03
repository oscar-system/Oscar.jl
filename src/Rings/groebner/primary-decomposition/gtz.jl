# ---------------------------------------------------------------------------
# The algorithm of Gianni, Trager and Zacharias for primary decomposition and
# minimal associated primes, see `primdecGTZ` and `minAssGTZ` of Singular's
# `primdec.lib`. For a maximal independent set `U` of the ideal `I` and the
# lcm `h` of the leading coefficients of a Groebner basis with respect to a block
# ordering with the variables `U` last,
#   I = (I : h^oo) ∩ (I + <h^k>),
# where the first ideal is the contraction of a zero-dimensional ideal of
# `K(U)[x]`, which is decomposed by `_pd_zero_dimensional_decomposition`, and
# the second ideal is treated recursively.
# ---------------------------------------------------------------------------

# One step of the algorithm. Returns the pairs `(Q, P)` of primary components
# and associated primes of the part `I : h^oo`, as well as the list of ideals
# that are left to be decomposed. In the mode `:primary` this is `I + <f^k>`, in
# the other mode (for associated primes) it is `I + <g>` for the irreducible
# factors `g` of `f`.
function _pd_gtz_step(I::MPolyIdeal; mode::Symbol=:primary)
  unit = ideal(base_ring(I), [one(base_ring(I))])
  R = base_ring(I)
  n = ngens(R)
  T = typeof(I)
  d = dim(I)
  if d == 0
    return _pd_zero_dimensional_decomposition(I), typeof(I)[]
  end
  if d == n
    # the zero ideal
    return Tuple{T,T}[(I, I)], T[]
  end
  U = _pd_independent_set(I)
  E = _pd_extension(R, U)
  J = _pd_extend_ideal(E, I)
  comps = _pd_zero_dimensional_decomposition(J)
  res = Tuple{T,T}[]
  for (Q, P) in comps
    push!(res, (_pd_contract(E, gens(Q)), _pd_contract(E, gens(P))))
  end
  h = _pd_lead_lcm(I, U)
  _, facs, f, k = _pd_split_by(I, h)
  if k == 0
    rest = T[]
  elseif mode == :primary
    rest = [I + ideal(R, [f^k])]
  else
    # every prime containing `I` and `f` contains some irreducible factor of `f`
    rest = [I + ideal(R, [g]) for g in facs]
  end
  return res, rest
end

@doc raw"""
    primary_decomposition_gtz(I::MPolyIdeal) -> Vector{Tuple{MPolyIdeal, MPolyIdeal}}

Compute a primary decomposition of `I` with the algorithm of Gianni, Trager and
Zacharias. The result is a vector of pairs `(Q, P)` where `Q` is a primary
component and `P` its radical (the associated prime). The decomposition is
irredundant: no component can be removed without changing the intersection.
The base field of the polynomial ring has to be `QQ` or a finite field.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [x*y, x*z]);

julia> length(Oscar.primary_decomposition_gtz(I))
2
```
"""
function primary_decomposition_gtz(I::MPolyIdeal)
  R = base_ring(I)
  T = typeof(I)
  I = _pd_clean(I)
  is_one(I) && return Tuple{T,T}[]
  res = Tuple{T,T}[]
  todo = [I]
  while !isempty(todo)
    cur = pop!(todo)
    comps, rest = _pd_gtz_step(cur; mode=:primary)
    append!(res, comps)
    append!(todo, filter(!is_one, rest))
  end
  return _pd_irredundant(res)
end

# Remove duplicates and redundant primary components, and merge components
# belonging to the same prime.
function _pd_irredundant(L::Vector{Tuple{T,T}}) where {T<:MPolyIdeal}
  isempty(L) && return L
  # merge components with the same prime
  merged = Tuple{T,T}[]
  for (Q, P) in L
    idx = findfirst(((_, P2),) -> _pd_ideals_equal(P, P2), merged)
    if idx === nothing
      push!(merged, (Q, P))
    else
      merged[idx] = (intersect(merged[idx][1], Q), P)
    end
  end
  # drop components whose removal does not change the intersection
  keep = trues(length(merged))
  for i in eachindex(merged)
    others = [merged[j][1] for j in eachindex(merged) if j != i && keep[j]]
    isempty(others) && continue
    inter = length(others) == 1 ? others[1] : intersect(others...)
    if is_subset(inter, merged[i][1])
      keep[i] = false
    end
  end
  return merged[keep]
end

@doc raw"""
    minimal_associated_gtz(I::MPolyIdeal) -> Vector{MPolyIdeal}

Compute the minimal associated primes of `I` with the algorithm of Gianni,
Trager and Zacharias.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> length(Oscar.minimal_associated_gtz(ideal(R, [x*y, x*z])))
2
```
"""
function minimal_associated_gtz(I::MPolyIdeal)
  T = typeof(I)
  I = _pd_clean(I)
  is_one(I) && return T[]
  primes = T[]
  todo = [I]
  while !isempty(todo)
    cur = pop!(todo)
    # only the part of `cur` that is not already covered by known primes is relevant
    comps, rest = _pd_gtz_step(cur; mode=:radical)
    for (_, P) in comps
      push!(primes, P)
    end
    append!(todo, filter(!is_one, rest))
  end
  return _pd_minimal_elements(primes)
end
