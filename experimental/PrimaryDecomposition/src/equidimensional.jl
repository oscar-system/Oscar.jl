# ---------------------------------------------------------------------------
# Equidimensional decompositions, see `equidim`, `equidimMax`, `equiRadical`
# and `prepareAss` of Singular's `primdec.lib`. Each GTZ step splits off an
# equidimensional part `I : h^oo` whose associated primes all have the
# dimension of the independent set `U`.
# ---------------------------------------------------------------------------

# Equidimensional pieces `(dimension, ideal)` with `I` the intersection of the
# ideals; the pieces are in order of production, not grouped.
function _pd_equidim_pieces(I::MPolyIdeal; radical::Bool=false)
  R = base_ring(I)
  n = ngens(R)
  unit = ideal(R, [one(R)])
  pieces = Tuple{Int,typeof(I)}[]
  todo = [_pd_clean(I)]
  while !isempty(todo)
    cur = pop!(todo)
    is_one(cur) && continue
    d = dim(cur)
    if d == 0 || d == n
      push!(pieces, (d, radical && d == 0 ? _pd_zero_radical(cur) : cur))
      continue
    end
    U = _pd_independent_set(cur)
    h = _pd_lead_lcm(cur, U)
    if radical
      E = _pd_extension(R, U)
      S = _pd_contract(E, gens(_pd_zero_radical(_pd_extend_ideal(E, cur))))
      _, facs, _, _ = _pd_split_by(cur, h)
      rest = [cur + ideal(R, [g]) for g in facs]
    else
      S, _, f, k = _pd_split_by(cur, h)
      rest = k == 0 ? typeof(cur)[] : [cur + ideal(R, [f^k])]
    end
    push!(pieces, (d, S))
    append!(todo, rest)
  end
  return pieces
end

# The pieces grouped by dimension (intersected) in ascending order of dimension.
function _pd_group_pieces(pieces)
  isempty(pieces) && return Tuple{Int,Any}[]
  res = Tuple{Int,typeof(pieces[1][2])}[]
  for d in sort!(unique!([p[1] for p in pieces]))
    ids = [p[2] for p in pieces if p[1] == d]
    push!(res, (d, length(ids) == 1 ? ids[1] : intersect(ids...)))
  end
  return res
end

# Radicals of the equidimensional components, in ascending dimension.
function _pd_equidim_radicals(I::MPolyIdeal)
  return [P for (_, P) in _pd_group_pieces(_pd_equidim_pieces(I; radical=true))]
end

@doc raw"""
    equidim(I::MPolyIdeal) -> Vector{MPolyIdeal}

Compute equidimensional ideals `a[1], ..., a[s]` of increasing dimension whose
intersection is `I`, where `a[s]` is the equidimensional locus of `I`, i.e. the
intersection of the primary components of `I` of the maximal dimension. For
the unit ideal the list consisting of the unit ideal is returned.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> length(Oscar.equidim(ideal(R, [x*y, x*z])))
2
```
"""
function equidim(I::MPolyIdeal)
  I = _pd_clean(I)
  is_one(I) && return [I]
  return [S for (_, S) in _pd_group_pieces(_pd_equidim_pieces(I))]
end

@doc raw"""
    equidim_max(I::MPolyIdeal) -> MPolyIdeal

Compute the equidimensional locus of maximal dimension of `I`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> Oscar.equidim_max(ideal(R, [x*y, x*z])) == ideal(R, [x])
true
```
"""
function equidim_max(I::MPolyIdeal)
  return last(equidim(I))
end

@doc raw"""
    equi_radical(I::MPolyIdeal) -> MPolyIdeal

Compute the intersection of the associated primes of `I` of maximal dimension.
"""
function equi_radical(I::MPolyIdeal)
  I = _pd_clean(I)
  is_one(I) && return I
  return last(_pd_equidim_radicals(I))
end

@doc raw"""
    prepare_ass(I::MPolyIdeal) -> Vector{MPolyIdeal}

Compute the radicals of the equidimensional components of `I`, in ascending
order of dimension.
"""
function prepare_ass(I::MPolyIdeal)
  I = _pd_clean(I)
  is_one(I) && return [I]
  return _pd_equidim_radicals(I)
end
