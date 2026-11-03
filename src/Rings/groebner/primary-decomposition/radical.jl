# ---------------------------------------------------------------------------
# Radical computation following Krick and Logar (with the modifications of
# Laplagne), see `radical` of Singular's `primdec.lib`: with `U` and `h` as in the
# GTZ algorithm, the radical of `I` is the intersection of the radical of
# `I : h^oo` (the contraction of the radical of a zero-dimensional ideal over
# `K(U)`) and the radical of `I + <h>`.
# ---------------------------------------------------------------------------

# The part `rad(I : h^oo)` together with the ideals `I + <g>` that remain, for the
# irreducible factors `g` of a minimal divisor of `h`.
function _pd_radical_step(I::MPolyIdeal)
  R = base_ring(I)
  n = ngens(R)
  d = dim(I)
  d == 0 && return _pd_zero_radical(I), typeof(I)[]
  d == n && return I, typeof(I)[]
  U = _pd_independent_set(I)
  E = _pd_extension(R, U)
  J = _pd_extend_ideal(E, I)
  Jr = _pd_zero_radical(J)
  S = _pd_contract(E, gens(Jr))
  h = _pd_lead_lcm(I, U)
  _, facs, _, _ = _pd_split_by(I, h)
  return S, [I + ideal(R, [g]) for g in facs]
end

@doc raw"""
    radical_kl(I::MPolyIdeal) -> MPolyIdeal

Compute the radical of `I` with the algorithm of Krick and Logar.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> Oscar.radical_kl(ideal(R, [x^2*y, x*y^2])) == ideal(R, [x*y])
true
```
"""
function radical_kl(I::MPolyIdeal)
  R = base_ring(I)
  I = _pd_clean(I)
  is_one(I) && return I
  parts = typeof(I)[]
  todo = [I]
  while !isempty(todo)
    cur = pop!(todo)
    S, rest = _pd_radical_step(cur)
    push!(parts, S)
    append!(todo, filter(!is_one, rest))
  end
  return length(parts) == 1 ? parts[1] : intersect(parts...)
end

@doc raw"""
    radical_equidim(I::MPolyIdeal) -> MPolyIdeal

Compute the radical of `I` as the intersection of the radicals of the
equidimensional components, in the spirit of `radicalEHV` in Singular's `primdec.lib`. The
equidimensional components are computed with the GTZ step instead of
Ext-modules.
"""
function radical_equidim(I::MPolyIdeal)
  parts = _pd_equidim_radicals(I)
  isempty(parts) && return _pd_clean(I)
  return length(parts) == 1 ? parts[1] : intersect(parts...)
end
