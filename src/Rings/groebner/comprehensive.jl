@doc raw"""
    comprehensive_groebner_system(I::MPolyIdeal;
      null::MPolyIdeal, nonnull::MPolyIdeal,
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      algorithm::Symbol = :KSW)

If `ordering` is global, return a comprehensive Gröbner system of $I$ with respect to `ordering`.

The output is a vector of triples $(E_i, N_i, G_i)$,
where $E_i$ and $N_i$ are ideals of the coefficient ring
and $G_i$ is an ideal of the base ring of $I$.
For every point $a$ of each segment $V(E_i)\setminus V(N_i)$,
substituting $a$ for the parameters in the generators of $G_i$ gives a Gröbner basis
of the same substitution in $I$.

!!! note
    The base ring of $I$ must be a multivariate polynomial ring,
    whose coefficient ring is again a multivariate polynomial ring over `QQ`
    or over a prime field created by `GF(p)`.
    The variables of the coefficient ring are the parameters.

!!! note
    The keywords `null` and `nonnull` are ideals of the coefficient ring,
    and the segments cover $V(\mathrm{null}) \setminus V(\mathrm{nonnull})$.
    By default, `null` is the zero ideal and `nonnull` is the unit ideal,
    so that the segments cover the whole parameter space.

The keyword `algorithm` can be set to
- `:KSW` (implementation of the Kapur-Sun-Wang algorithm [KSW10](@cite) in *Singular*).

# Examples
```jldoctest
julia> A, (a, b) = polynomial_ring(QQ, [:a, :b]);

julia> R, (x,) = polynomial_ring(A, [:x]);

julia> comprehensive_groebner_system(ideal(R, [a*x - b]))
3-element Vector{Tuple{MPolyIdeal{QQMPolyRingElem}, MPolyIdeal{QQMPolyRingElem}, MPolyIdeal{AbstractAlgebra.Generic.MPoly{QQMPolyRingElem}}}}:
 (Ideal (0), Ideal (a), Ideal (a*x - b))
 (Ideal (a), Ideal (b, a), Ideal (1))
 (Ideal (b, a), Ideal (1), Ideal (0))
```
"""
function comprehensive_groebner_system(I::MPolyIdeal;
    null::Union{MPolyIdeal,Nothing}=nothing, nonnull::Union{MPolyIdeal,Nothing}=nothing,
    ordering::MonomialOrdering=default_ordering(base_ring(I)),
    algorithm::Symbol=:KSW)
  @req algorithm in (:KSW,) "Unsupported algorithm"
  R = base_ring(I)
  A = coefficient_ring(R)
  @req A isa MPolyRing "the coefficient ring must be a multivariate polynomial ring"
  K = coefficient_ring(A)
  @req K isa QQField || (K isa FqField && degree(K) == 1) "the coefficients must lie in QQ or in a prime field GF(p)"
  @req base_ring(ordering) === R "the ordering must be defined on the base ring of the ideal"
  @req is_global(ordering) "the ordering must be global"
  null = something(null, ideal(A, [zero(A)]))
  nonnull = something(nonnull, ideal(A, [one(A)]))
  @req base_ring(null) === A && base_ring(nonnull) === A "null and nonnull must be ideals of the coefficient ring"
  # grobcov.lib needs the parameters in the coefficient field, so we compute in Frac(K[U])[X]
  FX, phi = change_base_ring(fraction_field(A), R)
  o = induce(gens(FX), ordering)
  branches = Singular.LibGrobcov.cgsdr(singular_generators(phi(I), o),
      "null", singular_generators(ideal(FX, gens(null)), o),
      "nonnull", singular_generators(ideal(FX, gens(nonnull)), o))
  # cgsdr returns polynomials over Frac(K[U]) whose coefficients lie in K[U]
  to_A(c) = divexact(numerator(c), denominator(c))
  function to_ideal(S, f, J)
    v = elem_type(S)[f(g) for g in gens(ideal(FX, J)) if !iszero(g)]
    return ideal(S, isempty(v) ? [zero(S)] : v)
  end
  to_A_ideal(J) = to_ideal(A, g -> to_A(constant_coefficient(g)), J)
  to_R_ideal(J) = to_ideal(R, g -> map_coefficients(to_A, g; parent=R), J)
  return [(to_A_ideal(E), to_A_ideal(N), to_R_ideal(G)) for (E, N, G) in branches]
end
