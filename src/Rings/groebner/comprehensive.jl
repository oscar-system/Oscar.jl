@doc raw"""
    comprehensive_groebner_system(I::MPolyIdeal{<:MPolyRingElem{T}};
      null::MPolyIdeal{T}, nonnull::MPolyIdeal{T},
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      algorithm::Symbol = :KSW) where {T <: MPolyRingElem}

If `ordering` is global, return a comprehensive Gröbner system of $I$ with respect to `ordering`.

The output is a vector of triples $(E_i, N_i, G_i)$,
where $E_i$ and $N_i$ are ideals of the coefficient ring
and $G_i$ is an ideal of the base ring of $I$.
For every point $a$ of each segment $V(E_i)\setminus V(N_i)$,
substituting $a$ for the parameters in the generators of $G_i$ gives a Gröbner basis
of the same substitution in $I$.

!!! note
    The base ring of $I$ must be a multivariate polynomial ring,
    whose coefficient ring is again a multivariate polynomial ring.
    The variables of the coefficient ring are the parameters.

!!! note
    The keywords `null` and `nonnull` are ideals of the coefficient ring,
    and the segments cover $V(\mathrm{null}) \setminus V(\mathrm{nonnull})$.
    By default, `null` is the zero ideal and `nonnull` is the unit ideal,
    so that the segments cover the whole parameter space.

The keyword `algorithm` can be set to
- `:KSW` (implementation of the Kapur-Sun-Wang algorithm [KSW10](@cite) in *Singular*).
  This implementation requires the coefficient ring of the parameters to be
  `QQ` or a prime field created by `GF(p)`.

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
function comprehensive_groebner_system(I::MPolyIdeal{<:MPolyRingElem{T}};
    null::MPolyIdeal{T}=ideal(coefficient_ring(base_ring(I)), 0),
    nonnull::MPolyIdeal{T}=ideal(coefficient_ring(base_ring(I)), 1),
    ordering::MonomialOrdering=default_ordering(base_ring(I)),
    algorithm::Symbol=:KSW) where {T<:MPolyRingElem}
  @req algorithm in (:KSW,) "Unsupported algorithm"
  R = base_ring(I)
  A = coefficient_ring(R)
  K = coefficient_ring(A)
  @req K isa QQField || (K isa FqField && degree(K) == 1) "the coefficients must lie in QQ or in a prime field GF(p)"
  @req base_ring(ordering) === R "the ordering must be defined on the base ring of the ideal"
  @req is_global(ordering) "the ordering must be global"
  @req base_ring(null) === A && base_ring(nonnull) === A "null and nonnull must be ideals of the coefficient ring"
  # grobcov.lib needs the parameters in the coefficient field, so we compute in Frac(K[U])[X]
  FX, phi = change_base_ring(fraction_field(A), R)
  o = induce(gens(FX), ordering)
  # cgsdr sets the global Singular options redSB and returnSB without restoring them
  redSB = Singular.libSingular.set_option("OPT_REDSB", false)
  returnSB = Singular.libSingular.set_option("OPT_RETURN_SB", false)
  branches = Singular.LibGrobcov.cgsdr(singular_generators(phi(I), o),
      "null", singular_generators(ideal(FX, gens(null)), o),
      "nonnull", singular_generators(ideal(FX, gens(nonnull)), o))
  Singular.libSingular.set_option("OPT_REDSB", redSB)
  Singular.libSingular.set_option("OPT_RETURN_SB", returnSB)
  # cgsdr returns polynomials over Frac(K[U]) whose coefficients lie in K[U]
  cgs = Tuple{MPolyIdeal{elem_type(A)}, MPolyIdeal{elem_type(A)}, MPolyIdeal{elem_type(R)}}[]
  for branch in branches
    E, N, G = ([map_coefficients(c -> divexact(numerator(c), denominator(c)), g; parent=R)
                for g in gens(ideal(FX, J))] for J in branch)
    push!(cgs, (ideal(A, constant_coefficient.(E)), ideal(A, constant_coefficient.(N)), ideal(R, G)))
  end
  return cgs
end
