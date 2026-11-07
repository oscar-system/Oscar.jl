export comprehensive_groebner_system

@doc raw"""
    comprehensive_groebner_system(I::MPolyIdeal;
      full::AbsAffineAlgebraicSet, hole::AbsAffineAlgebraicSet,
      ordering::MonomialOrdering = default_ordering(base_ring(I)))

If `ordering` is global, return a comprehensive Gröbner system of $I$ with respect to `ordering`.

The output is a vector of triples, with each triple containing a full variety $V(E_i)$,
a hole variety $V(N_i)$, and a vector of polynomials $G_i$.
For every point $a$ of each segment $V(E_i)\setminus V(N_i)$,
substituting $a$ for the parameters in $G_i$ gives a Gröbner basis
of the same substitution in $I$.

!!! note
    The base ring of $I$ must be a multivariate polynomial ring,
    whose coefficient ring is again a multivariate polynomial ring over `QQ`
    or over a prime field created by `GF(p)`.
    The variables of the coefficient ring are the parameters.

!!! note
    The keywords `full` and `hole` are affine algebraic sets in the parameter
    space, defined by ideals of the coefficient ring, and the segments must cover
    $\mathrm{full} \setminus \mathrm{hole}$.
    By default, `full` is the whole parameter space and `hole` is empty,
    so that the segments cover the whole parameter space.

Internally, the function calls `cgsdr` from the Singular library `grobcov.lib`,
which implements the Kapur-Sun-Wang algorithm of [KSW10](@cite).
Additionally, `cgsdr` ensures that the segments are disjoint
and that the Gröbner bases are reduced.

# Examples
```jldoctest
julia> A, (a, b) = polynomial_ring(QQ, [:a, :b]);

julia> R, (x,) = polynomial_ring(A, [:x]);

julia> comprehensive_groebner_system(ideal(R, [a*x - b]))
3-element Vector{Tuple{AffineAlgebraicSet{QQField, MPolyQuoRing{QQMPolyRingElem}}, AffineAlgebraicSet{QQField, MPolyQuoRing{QQMPolyRingElem}}, Vector{AbstractAlgebra.Generic.MPoly{QQMPolyRingElem}}}}:
 (V(0), V(a), [a*x - b])
 (V(a), V(b, a), [1])
 (V(b, a), V(1), [])
```
"""
function comprehensive_groebner_system(I::MPolyIdeal; full::Union{AbsAffineAlgebraicSet,Nothing}=nothing,
    hole::Union{AbsAffineAlgebraicSet,Nothing}=nothing,
    ordering::MonomialOrdering=default_ordering(base_ring(I)))
    R = base_ring(I)
    A = coefficient_ring(R)
    @req A isa MPolyRing "the coefficient ring must be a multivariate polynomial ring"
    K = coefficient_ring(A)
    @req K isa QQField || (K isa FqField && degree(K) == 1) "the coefficients must lie in QQ or in a prime field GF(p)"
    @req base_ring(ordering) === R "the ordering must be defined on the base ring of the ideal"
    @req is_global(ordering) "the ordering must be global"
    full = something(full, algebraic_set(zero(A)))
    hole = something(hole, algebraic_set(one(A)))
    @req ambient_coordinate_ring(full) === A && ambient_coordinate_ring(hole) === A "full and hole must be defined by ideals of the coefficient ring"
    # grobcov.lib needs the parameters in the coefficient field, so we compute in Frac(K[U])[X]
    FX, phi = change_base_ring(fraction_field(A), R)
    o = induce(gens(FX), ordering)
    branches = Singular.LibGrobcov.cgsdr(singular_generators(phi(I), o),
        "null", singular_generators(ideal(FX, gens(fat_ideal(full))), o),
        "nonnull", singular_generators(ideal(FX, gens(fat_ideal(hole))), o))
    # cgsdr returns polynomials over Frac(K[U]) whose coefficients lie in K[U]
    to_A(c) = divexact(numerator(c), denominator(c))
    function to_set(J)
        v = elem_type(A)[to_A(constant_coefficient(h)) for h in gens(ideal(FX, J)) if !iszero(h)]
        return algebraic_set(ideal(A, isempty(v) ? [zero(A)] : v))
    end
    to_R(J) = elem_type(R)[map_coefficients(to_A, g; parent=R) for g in gens(ideal(FX, J)) if !iszero(g)]
    return [(to_set(E), to_set(N), to_R(G)) for (E, N, G) in branches]
end
