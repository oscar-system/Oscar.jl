###########################################################################
# 1. Defining attributes of rational equivalence class of algebraic cycles
###########################################################################


@doc raw"""
    toric_variety(ac::RationalEquivalenceClass)

Return the normal toric variety of a rational
equivalence class of algebraic cycles.

# Examples
```jldoctest
julia> dP2 = del_pezzo_surface(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = toric_divisor(dP2, [1, 2, 3, 4, 5])
Torus-invariant, non-prime divisor on a normal toric variety

julia> ac = rational_equivalence_class(d)
Rational equivalence class on a normal toric variety represented by 6V(x3)+V(e1)+7V(e2)

julia> toric_variety(ac)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
toric_variety(ac::RationalEquivalenceClass) = ac.v


@doc raw"""
    polynomial(ac::RationalEquivalenceClass)

On a simplicial and complete toric variety, the Chow ring
is isomorphic to a certain quotient of the Cox ring. This
function returns the ring element corresponding to a given
rational equivalence class of algebraic cycles.

# Examples
```jldoctest
julia> dP2 = del_pezzo_surface(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = toric_divisor(dP2, [1, 2, 3, 4, 5])
Torus-invariant, non-prime divisor on a normal toric variety

julia> ac = rational_equivalence_class(d)
Rational equivalence class on a normal toric variety represented by 6V(x3)+V(e1)+7V(e2)

julia> polynomial(ac)
6*x3 + e1 + 7*e2
```
"""
polynomial(ac::RationalEquivalenceClass) = ac.p


@doc raw"""
    polynomial(ring::MPolyQuoRing, ac::RationalEquivalenceClass)

On a simplicial and complete toric variety, the Chow ring
is isomorphic to a certain quotient of the Cox ring. This
function returns the ring element corresponding to a given
rational equivalence class of algebraic cycles. The first
argument of this function allows to obtain this ring element
in a different ring. This allows to change the coefficient
ring if desired.

# Examples
```jldoctest
julia> dP2 = del_pezzo_surface(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = toric_divisor(dP2, [1, 2, 3, 4, 5])
Torus-invariant, non-prime divisor on a normal toric variety

julia> ac = rational_equivalence_class(d)
Rational equivalence class on a normal toric variety represented by 6V(x3)+V(e1)+7V(e2)

julia> R, _ = polynomial_ring(QQ, 5)
(Multivariate Polynomial Ring in x1, x2, x3, x4, x5 over Rational Field, QQMPolyRingElem[x1, x2, x3, x4, x5])

julia> (x1, x2, x3, x4, x5) = gens(R)
5-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3
 x4
 x5

julia> sr_and_linear_relation_ideal = ideal([x1*x3, x1*x5, x2*x4, x2*x5, x3*x4, x1 + x2 - x5, x2 + x3 - x4 - x5])
ideal(x1*x3, x1*x5, x2*x4, x2*x5, x3*x4, x1 + x2 - x5, x2 + x3 - x4 - x5)

julia> R_quo = quo(R, sr_and_linear_relation_ideal)[1]
Quotient of Multivariate Polynomial Ring in x1, x2, x3, x4, x5 over Rational Field by ideal(x1*x3, x1*x5, x2*x4, x2*x5, x3*x4, x1 + x2 - x5, x2 + x3 - x4 - x5)

julia> polynomial(R_quo, ac)
6*x3 + x4 + 7*x5
```
"""
function polynomial(ring::MPolyQuoRing, ac::RationalEquivalenceClass)
    p = polynomial(ac)
    if iszero(p)
        return zero(ring)
    end
    coeffs = [k for k in AbstractAlgebra.coefficients(p.f)]
    expos = matrix(ZZ, [k for k in AbstractAlgebra.exponent_vectors(p.f)])
    indets = gens(ring)
    monoms = [prod(indets[j]^expos[k, j] for j in 1:ncols(expos)) for k in 1:nrows(expos)]
    return sum(coeffs[k]*monoms[k] for k in 1:length(monoms))
end


###########################################################################
# 2. Representing rational equivalence classes of algebraic cycles
###########################################################################


@doc raw"""
    representative(ac::RationalEquivalenceClass)

Return a polynomial in the Cox ring mapping to `polynomial(ac)`.

# Examples
```jldoctest
julia> dP2 = del_pezzo_surface(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = toric_divisor(dP2, [1, 2, 3, 4, 5])
Torus-invariant, non-prime divisor on a normal toric variety

julia> ac = rational_equivalence_class(d)
Rational equivalence class on a normal toric variety represented by 6V(x3)+V(e1)+7V(e2)

julia> ac*ac
Rational equivalence class on a normal toric variety represented by 34V(x2,x3)

julia> representative(ac*ac)
34*x2*x3
```
"""
@attr MPolyDecRingElem{QQFieldElem, QQMPolyRingElem} function representative(ac::RationalEquivalenceClass)
    if is_trivial(ac)
        return zero(cox_ring(toric_variety(ac)))
    end
    coeffs = coefficients(ac)
    mapped_monomials = [map_gens_of_chow_ring_to_cox_ring(toric_variety(ac))[m] for m in monomials(polynomial(ac).f)]
    return sum([coeffs[i]*mapped_monomials[i] for i in 1:length(mapped_monomials)])
end


@doc raw"""
    coefficients(ac::RationalEquivalenceClass)

Return the coefficients of `polynomial(ac)`.

# Examples
```jldoctest
julia> dP2 = del_pezzo_surface(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = toric_divisor(dP2, [1, 2, 3, 4, 5])
Torus-invariant, non-prime divisor on a normal toric variety

julia> ac = rational_equivalence_class(d)
Rational equivalence class on a normal toric variety represented by 6V(x3)+V(e1)+7V(e2)

julia> coefficients(ac*ac)
1-element Vector{QQFieldElem}:
 -34
```
"""
@attr Vector{QQFieldElem} function coefficients(ac::RationalEquivalenceClass)
    if is_trivial(ac)
        return QQFieldElem[]
    end
    return [coefficient_ring(toric_variety(ac))(k) for k in AbstractAlgebra.coefficients(polynomial(ac).f)]
end


@doc raw"""
    components(ac::RationalEquivalenceClass)

Turn each monomial of `representative(ac)` into a
closed subvariety and return the list formed from these
subvarieties. Note that each of these subvarieties is
irreducible and their formal linear sum, with the coefficients
computed by the method `coefficients(ac::RationalEquivalenceClass)`,
defines an algebraic cycle, whose rational equivalence
class is identical to the one given to this method.

# Examples
```jldoctest
julia> dP2 = del_pezzo_surface(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = toric_divisor(dP2, [1, 2, 3, 4, 5])
Torus-invariant, non-prime divisor on a normal toric variety

julia> ac = rational_equivalence_class(d)
Rational equivalence class on a normal toric variety represented by 6V(x3)+V(e1)+7V(e2)

julia> length(components(ac*ac))
1
```
"""
@attr Vector{ClosedSubvarietyOfToricVariety} function components(ac::RationalEquivalenceClass)
    if is_trivial(ac)
        return ClosedSubvarietyOfToricVariety[]
    end
    variety = toric_variety(ac)
    gs = gens(cox_ring(toric_variety(ac)))
    mons = [m for m in monomials(representative(ac))]
    expos = [[e for e in AbstractAlgebra.exponent_vectors(m)][1] for m in mons]
    return [closed_subvariety_of_toric_variety(variety, [gs[k] for k in findall(!iszero, exps)]) for exps in expos]
end


###########################################################################
# 3. Other attributes of rational equivalence class of algebraic cycles
###########################################################################


@doc raw"""
    cohomology_class(ac::RationalEquivalenceClass)

Return the cohomology class of a rational
equilvalence class of algebraic cycles.

# Examples
```jldoctest
julia> dP2 = del_pezzo_surface(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = toric_divisor(dP2, [1, 2, 3, 4, 5])
Torus-invariant, non-prime divisor on a normal toric variety

julia> ac = rational_equivalence_class(d)
Rational equivalence class on a normal toric variety represented by 6V(x3)+V(e1)+7V(e2)

julia> cohomology_class(ac)
Cohomology class on a normal toric variety given by 6*x3 + e1 + 7*e2
```
"""
@attr CohomologyClass cohomology_class(ac::RationalEquivalenceClass) = CohomologyClass(toric_variety(ac), polynomial(cohomology_ring(toric_variety(ac)),ac))
