########################
# 1: Special attributes of toric varieties
########################

@doc raw"""
    cohomology_ring(v::AbstractNormalToricVariety)

Return the cohomology ring of the simplicial and complete toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> ngens(cohomology_ring(p2))
3
```
"""
@attr MPolyQuoRing function cohomology_ring(v::AbstractNormalToricVariety)
    @req is_simplicial(v) && is_complete(v) "The cohomology ring is only supported for simplicial and complete toric varieties"
    R, _ = polynomial_ring(coefficient_ring(v), coordinate_names(v), cached = false)
    weights = [1 for i in 1:ngens(R)]
    R = grade(R, weights)[1]
    linear_relations = ideal_of_linear_relations(R, v)
    stanley_reisner = stanley_reisner_ideal(R, v)
    return quo(R, linear_relations + stanley_reisner)[1]
end


@doc raw"""
    volume_form(v::NormalToricVariety)

Construct the volume form of the normal
toric toric variety `v`.

# Examples
```jldoctest
julia> polynomial(volume_form(projective_space(NormalToricVariety, 2)))
x3^2

julia> polynomial(volume_form(del_pezzo_surface(NormalToricVariety, 3)))
-e3^2

julia> polynomial(volume_form(hirzebruch_surface(NormalToricVariety, 5)))
1//5*x2^2
```
"""
@attr CohomologyClass function volume_form(v::NormalToricVariety)
    mc = ray_indices(maximal_cones(v))
    exponents = [ZZRingElem(mc[1, i]) for i in 1:length(mc[1, :])]
    indets = gens(cohomology_ring(v))
    poly = prod(indets[k]^exponents[k] for k in 1:length(exponents))
    @req !iszero(poly) && degree(poly)[1] == dim(v) "The volume class does not exist"
    return CohomologyClass(v, poly)
end


@attr function _intersection_form_via_exponents(v::NormalToricVariety)
    # extract the cohomology classes corresponding to the torus-invariant prime divisors
    generators = [cohomology_class(d) for d in torusinvariant_prime_divisors(v)]
    
    # find combinations of those classes that we have to integrate
    S, _ = polynomial_ring(QQ, ["g$(i)" for i in 1:length(generators)], cached=false)
    S = grade(S, [1 for i in 1:ngens(S)])[1]
    hc = homogeneous_component(S, [dim(v)])
    monoms = [hc[2](x) for x in gens(hc[1])]
    combinations = reduce(vcat, [[[ZZRingElem(l) for l in k] for k in AbstractAlgebra.exponent_vectors(m)] for m in monoms])
    
    # perform the integrals
    intersection_dict = Dict{ZZMatrix, QQFieldElem}()
    for expos in combinations
        cc = prod(generators[k]^expos[k] for k in 1:length(generators))
        intersection_dict[matrix(ZZ, [expos])] = integrate(cc)
    end
    
    # return the dictionary
    return intersection_dict
end


@doc raw"""
    intersection_form(v::NormalToricVariety)

Computes the intersection numbers among the cohomology classes
associated to the torusinvariant prime divisors of the normal toric toric variety `v`.

# Examples
```jldoctest
julia> F3 = hirzebruch_surface(NormalToricVariety, 3)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> length(intersection_form(F3))
10
```
"""
function intersection_form(v::NormalToricVariety)
    intersection_dict = _intersection_form_via_exponents(v)
    monoms_of_prime_divisors = gens(cox_ring(v))
    intersection_dict_for_user = Dict{MPolyRingElem, QQFieldElem}()
    for (expos, v) in intersection_dict
        monom = prod(monoms_of_prime_divisors[k]^expos[k] for k in 1:length(monoms_of_prime_divisors))
        intersection_dict_for_user[monom] = v
    end
    return intersection_dict_for_user
end
