########################
# 1: Special attributes of toric varieties
########################

@doc Markdown.doc"""
    cohomology_ring(v::AbstractNormalToricVariety)

Return the cohomology ring of the simplicial and complete toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> ngens(cohomology_ring(p2))
3
```
"""
function cohomology_ring(v::AbstractNormalToricVariety)
    if ((issimplicial(v) == false) || (iscomplete(v) == false))
        throw(ArgumentError("The cohomology ring is (currently) only supported for simplicial and complete toric varieties."))
    end
    R, _ = PolynomialRing(coefficient_ring(v), coordinate_names(v), cached = false)
    weights = [1 for i in 1:ngens(R)]
    R = grade(R, weights)[1]
    linear_relations = ideal_of_linear_relations(R, v)
    stanley_reisner = stanley_reisner_ideal(R,v)
    return quo(R, linear_relations + stanley_reisner)[1]
end
export cohomology_ring


@doc Markdown.doc"""
    volume_form(v::NormalToricVariety)

Construct the volume form of the normal
toric toric variety `v`.

# Examples
```jldoctest
julia> polynomial(volume_form(projective_space(NormalToricVariety, 2)))
x3^2

julia> polynomial(volume_form(del_pezzo(3)))
-e2^2

julia> polynomial(volume_form(hirzebruch_surface(5)))
1//5*x2^2
```
"""
@attr CohomologyClass function volume_form(v::NormalToricVariety)
    mc = ray_indices(maximal_cones(v))
    exponents = [fmpz(mc[1,i]) for i in 1:length(mc[1,:])]
    indets = gens(cohomology_ring(v))
    poly = fmpz(1) * prod([indets[k]^exponents[k] for k in 1:length(exponents)])
    if (iszero(poly) || degree(poly)[1] != dim(v))
        throw(ArgumentError("The volume class does not exist."))
    end
    return CohomologyClass(v, poly)
end
export volume_form


@attr function _intersection_form_via_exponents(v::NormalToricVariety)
    # extract the cohomology classes corresponding to the torus-invariant prime divisors
    generators = [CohomologyClass(d) for d in torusinvariant_prime_divisors(v)]
    
    # find combinations of those classes that we have to integrate
    S, _ = PolynomialRing(QQ, ["g$(i)" for i in 1:length(generators)], cached=false)
    S = grade(S, [1 for i in 1:ngens(S)])[1]
    hc = homogeneous_component(S, [dim(v)])
    monoms = [hc[2](x) for x in gens(hc[1])]
    combinations = reduce(vcat,[[[fmpz(l) for l in k] for k in exponent_vectors(m)] for m in monoms])
    
    # perform the integrals
    intersection_dict = Dict{fmpz_mat, fmpq}()
    for expos in combinations
        cc = prod([generators[k]^expos[k] for k in 1:length(generators)])
        intersection_dict[matrix(ZZ,[expos])] = integrate(cc)
    end
    
    # return the dictionary
    return intersection_dict
end


@doc Markdown.doc"""
    intersection_form(v::NormalToricVariety)

Computes the intersection numbers among the cohomology classes
associated to the torusinvariant prime divisors of the normal toric toric variety `v`.

# Examples
```jldoctest
julia> F3 = hirzebruch_surface(3)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> length(intersection_form(F3))
10
```
"""
function intersection_form(v::NormalToricVariety)
    intersection_dict = _intersection_form_via_exponents(v)
    monoms_of_prime_divisors = gens(cox_ring(v))    
    intersection_dict_for_user = Dict{MPolyElem, fmpq}()
    for (expos,v) in intersection_dict
        monom = prod([monoms_of_prime_divisors[k]^expos[k] for k in 1:length(monoms_of_prime_divisors)])
        intersection_dict_for_user[monom] = v
    end
    return intersection_dict_for_user
end
export intersection_form
