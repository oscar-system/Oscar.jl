############################
# Dimensions
############################


@doc raw"""
    dim(v::AbstractNormalToricVariety)

Return the dimension of the normal toric variety `v`.

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> dim(antv)
1
```
"""
@attr Int dim(v::AbstractNormalToricVariety) = pm_object(v).FAN_DIM


@doc raw"""
    dim_of_torusfactor(v::AbstractNormalToricVariety)

Return the dimension of the torus factor of the normal toric variety `v`.

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> dim_of_torusfactor(antv)
1
```
"""
@attr Int function dim_of_torusfactor(v::AbstractNormalToricVariety)
    if has_torusfactor(v) == false
        return 0
    end
    dimension_of_fan = pm_object(v).FAN_DIM::Int
    ambient_dimension = pm_object(v).FAN_AMBIENT_DIM::Int
    return ambient_dimension - dimension_of_fan
end


@doc raw"""
    euler_characteristic(v::AbstractNormalToricVariety)

Return the Euler characteristic of the normal toric variety `v`.

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> euler_characteristic(antv)
1
```
"""
@attr Int function euler_characteristic(v::AbstractNormalToricVariety)
    f_vector = Vector{Int}(pm_object(v).F_VECTOR)
    return f_vector[dim(v)]
end


###############################
# Setters for rings and ideals
###############################


@doc raw"""
    is_finalized(v::AbstractNormalToricVariety)

Checks if the Cox ring, the coordinate ring of the torus,
the cohomology_ring, the Chow ring, the Stanley-Reisner ideal,
the irrelevant ideal, the ideal of linear relations
or the toric ideal has been cached. If any of these has been
cached, then this function returns `true` and otherwise `false`.

# Examples
```jldoctest
julia> is_finalized(del_pezzo_surface(NormalToricVariety, 3))
false
```
"""
function is_finalized(v::AbstractNormalToricVariety)
    properties = [
        :cox_ring,
        :coordinate_ring_of_torus,
        :cohomology_ring,
        :chow_ring,
        :stanley_reisner_ideal,
        :irrelevant_ideal,
        :ideal_of_linear_relations,
        :toric_ideal,
    ]
    return any(p -> has_attribute(v, p), properties)
end


@doc raw"""
    set_coordinate_names(v::AbstractNormalToricVariety, coordinate_names::Vector{String})

Allows to set the names of the homogeneous coordinates as long as the toric variety in
question is not yet finalized (cf. [`is_finalized(v::AbstractNormalToricVariety)`](@ref)).

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> set_coordinate_names(antv, ["u"])

julia> coordinate_names(antv)
1-element Vector{String}:
 "u"
```
"""
function set_coordinate_names(v::AbstractNormalToricVariety, coordinate_names::Vector{String})
    if is_finalized(v)
        error("The coordinate names cannot be modified since the toric variety is finalized")
    end
    @req length(coordinate_names) == nrays(v) "The provided list of coordinate names must match the number of rays in the fan"
    set_attribute!(v, :coordinate_names, coordinate_names)
end


@doc raw"""
    set_coordinate_names_of_torus(v::AbstractNormalToricVariety, coordinate_names::Vector{String})

Allows to set the names of the coordinates of the torus.

# Examples
```jldoctest
julia> F3 = hirzebruch_surface(NormalToricVariety, 3);

julia> set_coordinate_names_of_torus(F3, ["u", "v"])

julia> coordinate_names_of_torus(F3)
2-element Vector{String}:
 "u"
 "v"
```
"""
function set_coordinate_names_of_torus(v::AbstractNormalToricVariety, coordinate_names::Vector{String})
    if is_finalized(v)
        error("The coordinate names of the torus cannot be modified since the toric variety is finalized")
    end
    @req length(coordinate_names) == ambient_dim(v) "The provided list of coordinate names must match the ambient dimension of the fan"
    set_attribute!(v, :coordinate_names_of_torus, coordinate_names)
end



########################################
# Cox ring
########################################


@doc raw"""
    coefficient_ring(v::AbstractNormalToricVariety)

Return the coefficient_ring `QQ` of the normal toric variety `v`.

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> coefficient_ring(antv) == QQ
true
```
"""
coefficient_ring(v::AbstractNormalToricVariety) = QQ


@doc raw"""
    coordinate_names(v::AbstractNormalToricVariety)

Return the names of the homogeneous coordinates of
the normal toric variety `v`. The default is `x1, ..., xn`.

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> coordinate_names(antv)
1-element Vector{String}:
 "x1"
```
"""
@attr Vector{String} function coordinate_names(v::AbstractNormalToricVariety)
    return ["x$(i)" for i in 1:rank(torusinvariant_weil_divisor_group(v))]
end


function _cox_ring_weights(v::AbstractNormalToricVariety)
    return get_attribute(v, :cox_ring_weights) do
        return [map_from_torusinvariant_weil_divisor_group_to_class_group(v)(x) for x in gens(torusinvariant_weil_divisor_group(v))]
    end
end


@doc raw"""
    cox_ring(R::MPolyRing, v::AbstractNormalToricVariety)

Computes the Cox ring of the normal toric variety `v`, in this case by adding
the Cox grading to the given ring `R`.
Note that [CLS11](@cite) refers to this ring as the "total coordinate ring".

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R, _ = polynomial_ring(QQ, 3);

julia> cox_ring(R, p2)
Multivariate polynomial ring in 3 variables over QQ graded by
  x1 -> [1]
  x2 -> [1]
  x3 -> [1]
```
"""
function cox_ring(R::MPolyRing, v::AbstractNormalToricVariety)
    weights = _cox_ring_weights(v)
    @req length(weights) == nvars(R) "Wrong number of variables"
    return grade(R, weights)[1]
end


@doc raw"""
    cox_ring(v::AbstractNormalToricVariety)

Computes the Cox ring of the normal toric variety `v`.
Note that [CLS11](@cite) refers to this ring as the "total coordinate ring".

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> set_coordinate_names(p2, ["y1", "y2", "y3"])

julia> cox_ring(p2)
Multivariate polynomial ring in 3 variables over QQ graded by
  y1 -> [1]
  y2 -> [1]
  y3 -> [1]
```
"""
@attr MPolyRing function cox_ring(v::AbstractNormalToricVariety)
    S, _ = polynomial_ring(coefficient_ring(v), coordinate_names(v), cached=false)
    return cox_ring(S, v)
end


########################################
# Stanley-Reisner ideal
########################################


@attr Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric} function _minimal_nonfaces(v::AbstractNormalToricVariety)
    I = ray_indices(maximal_cones(v))
    K = SimplicialComplex(I)
    return minimal_nonfaces(IncidenceMatrix, K)
end


@doc raw"""
    stanley_reisner_ideal(R::MPolyRing, v::AbstractNormalToricVariety)

Return the Stanley-Reisner ideal of a normal toric variety `v` as an ideal of
`R`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R, _ = polynomial_ring(QQ, 3);

julia> ngens(stanley_reisner_ideal(R, p2))
1
```
"""
function stanley_reisner_ideal(R::MPolyRing, v::AbstractNormalToricVariety)
    n = nrays(v)
    @req n == nvars(R) "Wrong number of variables"
    mnf = _minimal_nonfaces(v)
    Polymake.nrows(mnf) > 0 || return ideal([zero(R)])
    return ideal([ R([1], [Vector{Int}(mnf[i, :])]) for i in 1:Polymake.nrows(mnf) ])
end


@doc raw"""
    stanley_reisner_ideal(v::AbstractNormalToricVariety)

Return the Stanley-Reisner ideal of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> ngens(stanley_reisner_ideal(p2))
1
```
"""
@attr MPolyIdeal stanley_reisner_ideal(v::AbstractNormalToricVariety) = stanley_reisner_ideal(cox_ring(v), v)


########################################
# Irrelevant ideal
########################################


@attr Vector{Vector{Int}} function _irrelevant_ideal_monomials(v::AbstractNormalToricVariety)
    mc = ray_indices(maximal_cones(v))
    result = Vector{Vector{Int}}()
    onesv = ones(Int, Polymake.ncols(mc))
    for i in 1:Polymake.nrows(mc)
        push!(result, onesv - Vector{Int}(mc[i, :]))
    end
    return result
end


@doc raw"""
    irrelevant_ideal(R::MPolyRing, v::AbstractNormalToricVariety)

Return the irrelevant ideal of a normal toric variety `v` as an ideal in `R`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R, _ = polynomial_ring(QQ, 3);

julia> length(gens(irrelevant_ideal(R, p2)))
3
```
"""
function irrelevant_ideal(R::MPolyRing, v::AbstractNormalToricVariety)
    monoms = _irrelevant_ideal_monomials(v)
    @req nvars(R) == nrays(v) "Wrong number of variables in polynomial ring"
    return ideal([R([1], [x]) for x in monoms])
end


@doc raw"""
    irrelevant_ideal(v::AbstractNormalToricVariety)

Return the irrelevant ideal of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> length(gens(irrelevant_ideal(p2)))
3
```
"""
@attr MPolyIdeal function irrelevant_ideal(v::AbstractNormalToricVariety)
    R = cox_ring(v)
    return irrelevant_ideal(R, v)
end


########################################
# Ideal of linear relations
########################################


@doc raw"""
    ideal_of_linear_relations(R::MPolyRing, v::AbstractNormalToricVariety)

Return the ideal of linear relations of the simplicial and complete toric variety `v` in the ring R.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R, _ = polynomial_ring(QQ, 3);

julia> ngens(ideal_of_linear_relations(R, p2))
2
```
"""
function ideal_of_linear_relations(R::MPolyRing, v::AbstractNormalToricVariety)
    @req is_simplicial(v) "The ideal of linear relations is only supported for simplicial toric varieties"
    @req ngens(R) == nrays(v) "The given polynomial ring must have exactly as many indeterminates as rays for the toric variety"

    indeterminates = gens(R)
    d = rank(character_lattice(v))
    generators = [sum([rays(v)[j][i] * indeterminates[j] for j in 1:nrays(v)]) for i in 1:d]
    return ideal(generators)
end


@doc raw"""
    ideal_of_linear_relations(v::AbstractNormalToricVariety)

Return the ideal of linear relations of the simplicial and complete toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> ngens(ideal_of_linear_relations(p2))
2
```
"""
@attr MPolyIdeal function ideal_of_linear_relations(v::AbstractNormalToricVariety)
    R, _ = polynomial_ring(coefficient_ring(v), coordinate_names(v))
    weights = [1 for i in 1:ngens(R)]
    grade(R, weights)
    return ideal_of_linear_relations(R, v)
end


########################################
# Toric ideal
########################################


@doc raw"""
    toric_ideal(R::MPolyRing, antv::AffineNormalToricVariety)

Return the toric ideal defining the affine normal toric variety as an ideal in
`R`.

# Examples
Take the cone over the square at height one. The resulting toric variety has
one defining equation. In projective space this corresponds to
$\mathbb{P}^1\times\mathbb{P}^1$. Note that this cone is self-dual, the toric
ideal comes from the dual cone.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 1 0 1; 1 1 1])
Polyhedral cone in ambient dimension 3

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> R, _ = polynomial_ring(QQ, 4);

julia> toric_ideal(R, antv)
ideal(-x1*x2 + x3*x4)
```
"""
function toric_ideal(R::MPolyRing, antv::AffineNormalToricVariety)
    C = cone(pm_object(antv).WEIGHT_CONE)
    gens = pm_object(C).CONE_TORIC_IDEAL.BINOMIAL_GENERATORS
    return binomial_exponents_to_ideal(R, gens)
end


@doc raw"""
    toric_ideal(antv::AffineNormalToricVariety)

Return the toric ideal defining the affine normal toric variety.

# Examples
Take the cone over the square at height one. The resulting toric variety has
one defining equation. In projective space this corresponds to
$\mathbb{P}^1\times\mathbb{P}^1$. Note that this cone is self-dual, the toric
ideal comes from the dual cone.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 1 0 1; 1 1 1])
Polyhedral cone in ambient dimension 3

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> toric_ideal(antv)
ideal(-x1*x2 + x3*x4)
```
"""
@attr MPolyIdeal function toric_ideal(antv::AffineNormalToricVariety)
    C = cone(pm_object(antv).WEIGHT_CONE)
    n = length(hilbert_basis(C))
    R, _ = polynomial_ring(coefficient_ring(antv), n, cached=false)
    return toric_ideal(R, antv)
end


@attr MPolyIdeal function toric_ideal(ntv::NormalToricVariety)
    is_affine(ntv) || error("Cannot construct affine toric variety from non-affine input")
    return toric_ideal(affine_normal_toric_variety(ntv))
end


########################################
# Coordinate ring of torus
########################################


@doc raw"""
    coordinate_names_of_torus(v::AbstractNormalToricVariety)

Return the names of the coordinates of the torus of
the normal toric variety `v`. The default is `x1, ..., xn`.
"""
@attr Vector{String} function coordinate_names_of_torus(v::AbstractNormalToricVariety)
    return ["x$(i)" for i in 1:ambient_dim(v)]
end


@doc raw"""
    coordinate_ring_of_torus(R::MPolyRing, v::AbstractNormalToricVariety)

Computes the coordinate ring of the torus of the normal toric variety `v`
in the given polynomial ring `R`.
"""
function coordinate_ring_of_torus(R::MPolyRing, v::AbstractNormalToricVariety)
    n = length(coordinate_names_of_torus(v))
    @req ngens(R) >= 2 * n "The given ring must have at least $n indeterminates"
    relations = [gen(R, i) * gen(R, i+n) - one(coefficient_ring(R)) for i in 1:n]
    return quo(R, ideal(relations))[1]
end


@doc raw"""
    coordinate_ring_of_torus(v::AbstractNormalToricVariety)

Computes the coordinate ring of the torus of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> set_coordinate_names_of_torus(p2, ["y1", "y2"])

julia> coordinate_ring_of_torus(p2)
Quotient of Multivariate polynomial ring in 4 variables over QQ by ideal(y1*y1_ - 1, y2*y2_ - 1)
```
"""
@attr MPolyQuoRing function coordinate_ring_of_torus(v::AbstractNormalToricVariety)
    S, _ = polynomial_ring(coefficient_ring(v), vcat(coordinate_names_of_torus(v), [x*"_" for x in coordinate_names_of_torus(v)]), cached=false)
    return coordinate_ring_of_torus(S, v)
end


#########################################
# Turn characters into rational functions
#########################################


@doc raw"""
    character_to_rational_function(v::AbstractNormalToricVariety, character::Vector{ZZRingElem})

Computes the rational function corresponding to a character of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> character_to_rational_function(p2, [-1, 2])
x2^2*x1_
```
"""
function character_to_rational_function(v::AbstractNormalToricVariety, character::Vector{ZZRingElem})
    S, _ = polynomial_ring(coefficient_ring(v), vcat(coordinate_names_of_torus(v), [x*"_" for x in coordinate_names_of_torus(v)]), cached=false)
    return character_to_rational_function(S, v::AbstractNormalToricVariety, character::Vector{ZZRingElem})
end
character_to_rational_function(v::AbstractNormalToricVariety, character::Vector{Int}) = character_to_rational_function(v, [ZZRingElem(k) for k in character])


@doc raw"""
    character_to_rational_function(R::MPolyRing, v::AbstractNormalToricVariety, character::Vector{ZZRingElem})

Computes the rational function corresponding to a character of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R, _ = polynomial_ring(QQ, 4);

julia> character_to_rational_function(R, p2, [-1, 2])
x2^2*x3
```
"""
function character_to_rational_function(R::MPolyRing, v::AbstractNormalToricVariety, character::Vector{ZZRingElem})
    if ambient_dim(v) != length(character)
        throw(ArgumentError("A character consist of as many integers as the ambient dimension of the variety"))
    end
    generators = gens(coordinate_ring_of_torus(R, v))
    rational_function = one(coefficient_ring(R))
    for i in 1:length(character)
        if character[i] < 0
            rational_function *= generators[i+ambient_dim(v)]^(-1*character[i])
        else
            rational_function *= generators[i]^(character[i])
        end
    end
    return rational_function
end
character_to_rational_function(R::MPolyRing, v::AbstractNormalToricVariety, character::Vector{Int}) = character_to_rational_function(R, v, [ZZRingElem(k) for k in character])


##############################################
# Characters, Weil divisor and the class group
##############################################


@doc raw"""
    character_lattice(v::AbstractNormalToricVariety)

Return the character lattice of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> character_lattice(p2)
GrpAb: Z^2
```
"""
@attr GrpAbFinGen character_lattice(v::AbstractNormalToricVariety) = free_abelian_group(ambient_dim(v))


@doc raw"""
    torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)

Return the torusinvariant divisor group of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> torusinvariant_weil_divisor_group(p2)
GrpAb: Z^3
```
"""
@attr GrpAbFinGen torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety) = free_abelian_group(nrays(v))


@doc raw"""
    map_from_character_lattice_to_torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)

Return the map from the character lattice to the group of principal divisors of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> map_from_character_lattice_to_torusinvariant_weil_divisor_group(p2)
Map with following data
Domain:
=======
Abelian group with structure: Z^2
Codomain:
=========
Abelian group with structure: Z^3
```
"""
@attr GrpAbFinGenMap function map_from_character_lattice_to_torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)
    mat = transpose(matrix(ZZ, rays(v)))
    return hom(character_lattice(v), torusinvariant_weil_divisor_group(v), mat)
end


@doc raw"""
    torusinvariant_prime_divisors(v::AbstractNormalToricVariety)

Return the list of all torus invariant prime divisors in a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> torusinvariant_prime_divisors(p2)
3-element Vector{ToricDivisor}:
 Torus-invariant, prime divisor on a normal toric variety
 Torus-invariant, prime divisor on a normal toric variety
 Torus-invariant, prime divisor on a normal toric variety
```
"""
@attr Vector{ToricDivisor} function torusinvariant_prime_divisors(v::AbstractNormalToricVariety)
    ti_divisors = torusinvariant_weil_divisor_group(v)
    prime_divisors = ToricDivisor[]
    for i in 1:rank(ti_divisors)
        coeffs = zeros(Int, rank(ti_divisors))
        coeffs[i] = 1
        push!(prime_divisors, toric_divisor(v, coeffs))
    end
    return prime_divisors
end


@doc raw"""
    class_group(v::AbstractNormalToricVariety)

Return the class group of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> class_group(p2)
GrpAb: Z
```
"""
@attr GrpAbFinGen class_group(v::AbstractNormalToricVariety) = codomain(map_from_torusinvariant_weil_divisor_group_to_class_group(v))


@doc raw"""
    map_from_torusinvariant_weil_divisor_group_to_class_group(v::AbstractNormalToricVariety)

Return the map from the group of Weil divisors to the class of group of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> map_from_torusinvariant_weil_divisor_group_to_class_group(p2)
Map with following data
Domain:
=======
Abelian group with structure: Z^3
Codomain:
=========
Abelian group with structure: Z
```
"""
@attr GrpAbFinGenMap function map_from_torusinvariant_weil_divisor_group_to_class_group(v::AbstractNormalToricVariety)
    map1 = cokernel(map_from_character_lattice_to_torusinvariant_weil_divisor_group(v))[2]
    map2 = inv(snf(codomain(map1))[2])
    return map1*map2
end


@doc raw"""
    map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)

Return the embedding of the group of Cartier divisors into the group of
torus-invariant Weil divisors of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(p2)
Map with following data
Domain:
=======
Abelian group with structure: Z^3
Codomain:
=========
Abelian group with structure: Z^3
```
"""
@attr Map{GrpAbFinGen, GrpAbFinGen} function map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)
    # check input
    if has_torusfactor(v)
        throw(ArgumentError("Group of the torus-invariant Cartier divisors can only be computed if the variety has no torus factor"))
    end
    
    # identify fan_rays and cones
    fan_rays = transpose(matrix(ZZ, rays(v)))
    max_cones = ray_indices(maximal_cones(v))
    number_of_rays = ncols(fan_rays)
    number_of_cones = size(max_cones)[1]
    
    # compute quantities needed to construct the matrices
    rc = rank(character_lattice(v))
    number_ray_is_part_of_max_cones = [length(max_cones[:, k].s) for k in 1:number_of_rays]
    s = sum(number_ray_is_part_of_max_cones)
    cones_ray_is_part_of = [filter(x -> max_cones[x, r], 1:number_of_cones) for r in 1:number_of_rays]
    
    # compute the matrix for the scalar products
    map_for_scalar_products = zero_matrix(ZZ, number_of_cones * rc, s)
    col = 1
    for i in 1:number_of_rays
        for j in cones_ray_is_part_of[i]
            map_for_scalar_products[(j-1)*rc+1:j*rc, col] = [ZZRingElem(c) for c in fan_rays[:, i]]
            col += 1
        end
    end
    
    # compute the matrix for differences
    map_for_difference_of_elements = zero_matrix(ZZ, s, s-number_of_rays)
    row = 1
    col = 1
    for i in 1:number_of_rays
        ncol = number_ray_is_part_of_max_cones[i]-1
        map_for_difference_of_elements[row, col:(col+ncol-1)] = fill(ZZRingElem(1), ncol)
        map_for_difference_of_elements[(row+1):(row+ncol), col:(col+ncol-1)] = -identity_matrix(ZZ, ncol)
        row += ncol + 1
        col += ncol
    end
    
    # compute the matrix for mapping to torusinvariant Weil divisors
    map_to_weil_divisors = zero_matrix(ZZ, number_of_cones * rc, rank(torusinvariant_weil_divisor_group(v)))
    for i in 1:number_of_rays
        map_to_weil_divisors[(cones_ray_is_part_of[i][1]-1)*rc+1:cones_ray_is_part_of[i][1]*rc, i] = [ZZRingElem(-c) for c in fan_rays[:, i]]
    end
    
    # compute the total map
    mapping_matrix = map_for_scalar_products * map_for_difference_of_elements
    source = free_abelian_group(nrows(mapping_matrix))
    target = free_abelian_group(ncols(mapping_matrix))
    total_map = hom(source, target, mapping_matrix)

    # identify the embedding of the cartier_data_group
    ker = kernel(total_map)
    embedding = snf(ker[1])[2] * ker[2] * hom(codomain(ker[2]), torusinvariant_weil_divisor_group(v), map_to_weil_divisors)
    
    # return the image of this embedding
    return image(embedding)[2]
end


@doc raw"""
    torusinvariant_cartier_divisor_group(v::AbstractNormalToricVariety)

Return the Cartier divisor group of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> torusinvariant_cartier_divisor_group(p2)
GrpAb: Z^3
```
"""
@attr GrpAbFinGen function torusinvariant_cartier_divisor_group(v::AbstractNormalToricVariety)
    return domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v))
end


@doc raw"""
    map_from_torusinvariant_cartier_divisor_group_to_picard_group(v::AbstractNormalToricVariety)

Return the map from the Cartier divisors to the Picard group
of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> map_from_torusinvariant_cartier_divisor_group_to_picard_group(p2)
Map with following data
Domain:
=======
Abelian group with structure: Z^3
Codomain:
=========
Abelian group with structure: Z
```
"""
@attr GrpAbFinGenMap function map_from_torusinvariant_cartier_divisor_group_to_picard_group(v::AbstractNormalToricVariety)
    # check input
    if has_torusfactor(v)
        throw(ArgumentError("Group of the torus-invariant Cartier divisors can only be computed if the variety has no torus factor"))
    end
    
    # compute mapping
    map1 = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v)
    map2 = map_from_torusinvariant_weil_divisor_group_to_class_group(v)
    return restrict_codomain(map1*map2)
end


@doc raw"""
    picard_group(v::AbstractNormalToricVariety)

Return the Picard group of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> picard_group(p2)
GrpAb: Z
```
"""
@attr GrpAbFinGen function picard_group(v::AbstractNormalToricVariety)
    return codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(v))
end


############################
# Cones and fans
############################

"""
    nef_cone(v::NormalToricVariety)

Return the nef cone of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> nef = nef_cone(p2)
Polyhedral cone in ambient dimension 1

julia> dim(nef)
1
```
"""
@attr Cone nef_cone(v::NormalToricVariety) = cone(pm_object(v).NEF_CONE)


"""
    mori_cone(v::NormalToricVariety)

Return the mori cone of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> mori = mori_cone(p2)
Polyhedral cone in ambient dimension 1

julia> dim(mori)
1
```
"""
@attr Cone mori_cone(v::NormalToricVariety) = cone(pm_object(v).MORI_CONE)


@doc raw"""
    fan(v::AbstractNormalToricVariety)

Return the fan of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> fan(p2)
Polyhedral fan in ambient dimension 2
```
"""
@attr PolyhedralFan{QQFieldElem} fan(v::AbstractNormalToricVariety) = PolyhedralFan{QQFieldElem}(pm_object(v))


@doc raw"""
    cone(v::AffineNormalToricVariety)

Return the cone of the affine normal toric variety `v`.

# Examples
```jldoctest
julia> cone(affine_normal_toric_variety(Oscar.positive_hull([1 1; -1 1])))
Polyhedral cone in ambient dimension 2
```
"""
@attr Cone cone(v::AffineNormalToricVariety) = maximal_cones(v)[1]


@doc raw"""
    dual_cone(v::AffineNormalToricVariety)

Return the dual cone of the affine normal toric variety `v`.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> dual_cone(antv)
Polyhedral cone in ambient dimension 2

julia> polarize(cone(antv)) == dual_cone(antv)
true
```
"""
@attr Cone dual_cone(v::AffineNormalToricVariety) = polarize(cone(v))


@doc raw"""
    hilbert_basis(v::AffineNormalToricVariety)

For an affine toric variety ``v``, this returns the Hilbert
basis of the cone dual to the cone of ``v``.

# Examples
```jldoctest
julia> C = positive_hull([-1 1; 1 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> hilbert_basis(antv)
[-1   1]
[ 1   1]
[ 0   1]
```
"""
@attr ZZMatrix hilbert_basis(v::AffineNormalToricVariety) = matrix(ZZ, hilbert_basis(dual_cone(v)))


_variable_ray_correspondence(v::AbstractNormalToricVariety) = Dict{RayVector, MPolyRingElem}(zip(rays(v), gens(cox_ring(v))))
_ray_variable_correspondence(v::AbstractNormalToricVariety) = Dict{MPolyRingElem, RayVector}(zip(gens(cox_ring(v)), rays(v)))


############################
# Affine covering
############################


@doc raw"""
    affine_open_covering(v::AbstractNormalToricVariety)

Compute an affine open cover of the normal toric variety `v`, i.e. returns a list of affine toric varieties.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> affine_open_covering(p2)
3-element Vector{AffineNormalToricVariety}:
 Normal, affine toric variety
 Normal, affine toric variety
 Normal, affine toric variety
```
"""
@attr Vector{AffineNormalToricVariety} function affine_open_covering(v::AbstractNormalToricVariety)
    charts = Vector{AffineNormalToricVariety}(undef, pm_object(v).N_MAXIMAL_CONES)
    for i in 1:pm_object(v).N_MAXIMAL_CONES
        charts[i] = affine_normal_toric_variety(cone(Polymake.fan.cone(pm_object(v), i-1)))
    end
    return charts
end
