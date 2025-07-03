############################
# Dimensions
############################


@doc raw"""
    dim(v::NormalToricVarietyType)

Return the dimension of the normal toric variety `v`.

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> dim(antv)
1
```
"""
@attr Int dim(v::NormalToricVarietyType) = pm_object(v).FAN_DIM


@doc raw"""
    dim_of_torusfactor(v::NormalToricVarietyType)

Return the dimension of the torus factor of the normal toric variety `v`.

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> dim_of_torusfactor(antv)
1
```
"""
@attr Int function dim_of_torusfactor(v::NormalToricVarietyType)
    if has_torusfactor(v) == false
        return 0
    end
    dimension_of_fan = pm_object(v).FAN_DIM::Int
    ambient_dimension = pm_object(v).FAN_AMBIENT_DIM::Int
    return ambient_dimension - dimension_of_fan
end


@doc raw"""
    euler_characteristic(v::NormalToricVarietyType)

Return the Euler characteristic of the normal toric variety `v`.

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> euler_characteristic(antv)
1
```
"""
@attr Int function euler_characteristic(v::NormalToricVarietyType)
    f_vector = Vector{Int}(pm_object(v).F_VECTOR)
    return f_vector[dim(v)]
end


###############################
# Setters for rings and ideals
###############################


@doc raw"""
    is_finalized(v::NormalToricVarietyType)

Check if the Cox ring, the coordinate ring of the torus,
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
function is_finalized(v::NormalToricVarietyType)
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
    set_coordinate_names(v::NormalToricVarietyType, coordinate_names::AbstractVector{<:VarName})

Allows to set the names of the homogeneous coordinates as long as the toric variety in
question is not yet finalized (cf. [`is_finalized(v::NormalToricVarietyType)`](@ref)).

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> set_coordinate_names(antv, [:u])

julia> coordinate_names(antv)
1-element Vector{String}:
 "u"

julia> set_coordinate_names(antv, ["v"])

julia> coordinate_names(antv)
1-element Vector{String}:
 "v"

julia> set_coordinate_names(antv, ['w'])

julia> coordinate_names(antv)
1-element Vector{String}:
 "w"
```
"""
function set_coordinate_names(v::NormalToricVarietyType, coordinate_names::AbstractVector{<:VarName})
    if is_finalized(v)
        error("The coordinate names cannot be modified since the toric variety is finalized")
    end
    @req length(coordinate_names) == n_rays(v) "The provided list of coordinate names must match the number of rays in the fan"
    set_attribute!(v, :coordinate_names, string.(coordinate_names))
end


@doc raw"""
    set_coordinate_names_of_torus(v::NormalToricVarietyType, coordinate_names::AbstractVector{<:VarName})

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
function set_coordinate_names_of_torus(v::NormalToricVarietyType, coordinate_names::AbstractVector{<:VarName})
    if is_finalized(v)
        error("The coordinate names of the torus cannot be modified since the toric variety is finalized")
    end
    @req length(coordinate_names) == ambient_dim(v) "The provided list of coordinate names must match the ambient dimension of the fan"
    set_attribute!(v, :coordinate_names_of_torus, string.(coordinate_names))
end



########################################
# Cox ring
########################################


@doc raw"""
    coefficient_ring(v::NormalToricVarietyType)

Return the coefficient_ring `QQ` of the normal toric variety `v`.

# Examples
```jldoctest
julia> C = Oscar.positive_hull([1 0]);

julia> antv = affine_normal_toric_variety(C);

julia> coefficient_ring(antv) == QQ
true
```
"""
coefficient_ring(v::NormalToricVarietyType) = QQ


@doc raw"""
    coordinate_names(v::NormalToricVarietyType)

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
@attr Vector{String} function coordinate_names(v::NormalToricVarietyType)
  if has_attribute(v, :cox_ring)
    return string.(symbols(cox_ring(v)))
  end
  return ["x$(i)" for i in 1:torsion_free_rank(torusinvariant_weil_divisor_group(v))]
end


@attr Vector{FinGenAbGroupElem} function _cox_ring_weights(v::NormalToricVarietyType)
  f = map_from_torusinvariant_weil_divisor_group_to_class_group(v)
  return [f(x) for x in gens(torusinvariant_weil_divisor_group(v))]
end


@doc raw"""
    cox_ring(R::MPolyRing, v::NormalToricVarietyType)

Compute the Cox ring of the normal toric variety `v`, in this case by adding
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
function cox_ring(R::MPolyRing, v::NormalToricVarietyType)
    weights = _cox_ring_weights(v)
    @req length(weights) == nvars(R) "Wrong number of variables"
    return grade(R, weights)[1]
end


@doc raw"""
    cox_ring(v::NormalToricVarietyType)

Compute the Cox ring of the normal toric variety `v`.
Note that [CLS11](@cite) refers to this ring as the "total coordinate ring".
For uniformity with schemes, we also support the function
`coordinate_ring` to refer to the Cox ring.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> set_coordinate_names(p2, ["y1", "y2", "y3"])

julia> cox_ring(p2)
Multivariate polynomial ring in 3 variables over QQ graded by
  y1 -> [1]
  y2 -> [1]
  y3 -> [1]

julia> cox_ring(p2) == coordinate_ring(p2)
true
```
"""
@attr MPolyRing function cox_ring(v::NormalToricVarietyType)
    S, _ = polynomial_ring(coefficient_ring(v), coordinate_names(v); cached=false)
    return cox_ring(S, v)
end

coordinate_ring(v::NormalToricVarietyType) = cox_ring(v)

########################################
# Stanley-Reisner ideal
########################################


@attr Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric} function _minimal_nonfaces(v::NormalToricVarietyType)
    I = ray_indices(maximal_cones(v))
    K = simplicial_complex(I)
    return minimal_nonfaces(IncidenceMatrix, K)
end


@doc raw"""
    stanley_reisner_ideal(R::MPolyRing, v::NormalToricVarietyType)

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
function stanley_reisner_ideal(R::MPolyRing, v::NormalToricVarietyType)
    n = n_rays(v)
    @req n == nvars(R) "Wrong number of variables"
    mnf = _minimal_nonfaces(v)
    Polymake.nrows(mnf) > 0 || return ideal([zero(R)])
    return ideal([ R([1], [Vector{Int}(mnf[i, :])]) for i in 1:Polymake.nrows(mnf) ])
end


@doc raw"""
    stanley_reisner_ideal(v::NormalToricVarietyType)

Return the Stanley-Reisner ideal of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> ngens(stanley_reisner_ideal(p2))
1
```
"""
@attr MPolyIdeal stanley_reisner_ideal(v::NormalToricVarietyType) = stanley_reisner_ideal(cox_ring(v), v)


########################################
# Irrelevant ideal
########################################


@attr Vector{Vector{Int}} function _irrelevant_ideal_monomials(v::NormalToricVarietyType)
    mc = ray_indices(maximal_cones(v))
    result = Vector{Vector{Int}}()
    onesv = ones(Int, Polymake.ncols(mc))
    for i in 1:Polymake.nrows(mc)
        push!(result, onesv - Vector{Int}(mc[i, :]))
    end
    return result
end


@doc raw"""
    irrelevant_ideal(R::MPolyRing, v::NormalToricVarietyType)

Return the irrelevant ideal of a normal toric variety `v` as an ideal in `R`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R, _ = polynomial_ring(QQ, 3);

julia> length(gens(irrelevant_ideal(R, p2)))
3
```
"""
function irrelevant_ideal(R::MPolyRing, v::NormalToricVarietyType)
    monoms = _irrelevant_ideal_monomials(v)
    @req nvars(R) == n_rays(v) "Wrong number of variables in polynomial ring"
    return ideal([R([1], [x]) for x in monoms])
end


@doc raw"""
    irrelevant_ideal(v::NormalToricVarietyType)

Return the irrelevant ideal of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> length(gens(irrelevant_ideal(p2)))
3
```
"""
@attr MPolyIdeal function irrelevant_ideal(v::NormalToricVarietyType)
    R = cox_ring(v)
    return irrelevant_ideal(R, v)
end


########################################
# Ideal of linear relations
########################################


@doc raw"""
    ideal_of_linear_relations(R::MPolyRing, v::NormalToricVarietyType)

Return the ideal of linear relations of the simplicial and complete toric variety `v` in the ring R.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R, _ = polynomial_ring(QQ, 3);

julia> ngens(ideal_of_linear_relations(R, p2))
2
```
"""
function ideal_of_linear_relations(R::MPolyRing, v::NormalToricVarietyType)
    @req is_simplicial(v) "The ideal of linear relations is only supported for simplicial toric varieties"
    @req ngens(R) == n_rays(v) "The given polynomial ring must have exactly as many indeterminates as rays for the toric variety"
    return ideal(transpose(matrix(ZZ, rays(v))) * gens(R))
end


@doc raw"""
    ideal_of_linear_relations(v::NormalToricVarietyType)

Return the ideal of linear relations of the simplicial and complete toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> ngens(ideal_of_linear_relations(p2))
2
```
"""
@attr MPolyIdeal function ideal_of_linear_relations(v::NormalToricVarietyType)
    R, _ = graded_polynomial_ring(coefficient_ring(v), coordinate_names(v); cached=false)
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
Normal toric variety

julia> R, _ = polynomial_ring(QQ, 4);

julia> toric_ideal(R, antv)
Ideal generated by
  -x1*x2 + x3*x4
```
"""
function toric_ideal(R::MPolyRing, antv::AffineNormalToricVariety)
  C = cone(pm_object(antv).WEIGHT_CONE)
  HB = matrix(ZZ, hilbert_basis(C))
  # If the cone is smooth, the Hilbert basis is itself a basis and the toric
  # ideal is zero. We catch this here as 4ti2 is unhappy about some of the
  # smooth cones.
  if rank(HB) == nrows(HB)
    return ideal(R, [R(0)])
  else
    gens = pm_object(C).CONE_TORIC_IDEAL.BINOMIAL_GENERATORS
    return binomial_exponents_to_ideal(R, gens)
  end
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
Normal toric variety

julia> toric_ideal(antv)
Ideal generated by
  -x1*x2 + x3*x4
```
"""
@attr MPolyIdeal function toric_ideal(antv::AffineNormalToricVariety)
    C = cone(pm_object(antv).WEIGHT_CONE)
    n = length(hilbert_basis(C))
    R, _ = polynomial_ring(coefficient_ring(antv), n; cached=false)
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
    coordinate_names_of_torus(v::NormalToricVarietyType)

Return the names of the coordinates of the torus of
the normal toric variety `v`. The default is `x1, ..., xn`.
"""
@attr Vector{String} function coordinate_names_of_torus(v::NormalToricVarietyType)
    return ["x$(i)" for i in 1:ambient_dim(v)]
end


@doc raw"""
    coordinate_ring_of_torus(R::MPolyRing, v::NormalToricVarietyType)

Compute the coordinate ring of the torus of the normal toric variety `v`
in the given polynomial ring `R`.
"""
function coordinate_ring_of_torus(R::MPolyRing, v::NormalToricVarietyType)
    n = length(coordinate_names_of_torus(v))
    @req ngens(R) >= 2 * n "The given ring must have at least $n indeterminates"
    relations = [gen(R, i) * gen(R, i+n) - one(coefficient_ring(R)) for i in 1:n]
    return quo(R, ideal(relations))[1]
end


@doc raw"""
    coordinate_ring_of_torus(v::NormalToricVarietyType)

Compute the coordinate ring of the torus of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> set_coordinate_names_of_torus(p2, ["y1", "y2"])

julia> coordinate_ring_of_torus(p2)
Quotient
  of multivariate polynomial ring in 4 variables y1, y2, y1_, y2_
    over rational field
  by ideal (y1*y1_ - 1, y2*y2_ - 1)
```
"""
@attr MPolyQuoRing function coordinate_ring_of_torus(v::NormalToricVarietyType)
    S, _ = polynomial_ring(coefficient_ring(v), vcat(coordinate_names_of_torus(v), [x*"_" for x in coordinate_names_of_torus(v)]); cached=false)
    return coordinate_ring_of_torus(S, v)
end


#########################################
# Turn characters into rational functions
#########################################


@doc raw"""
    character_to_rational_function(v::NormalToricVarietyType, character::Vector{ZZRingElem})

Compute the rational function corresponding to a character of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> character_to_rational_function(p2, [-1, 2])
x2^2*x1_
```
"""
function character_to_rational_function(v::NormalToricVarietyType, character::Vector{ZZRingElem})
    S, _ = polynomial_ring(coefficient_ring(v), vcat(coordinate_names_of_torus(v), [x*"_" for x in coordinate_names_of_torus(v)]); cached=false)
    return character_to_rational_function(S, v::NormalToricVarietyType, character::Vector{ZZRingElem})
end
character_to_rational_function(v::NormalToricVarietyType, character::Vector{Int}) = character_to_rational_function(v, [ZZRingElem(k) for k in character])


@doc raw"""
    character_to_rational_function(R::MPolyRing, v::NormalToricVarietyType, character::Vector{ZZRingElem})

Compute the rational function corresponding to a character of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R, _ = polynomial_ring(QQ, 4);

julia> character_to_rational_function(R, p2, [-1, 2])
x2^2*x3
```
"""
function character_to_rational_function(R::MPolyRing, v::NormalToricVarietyType, character::Vector{ZZRingElem})
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
character_to_rational_function(R::MPolyRing, v::NormalToricVarietyType, character::Vector{Int}) = character_to_rational_function(R, v, [ZZRingElem(k) for k in character])


##############################################
# Characters, Weil divisor and the class group
##############################################


@doc raw"""
    character_lattice(v::NormalToricVarietyType)

Return the character lattice of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> character_lattice(p2)
Z^2
```
"""
@attr FinGenAbGroup character_lattice(v::NormalToricVarietyType) = free_abelian_group(ambient_dim(v))


@doc raw"""
    lattice_of_one_parameter_subgroups(v::NormalToricVarietyType)

Return the lattice of one parameter subgroups of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> lattice_of_one_parameter_subgroups(p2)
Z^2
```
"""
@attr FinGenAbGroup lattice_of_one_parameter_subgroups(v::NormalToricVarietyType) = free_abelian_group(ambient_dim(v))


@doc raw"""
    torusinvariant_weil_divisor_group(v::NormalToricVarietyType)

Return the torusinvariant divisor group of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> torusinvariant_weil_divisor_group(p2)
Z^3
```
"""
@attr FinGenAbGroup torusinvariant_weil_divisor_group(v::NormalToricVarietyType) = free_abelian_group(n_rays(v))


@doc raw"""
    map_from_character_lattice_to_torusinvariant_weil_divisor_group(v::NormalToricVarietyType)

Return the map from the character lattice to the group of principal divisors of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> map_from_character_lattice_to_torusinvariant_weil_divisor_group(p2)
Map
  from Z^2
  to Z^3
```
"""
@attr FinGenAbGroupHom function map_from_character_lattice_to_torusinvariant_weil_divisor_group(v::NormalToricVarietyType)
    mat = transpose(matrix(ZZ, rays(v)))
    return hom(character_lattice(v), torusinvariant_weil_divisor_group(v), mat)
end


@doc raw"""
    torusinvariant_prime_divisors(v::NormalToricVarietyType)

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
@attr Vector{ToricDivisor} function torusinvariant_prime_divisors(v::NormalToricVarietyType)
    ti_divisors = torusinvariant_weil_divisor_group(v)
    prime_divisors = ToricDivisor[]
    for i in 1:torsion_free_rank(ti_divisors)
        coeffs = zeros(Int, torsion_free_rank(ti_divisors))
        coeffs[i] = 1
        push!(prime_divisors, toric_divisor(v, coeffs))
    end
    return prime_divisors
end


@doc raw"""
    class_group(v::NormalToricVarietyType)

Return the class group of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> class_group(p2)
Z
```
"""
class_group(v::NormalToricVarietyType) = class_group_with_map(v)[1]


@doc raw"""
    class_group_with_map(v::NormalToricVarietyType)

Return the class group of the normal toric variety `v` together with the map from the torus
invariant Weil Divisor group into the class group.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> class_group_with_map(p2)
(Z, Map: Z^3 -> Z)
```
"""
@attr Tuple{FinGenAbGroup, FinGenAbGroupHom} function class_group_with_map(v::NormalToricVarietyType)
  map1 = cokernel(map_from_character_lattice_to_torusinvariant_weil_divisor_group(v))[2]
  map2 = inv(snf(codomain(map1))[2])
  map_into_cl = map1*map2  
  # Until August 2025, class_group was an attributed and has been serialized. The following tries to ensure consistency in this case.
  # If there were no .mrdi files with :class_group serialized, then it should be possible to remove the following if block without causing any errors.
  if has_attribute(v, :class_group)
    cg = get_attribute(v, :class_group)
    map_into_cl = map_into_cl * hom(codomain(map2), cg, gens(cg))
  end
  return codomain(map_into_cl), map_into_cl
end

@doc raw"""
    map_from_torusinvariant_weil_divisor_group_to_class_group(v::NormalToricVarietyType)

Return the map from the group of Weil divisors to the class of group of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> map_from_torusinvariant_weil_divisor_group_to_class_group(p2)
Map
  from Z^3
  to Z
```
"""
map_from_torusinvariant_weil_divisor_group_to_class_group(v::NormalToricVarietyType) = class_group_with_map(v)[2]


@doc raw"""
    map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v::NormalToricVarietyType)

Return the embedding of the group of Cartier divisors into the group of
torus-invariant Weil divisors of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(p2)
Map
  from Z^3
  to Z^3
```
"""
@attr Map{FinGenAbGroup, FinGenAbGroup} function map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v::NormalToricVarietyType)
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
    rc = torsion_free_rank(character_lattice(v))
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
    map_to_weil_divisors = zero_matrix(ZZ, number_of_cones * rc, torsion_free_rank(torusinvariant_weil_divisor_group(v)))
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
    torusinvariant_cartier_divisor_group(v::NormalToricVarietyType)

Return the Cartier divisor group of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> torusinvariant_cartier_divisor_group(p2)
Z^3
```
"""
@attr FinGenAbGroup function torusinvariant_cartier_divisor_group(v::NormalToricVarietyType)
    return domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v))
end

@doc raw"""
    map_from_torusinvariant_cartier_divisor_group_to_class_group(v::NormalToricVarietyType)

Return the map from the Cartier divisors to the class group
of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> map_from_torusinvariant_cartier_divisor_group_to_class_group(p2)
Map
  from Z^3
  to Z
```
"""
@attr FinGenAbGroupHom function map_from_torusinvariant_cartier_divisor_group_to_class_group(v::NormalToricVarietyType)
    # check input
    @req !has_torusfactor(v) "Group of the torus-invariant Cartier divisors can only be computed if the variety has no torus factor"

    f = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v)
    g = map_from_torusinvariant_weil_divisor_group_to_class_group(v)
    return f * g
end


@doc raw"""
    map_from_torusinvariant_cartier_divisor_group_to_picard_group(v::NormalToricVarietyType)

Return the map from the Cartier divisors to the Picard group
of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> map_from_torusinvariant_cartier_divisor_group_to_picard_group(p2)
Map
  from Z^3
  to Z
```
"""
@attr FinGenAbGroupHom function map_from_torusinvariant_cartier_divisor_group_to_picard_group(v::NormalToricVarietyType)
    # check input
    if has_torusfactor(v)
        throw(ArgumentError("Group of the torus-invariant Cartier divisors can only be computed if the variety has no torus factor"))
    end
    
    f = restrict_codomain(map_from_torusinvariant_cartier_divisor_group_to_class_group(v))
    return f * inv(snf(codomain(f))[2])
end


@doc raw"""
    picard_group(v::NormalToricVarietyType)

Return the Picard group of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> picard_group(p2)
Z
```
"""
@attr FinGenAbGroup function picard_group(v::NormalToricVarietyType)
    return codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(v))
end

@doc raw"""
    map_from_picard_group_to_class_group(v::NormalToricVarietyType)

Return the embedding of the Picard group into the class group of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> map_from_picard_group_to_class_group(p2)
Map
  from Z
  to Z
```
"""
@attr FinGenAbGroupHom function map_from_picard_group_to_class_group(v::NormalToricVarietyType)
    f = image(map_from_torusinvariant_cartier_divisor_group_to_class_group(v))[2]
    g = snf(domain(f))[2] * f
    cl = class_group_with_map(v)[1]
    return hom(picard_group(v), cl, matrix(g))
end


############################
# Gorenstein and Picard Index
############################

@doc raw"""
    gorenstein_index(v::NormalToricVarietyType)

Return the Gorenstein index of a $\mathbb{Q}$-Gorenstein normal toric variety `v`. 
This is the smallest positive integer $l$ such that $-l K$ is Cartier, where $K$
is a canonical divisor on `v`. See exercise 8.3.10 and 8.3.11 in [CLS11](@cite) for more details.

# Examples
```jldoctest
julia> gorenstein_index(weighted_projective_space(NormalToricVariety, [2,3,5]))
3
```
"""
@attr ZZRingElem function gorenstein_index(v::NormalToricVarietyType)
    @req is_q_gorenstein(v) "gorenstein index can only be computed for Q-gorenstein varieties"
    c = divisor_class(canonical_divisor_class(v))
    f = cokernel(map_from_picard_group_to_class_group(v))[2]
    order(f(c))
end


@doc raw"""
    picard_index(v::NormalToricVarietyType)

Return the index of the Picard group in the class group of a simplicial normal 
toric variety `v`. Here, the Picard group embeds as the group of Cartier divisor
classes into the class group via `map_from_picard_group_to_class_group`. See 
[HHS11](@cite) for more details.

# Examples
```jldoctest
julia> picard_index(weighted_projective_space(NormalToricVariety, [2,3,5]))
30
```
"""
@attr ZZRingElem function picard_index(v::NormalToricVarietyType) 
    @req is_simplicial(v) "picard index can only be computed for simplicial varieties"
    return order(cokernel(map_from_picard_group_to_class_group(v))[1]) 
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
Normal toric variety

julia> nef = nef_cone(p2)
Polyhedral cone in ambient dimension 1

julia> dim(nef)
1
```
"""
@attr Cone{QQFieldElem} function nef_cone(v::NormalToricVariety)
  result = cone(pm_object(v).NEF_CONE)
  oscar_projection = map_from_torusinvariant_weil_divisor_group_to_class_group(v).map
  polymake_lift = matrix(ZZ, pm_object(v).RATIONAL_DIVISOR_CLASS_GROUP.LIFTING)
  A = convert(Polymake.PolymakeType, transpose(polymake_lift * oscar_projection))
  return transform(result, A)
end


"""
    mori_cone(v::NormalToricVariety)

Return the mori cone of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> mori = mori_cone(p2)
Polyhedral cone in ambient dimension 1

julia> dim(mori)
1
```
"""
@attr Cone{QQFieldElem} mori_cone(v::NormalToricVariety) = polarize(nef_cone(v))


@doc raw"""
    polyhedral_fan(v::NormalToricVarietyType)

Return the fan of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> polyhedral_fan(p2)
Polyhedral fan in ambient dimension 2
```
"""
function polyhedral_fan(v::NormalToricVarietyType)
  result = Base.deepcopy(pm_object(v))
  Polymake.cast!(result, Polymake.BigObjectType("fan::PolyhedralFan<Rational>"))
  return PolyhedralFan{QQFieldElem}(result, QQ)
end


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
    weight_cone(v::AffineNormalToricVariety)

Return the dual cone of the affine normal toric variety `v`.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal toric variety

julia> weight_cone(antv)
Polyhedral cone in ambient dimension 2

julia> polarize(cone(antv)) == weight_cone(antv)
true
```
"""
@attr Cone weight_cone(v::AffineNormalToricVariety) = polarize(cone(v))


@doc raw"""
    hilbert_basis(v::AffineNormalToricVariety)

For an affine toric variety ``v``, this returns the Hilbert
basis of the cone dual to the cone of ``v``.

# Examples
```jldoctest
julia> C = positive_hull([-1 1; 1 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal toric variety

julia> hilbert_basis(antv)
[-1   1]
[ 1   1]
[ 0   1]
```
"""
@attr ZZMatrix function hilbert_basis(v::AffineNormalToricVariety)
  @req is_pointed(weight_cone(v)) "Weight cone is not pointed"
  @req is_fulldimensional(weight_cone(v)) "Weight cone is not full dimensional"
  return matrix(ZZ, hilbert_basis(weight_cone(v)))
end


_variable_ray_correspondence(v::NormalToricVarietyType) = Dict{RayVector, MPolyRingElem}(zip(rays(v), gens(cox_ring(v))))
_ray_variable_correspondence(v::NormalToricVarietyType) = Dict{MPolyRingElem, RayVector}(zip(gens(cox_ring(v)), rays(v)))


############################
# Affine covering
############################


@doc raw"""
    affine_open_covering(v::NormalToricVarietyType)

Compute an affine open cover of the normal toric variety `v`, i.e. returns a list of affine toric varieties.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> affine_open_covering(p2)
3-element Vector{AffineNormalToricVariety}:
 Normal toric variety
 Normal toric variety
 Normal toric variety
```
"""
@attr Vector{AffineNormalToricVariety} function affine_open_covering(v::NormalToricVarietyType)
    charts = Vector{AffineNormalToricVariety}(undef, pm_object(v).N_MAXIMAL_CONES)
    for i in 1:pm_object(v).N_MAXIMAL_CONES
        charts[i] = affine_normal_toric_variety(cone(Polymake.fan.cone(pm_object(v), i-1)))
    end
    return charts
end
