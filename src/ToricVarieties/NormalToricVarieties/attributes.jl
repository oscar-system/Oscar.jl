############################
# Dimensions
############################


@doc Markdown.doc"""
    dim(v::AbstractNormalToricVariety)

Return the dimension of the normal toric variety `v`.
"""
@attr Int function dim(v::AbstractNormalToricVariety)
    return pm_object(v).FAN_DIM
end
export dim


@doc Markdown.doc"""
    dim_of_torusfactor(v::AbstractNormalToricVariety)

Return the dimension of the torus factor of the normal toric variety `v`.
"""
@attr Int function dim_of_torusfactor(v::AbstractNormalToricVariety)
    if hastorusfactor(v) == false
        return 0
    end
    dimension_of_fan = pm_object(v).FAN_DIM::Int
    ambient_dimension = pm_object(v).FAN_AMBIENT_DIM::Int
    return ambient_dimension - dimension_of_fan
end
export dim_of_torusfactor


@doc Markdown.doc"""
    euler_characteristic(v::AbstractNormalToricVariety)

Return the Euler characteristic of the normal toric variety `v`.
"""
@attr Int function euler_characteristic(v::AbstractNormalToricVariety)
    f_vector = Vector{Int}(pm_object(v).F_VECTOR)
    return f_vector[dim(v)]
end
export euler_characteristic


############################
# Rings and ideals
############################


@doc Markdown.doc"""
    set_coefficient_ring(v::AbstractNormalToricVariety, coefficient_ring::AbstractAlgebra.Ring)

Allows to set the coefficient_ring. If the Cox ring of the variety has
already been computed, we do not allow this to be changed.
In this case an error is triggered.
"""
function set_coefficient_ring(v::AbstractNormalToricVariety, coefficient_ring::AbstractAlgebra.Ring)
    set_attribute!(v, :coefficient_ring, coefficient_ring)
end
export set_coefficient_ring


@doc Markdown.doc"""
    coefficient_ring(v::AbstractNormalToricVariety)

This method returns the coefficient_ring of the normal toric variety `v`.
The default is the ring `QQ`.
"""
coefficient_ring(v::AbstractNormalToricVariety) = get_attribute!(v, :coefficient_ring, QQ)


@doc Markdown.doc"""
    set_coordinate_names(v::AbstractNormalToricVariety, coordinate_names::Vector{String})

Allows to set the names of the homogeneous coordinates.
"""
function set_coordinate_names(v::AbstractNormalToricVariety, coordinate_names::Vector{String})
    if length(coordinate_names) != nrays(v)
        throw(ArgumentError("The provided list of coordinate names must match the number of rays in the fan."))
    end
    set_attribute!(v, :coordinate_names, coordinate_names)
end
export set_coordinate_names


@doc Markdown.doc"""
    coordinate_names(v::AbstractNormalToricVariety)

This method returns the names of the homogeneous coordinates of 
the normal toric variety `v`. The default is `x1,...,xn`.
"""
coordinate_names(v::AbstractNormalToricVariety) = get_attribute!(v, :coordinate_names, ["x$(i)" for i in 1:rank(torusinvariant_weil_divisor_group(v))])


function _cox_ring_weights(v::AbstractNormalToricVariety)
    return get_attribute(v, :cox_ring_weights) do
        return [map_from_torusinvariant_weil_divisor_group_to_class_group(v)(x) for x in gens(torusinvariant_weil_divisor_group(v))]
    end
end

@doc Markdown.doc"""
    cox_ring(v::AbstractNormalToricVariety)

Computes the Cox ring of the normal toric variety `v`.
Note that [CLS11](@cite) refers to this ring as the "total coordinate ring".

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> set_coordinate_names(p2, ["y1", "y2", "y3"])

julia> set_coefficient_ring(p2, ZZ)

julia> cox_ring(p2)
Multivariate Polynomial Ring in y1, y2, y3 over Integer Ring graded by 
  y1 -> [1]
  y2 -> [1]
  y3 -> [1]
```
"""
function cox_ring(v::AbstractNormalToricVariety)
    S, _ = PolynomialRing(coefficient_ring(v), coordinate_names(v), cached=false)
    return cox_ring(S, v)
end


@doc Markdown.doc"""
    cox_ring(R::MPolyRing, v::AbstractNormalToricVariety)

Computes the Cox ring of the normal toric variety `v`, in this case by adding
the Cox grading to the given ring `R`.
Note that [CLS11](@cite) refers to this ring as the "total coordinate ring".

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R,_ = PolynomialRing(QQ, 3);

julia> cox_ring(R, p2)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field graded by 
  x1 -> [1]
  x2 -> [1]
  x3 -> [1]
```
"""
function cox_ring(R::MPolyRing, v::AbstractNormalToricVariety)
    weights = _cox_ring_weights(v)
    length(weights) == nvars(R) || throw(ArgumentError("Wrong number of variables"))
    return grade(R, weights)[1]
end
export cox_ring


@attr Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric} function _minimal_nonfaces(v::AbstractNormalToricVariety)
    I = ray_indices(maximal_cones(v))
    K = SimplicialComplex(I)
    return minimal_nonfaces(IncidenceMatrix, K)
end

@doc Markdown.doc"""
    stanley_reisner_ideal(R::MPolyRing, v::AbstractNormalToricVariety)

Return the Stanley-Reisner ideal of a normal toric variety `v` as an ideal of
`R`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R,_ = PolynomialRing(QQ, 3);

julia> ngens(stanley_reisner_ideal(R, p2))
1
```
"""
function stanley_reisner_ideal(R::MPolyRing, v::AbstractNormalToricVariety)
    n = nrays(v)
    n == nvars(R) || throw(ArgumentError("Wrong number of variables"))
    mnf = _minimal_nonfaces(v)
    return ideal([ R([1], [Vector{Int}(mnf[i,:])]) for i in 1:Polymake.nrows(mnf) ])
end

@doc Markdown.doc"""
    stanley_reisner_ideal(v::AbstractNormalToricVariety)

Return the Stanley-Reisner ideal of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> ngens(stanley_reisner_ideal(p2))
1
```
"""
stanley_reisner_ideal(v::AbstractNormalToricVariety) = stanley_reisner_ideal(cox_ring(v), v)
export stanley_reisner_ideal


@attr Vector{Vector{Int}} function _irrelevant_ideal_monomials(v::AbstractNormalToricVariety)
    mc = ray_indices(maximal_cones(v))
    result = Vector{Vector{Int}}()
    onesv = ones(Int, Polymake.ncols(mc))
    for i in 1:Polymake.nrows(mc)
        push!(result, onesv - Vector{Int}(mc[i,:]))
    end
    return result
end


@doc Markdown.doc"""
    irrelevant_ideal(v::AbstractNormalToricVariety)

Return the irrelevant ideal of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> length(gens(irrelevant_ideal(p2)))
3
```
"""
function irrelevant_ideal(v::AbstractNormalToricVariety)
    R = cox_ring(v)
    return irrelevant_ideal(R, v)
end

@doc Markdown.doc"""
    irrelevant_ideal(R::MPolyRing, v::AbstractNormalToricVariety)

Return the irrelevant ideal of a normal toric variety `v` as an ideal in `R`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> R,_ = PolynomialRing(QQ, 3);

julia> length(gens(irrelevant_ideal(R, p2)))
3
```
"""
function irrelevant_ideal(R::MPolyRing, v::AbstractNormalToricVariety)
    monoms = _irrelevant_ideal_monomials(v)
    nvars(R) == nrays(v) || throw(ArgumentError("Wrong number of variables in polynomial ring."))
    return ideal([R([1], [x]) for x in monoms])
end
export irrelevant_ideal



@doc Markdown.doc"""
    toric_ideal(antv::AffineNormalToricVariety)

Return the toric ideal defining the affine normal toric variety.

# Examples
Take the cone over the square at height one. The resulting toric variety has
one defining equation. In projective space this corresponds to
$\mathbb{P}^1\times\mathbb{P}^1$. Note that this cone is self-dual, the toric
ideal comes from the dual cone.
```jldoctest
julia> C = positive_hull([1 0 0; 1 1 0; 1 0 1; 1 1 1])
A polyhedral cone in ambient dimension 3

julia> antv = AffineNormalToricVariety(C)
A normal, affine toric variety

julia> toric_ideal(antv)
ideal(-x1*x2 + x3*x4)
```
"""
function toric_ideal(antv::AffineNormalToricVariety)
    cone = Cone(pm_object(antv).WEIGHT_CONE)
    n = length(hilbert_basis(cone))
    R,_ = PolynomialRing(coefficient_ring(antv), n, cached=false)
    return toric_ideal(R, antv)
end

@doc Markdown.doc"""
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
A polyhedral cone in ambient dimension 3

julia> antv = AffineNormalToricVariety(C)
A normal, affine toric variety

julia> R,_ = PolynomialRing(QQ, 4);

julia> toric_ideal(R, antv)
ideal(-x1*x2 + x3*x4)
```
"""
function toric_ideal(R::MPolyRing, antv::AffineNormalToricVariety)
    cone = Cone(pm_object(antv).WEIGHT_CONE)
    gens = pm_object(cone).CONE_TORIC_IDEAL.BINOMIAL_GENERATORS
    return binomial_exponents_to_ideal(R, gens)
end

function toric_ideal(ntv::NormalToricVariety)
    isaffine(ntv) || error("Cannot construct affine toric variety from non-affine input")
    return toric_ideal(AffineNormalToricVariety(ntv))
end
export toric_ideal


############################
# Characters, Weil divisor and the class group
############################


@doc Markdown.doc"""
    character_lattice(v::AbstractNormalToricVariety)

Return the character lattice of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> character_lattice(p2)
GrpAb: Z^2
```
"""
@attr GrpAbFinGen function character_lattice(v::AbstractNormalToricVariety)
    return free_abelian_group(ambient_dim(v))
end
export character_lattice


@doc Markdown.doc"""
    torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)

Return the torusinvariant divisor group of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> torusinvariant_weil_divisor_group(p2)
GrpAb: Z^3
```
"""
@attr GrpAbFinGen function torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)
    return free_abelian_group(nrays(v))
end
export torusinvariant_weil_divisor_group


@doc Markdown.doc"""
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
export map_from_character_lattice_to_torusinvariant_weil_divisor_group


@doc Markdown.doc"""
    torusinvariant_prime_divisors(v::AbstractNormalToricVariety)

Return the list of all torus invariant prime divisors in a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> torusinvariant_prime_divisors(p2)
3-element Vector{ToricDivisor}:
 A torus-invariant, prime divisor on a normal toric variety
 A torus-invariant, prime divisor on a normal toric variety
 A torus-invariant, prime divisor on a normal toric variety
```
"""
@attr Vector{ToricDivisor} function torusinvariant_prime_divisors(v::AbstractNormalToricVariety)
    ti_divisors = torusinvariant_weil_divisor_group(v)
    prime_divisors = ToricDivisor[]
    for i in 1:rank(ti_divisors)
        coeffs = zeros(Int, rank(ti_divisors))
        coeffs[i] = 1
        push!(prime_divisors, ToricDivisor(v,coeffs))
    end
    return prime_divisors
end
export torusinvariant_prime_divisors


@doc Markdown.doc"""
    class_group(v::AbstractNormalToricVariety)

Return the class group of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> class_group(p2)
GrpAb: Z
```
"""
@attr GrpAbFinGen function class_group(v::AbstractNormalToricVariety)
    return codomain(map_from_torusinvariant_weil_divisor_group_to_class_group(v))
end
export class_group


@doc Markdown.doc"""
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
export map_from_torusinvariant_weil_divisor_group_to_class_group


@doc Markdown.doc"""
    map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)

Return the embedding of the group of Cartier divisors into the group of
torus-invariant Weil divisors of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(p2)
Identity map with

Domain:
=======
GrpAb: Z^3
```
"""
@attr Map{GrpAbFinGen, GrpAbFinGen} function map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v::AbstractNormalToricVariety)
    # check input
    if hastorusfactor(v)
        throw(ArgumentError("Group of the torus-invariant Cartier divisors can only be computed if the variety has no torus factor."))
    end

    # identify fan_rays and cones
    fan_rays = transpose(matrix(ZZ, rays(v)))
    max_cones = ray_indices(maximal_cones(v))
    number_of_rays = ncols(fan_rays)
    number_of_cones = size(max_cones)[1]

    # compute quantities needed to construct the matrices
    rc = rank(character_lattice(v))
    number_ray_is_part_of_max_cones = [length(max_cones[:,k].s) for k in 1:number_of_rays]
    s = sum(number_ray_is_part_of_max_cones)
    cones_ray_is_part_of = [filter(x -> max_cones[x,r], 1:number_of_cones) for r in 1:number_of_rays]

    # compute the matrix for the scalar products
    map_for_scalar_products = zero_matrix(ZZ, number_of_cones * rc, s)
    col = 1
    for i in 1:number_of_rays
        for j in cones_ray_is_part_of[i]
            map_for_scalar_products[(j-1)*rc+1:j*rc, col] = [fmpz(c) for c in fan_rays[:,i]]
            col += 1
        end
    end
    
    # compute the matrix for differences
    map_for_difference_of_elements = zero_matrix(ZZ, s, s-number_of_rays)
    row = 1
    col = 1
    for i in 1:number_of_rays
        ncol = number_ray_is_part_of_max_cones[i]-1
        map_for_difference_of_elements[row, col:(col+ncol-1)] = fill(fmpz(1), ncol)
        map_for_difference_of_elements[(row+1):(row+ncol), col:(col+ncol-1)] = -identity_matrix(ZZ, ncol)
        row += ncol + 1
        col += ncol
    end
    
    # compute the matrix for mapping to torusinvariant Weil divisors
    map_to_weil_divisors = zero_matrix(ZZ, number_of_cones * rc, rank(torusinvariant_weil_divisor_group(v)))
    for i in 1:number_of_rays
        map_to_weil_divisors[(cones_ray_is_part_of[i][1]-1)*rc+1:cones_ray_is_part_of[i][1]*rc, i] = [fmpz(-c) for c in fan_rays[:,i]]
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
export map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group

@doc Markdown.doc"""
    torusinvariant_cartier_divisor_group(v::AbstractNormalToricVariety)

Return the Cartier divisor group of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> torusinvariant_cartier_divisor_group(p2)
GrpAb: Z^3
```
"""
@attr GrpAbFinGen function torusinvariant_cartier_divisor_group(v::AbstractNormalToricVariety)
    return domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v))
end
export torusinvariant_cartier_divisor_group


@doc Markdown.doc"""
    map_from_torusinvariant_cartier_divisor_group_to_picard_group(v::AbstractNormalToricVariety)

Return the map from the Cartier divisors to the Picard group
of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

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
    if hastorusfactor(v)
        throw(ArgumentError("Group of the torus-invariant Cartier divisors can only be computed if the variety has no torus factor."))
    end
    
    # compute mapping
    map1 = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v)
    map2 = map_from_torusinvariant_weil_divisor_group_to_class_group(v)
    return restrict_codomain(map1*map2)
end
export map_from_torusinvariant_cartier_divisor_group_to_picard_group


@doc Markdown.doc"""
    picard_group(v::AbstractNormalToricVariety)

Return the Picard group of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> picard_group(p2)
GrpAb: Z
```
"""
@attr GrpAbFinGen function picard_group(v::AbstractNormalToricVariety)
    return codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(v))
end
export picard_group


############################
# Cones and fans
############################


"""
    nef_cone(v::NormalToricVariety)

Return the nef cone of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> nef = nef_cone(p2)
A polyhedral cone in ambient dimension 1

julia> dim(nef)
1
```
"""
@attr Cone function nef_cone(v::NormalToricVariety)
    return Cone(pm_object(v).NEF_CONE)
end
export nef_cone


"""
    mori_cone(v::NormalToricVariety)

Return the mori cone of the normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> mori = mori_cone(p2)
A polyhedral cone in ambient dimension 1

julia> dim(mori)
1
```
"""
@attr Cone function mori_cone(v::NormalToricVariety)
    return Cone(pm_object(v).MORI_CONE)
end
export mori_cone


@doc Markdown.doc"""
    fan(v::AbstractNormalToricVariety)

Return the fan of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> fan(p2)
A polyhedral fan in ambient dimension 2
```
"""
@attr PolyhedralFan function fan(v::AbstractNormalToricVariety)
    return PolyhedralFan(pm_object(v))
end
export fan


@doc Markdown.doc"""
    cone(v::AffineNormalToricVariety)

Return the cone of the affine normal toric variety `v`.

# Examples
```jldoctest
julia> cone(AffineNormalToricVariety(Oscar.positive_hull([1 1; -1 1])))
A polyhedral cone in ambient dimension 2
```
"""
@attr Cone function cone(v::AffineNormalToricVariety)
    return maximal_cones(v)[1]
end
export cone


############################
# Affine covering
############################


@doc Markdown.doc"""
    affine_open_covering(v::AbstractNormalToricVariety)

Compute an affine open cover of the normal toric variety `v`, i.e. returns a list of affine toric varieties.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> affine_open_covering(p2)
3-element Vector{AffineNormalToricVariety}:
 A normal, affine toric variety
 A normal, affine toric variety
 A normal, affine toric variety
```
"""
@attr Vector{AffineNormalToricVariety} function affine_open_covering(v::AbstractNormalToricVariety)
    charts = Vector{AffineNormalToricVariety}(undef, pm_object(v).N_MAXIMAL_CONES)
    for i in 1:pm_object(v).N_MAXIMAL_CONES
        charts[i] = AffineNormalToricVariety(Cone(Polymake.fan.cone(pm_object(v), i-1)))
    end
    return charts
end
export affine_open_covering
