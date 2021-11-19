############################
# Dimensions
############################


@doc Markdown.doc"""
    dim(v::AbstractNormalToricVariety)

Computes the dimension of the normal toric variety `v`.
"""
function dim(v::AbstractNormalToricVariety)
    return pm_object(v).FAN_DIM::Int
end
export dim


@doc Markdown.doc"""
    dim_of_torusfactor(v::AbstractNormalToricVariety)

Computes the dimension of the torus factor of the normal toric variety `v`.
"""
function dim_of_torusfactor(v::AbstractNormalToricVariety)

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

Computes the Euler characteristic of the normal toric variety `v`.
"""
function euler_characteristic(v::AbstractNormalToricVariety)
    f_vector = Vector{Int}(pm_object(v).F_VECTOR)
    return f_vector[dim(v)]
end
export euler_characteristic


############################
# Rings and ideals
############################


@doc Markdown.doc"""
    cox_ring(v::AbstractNormalToricVariety)

Computes the Cox ring of the normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> cox_ring(p2)
(Multivariate Polynomial Ring in x[1], x[2], x[3] over Rational Field graded by 
  x[1] -> [0 0 1]
  x[2] -> [0 0 1]
  x[3] -> [0 0 1], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3]])
```
"""
function cox_ring(v::AbstractNormalToricVariety)
    Qx, x = PolynomialRing(QQ, :x=>1:pm_object(v).N_RAYS)
    return grade(Qx,gens(class_group(v)))[1]
end
export cox_ring


@doc Markdown.doc"""
    stanley_reisner_ideal(v::AbstractNormalToricVariety)

Computes the Stanley-Reisner ideal of a normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> ngens(stanley_reisner_ideal(P2))
1
```
"""
function stanley_reisner_ideal(v::AbstractNormalToricVariety)
    collections = primitive_collections(fan(v))
    vars = Hecke.gens(cox_ring(v))
    SR_generators = [prod(vars[I]) for I in collections]
    return ideal(SR_generators)
end
export stanley_reisner_ideal


@doc Markdown.doc"""
    irrelevant_ideal(v::AbstractNormalToricVariety)

Computes the irrelevant ideal of a normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> length(irrelevant_ideal(p2).gens)
1
```
"""
function irrelevant_ideal(v::AbstractNormalToricVariety)
    # prepare maximal cone presentation
    max_cones = [findall(x->x!=0, l) for l in eachrow(pm_object(v).MAXIMAL_CONES)]
    n_ray = size(pm_object(v).RAYS, 1)
    maximal_cones = Vector{Int}[]
    for c in max_cones
        buffer = zeros(Int, n_ray)
        for k in c
            buffer[k] = 1
        end
        push!(maximal_cones, buffer)
    end
    # compute generators
    indeterminates = Hecke.gens(cox_ring(v))
    gens = MPolyElem_dec{fmpq, fmpq_mpoly}[]
    for i in 1:length(maximal_cones)
        monom = indeterminates[1]^(1 - maximal_cones[i][1])
        for j in 2:length(maximal_cones[i])
            monom = monom * indeterminates[j]^(1 - maximal_cones[i][j])
        end
        push!(gens, monom)
    end
    # return the ideal
    return ideal(gens)
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
A normal toric variety corresponding to a polyhedral fan in ambient dimension 3

julia> toric_ideal(antv)
ideal(-x[1]*x[2] + x[3]*x[4])
```
"""
function toric_ideal(antv::AffineNormalToricVariety)
    cone = Cone(pm_object(antv).WEIGHT_CONE)
    return toric_ideal(hilbert_basis(cone))
end
export toric_ideal
toric_ideal(ntv::NormalToricVariety) = toric_ideal(AffineNormalToricVariety(ntv))


############################
# Characters, Weil divisor and the class group
############################


@doc Markdown.doc"""
    character_lattice(v::AbstractNormalToricVariety)

Computes the character lattice of a normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> character_lattice(p2)
GrpAb: Z^2
```
"""
function character_lattice(v::AbstractNormalToricVariety)
    return abelian_group([0 for i in 1:pm_object(v).FAN_DIM])
end
export character_lattice


@doc Markdown.doc"""
    torusinvariant_divisor_group(v::AbstractNormalToricVariety)

Computes the torusinvariant divisor group of a normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> torusinvariant_divisor_group(p2)
GrpAb: Z^3
```
"""
function torusinvariant_divisor_group(v::AbstractNormalToricVariety)
    return abelian_group([0 for i in 1:pm_object(v).N_RAYS])
end
export torusinvariant_divisor_group


@doc Markdown.doc"""
    map_from_character_to_principal_divisors(v::AbstractNormalToricVariety)

Computes the map from the character lattice to the group of principal divisors of a normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> map_from_character_to_principal_divisors(p2)
Map with following data
Domain:
=======
Abelian group with structure: Z^2
Codomain:
=========
Abelian group with structure: Z^3
```
"""
function map_from_character_to_principal_divisors(v::AbstractNormalToricVariety)
    mat = transpose(Matrix{Int}(Polymake.common.primitive(pm_object(v).RAYS)))
    matrix = AbstractAlgebra.matrix(ZZ, mat)
    return hom(character_lattice(v), torusinvariant_divisor_group(v), matrix)
end
export map_from_character_to_principal_divisors


@doc Markdown.doc"""
    torusinvariant_prime_divisors(v::AbstractNormalToricVariety)

Computes the list of all torus invariant prime divisors in a normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> torusinvariant_prime_divisors(p2)
Free module of rank 3 over Integer Ring
```
"""
function torusinvariant_prime_divisors(v::AbstractNormalToricVariety)
    ti_divisors = torusinvariant_divisor_group(v)
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

Computes the class group of the normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> class_group(p2)
(General) abelian group with relation matrix
[1 0 -1 0 1 -1]
```
"""
function class_group(v::AbstractNormalToricVariety)
    return cokernel(map_from_character_to_principal_divisors(v))[1]
end
export class_group


@doc Markdown.doc"""
    map_from_weil_divisors_to_class_group(v::AbstractNormalToricVariety)

Computes the map from the group of Weil divisors to the class of group of a normal toric variety `v`.

# Examples
```jdoctest
julia> map_from_weil_divisors_to_class_group(p2)
Map with following data
Domain:
=======
Abelian group with structure: Z^3
Codomain:
=========
(General) abelian group with relation matrix
[0 0 0 0 0 0 0 0 0 1 0 -1 0 1 -1]
with structure of Abelian group with structure: Z
```
"""
function map_from_weil_divisors_to_class_group(v::AbstractNormalToricVariety)
    return cokernel(map_from_character_to_principal_divisors(v))[2]
end
export map_from_weil_divisors_to_class_group


@doc Markdown.doc"""
    map_from_cartier_divisor_group_to_torus_invariant_divisor_group(v::AbstractNormalToricVariety)

Computes the embedding of the group of Cartier divisors into the group of
torus-invariant Weil divisors of an abstract normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> map_from_cartier_divisor_group_to_torus_invariant_divisor_group(p2)
Map with following data
Domain:
=======
Abelian group with structure: Z^3
Codomain:
=========
Abelian group with structure: Z^3
```
"""
function map_from_cartier_divisor_group_to_torus_invariant_divisor_group(v::AbstractNormalToricVariety)
    # check input
    if hastorusfactor(v)
        throw(ArgumentError("Group of the torus-invariant Cartier divisors can only be computed if the variety has no torus factor."))
    end
    
    # identify rays and cones
    rays = Polymake.common.primitive(pm_object(v).RAYS)
    max_cones = incidence_matrix(maximal_cones(fan(v)))
    number_of_rays = size(rays)[1]
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
            map_for_scalar_products[(j-1)*rc+1:j*rc, col] = [fmpz(c) for c in rays[i,:]]
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
    map_to_weil_divisors = zero_matrix(ZZ, number_of_cones * rc, rank(torusinvariant_divisor_group(v)))
    for i in 1:number_of_rays
        map_to_weil_divisors[(cones_ray_is_part_of[i][1]-1)*rc+1:cones_ray_is_part_of[i][1]*rc, i] = [fmpz(-c) for c in rays[i,:]]
    end
    
    # compute the total map
    mapping_matrix = map_for_scalar_products * map_for_difference_of_elements
    source = abelian_group(zeros(Int, nrows(mapping_matrix)))
    target = abelian_group(zeros(Int, ncols(mapping_matrix)))
    total_map = hom(source, target, mapping_matrix)
    
    # identify the cartier_data_group and its embedding
    ker = kernel(total_map)
    cartier_data_group = snf(ker[1])[1]
    cartier_data_group_embedding = hom(cartier_data_group, source, ker[1].snf_map.map * ker[2].map)
    
    # construct the embedding of the Cartier divisor group into the group of torusinvariant divisors
    return cartier_divisor_group_embedding = image(hom(cartier_data_group, torusinvariant_divisor_group(v), cartier_data_group_embedding.map * map_to_weil_divisors))[2]
end
export map_from_cartier_divisor_group_to_torus_invariant_divisor_group

@doc Markdown.doc"""
    cartier_divisor_group(v::AbstractNormalToricVariety)

Computes the Cartier divisor group of an abstract normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> cartier_divisor_group(p2)
GrpAb: Z^3
```
"""
function cartier_divisor_group(v::AbstractNormalToricVariety)
    return domain(map_from_cartier_divisor_group_to_torus_invariant_divisor_group(v))
end
export cartier_divisor_group


@doc Markdown.doc"""
    picard_group(v::AbstractNormalToricVariety)

Computes the Picard group of an abstract normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> picard_group(p2)
GrpAb: Z
```
"""
function picard_group(v::AbstractNormalToricVariety)
    map1 = map_from_cartier_divisor_group_to_torus_invariant_divisor_group(v)
    map2 = map_from_weil_divisors_to_class_group(v)
    return snf(image(hom(domain(map1), codomain(map2), map1.map * map2.map))[1])[1]
end
export picard_group


############################
# Cones and fans
############################


"""
    nef_cone(v::NormalToricVariety)

Computes the nef cone of the normal toric variety `v`.

# Examples
```jldoctest
julia> pp = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> nef = nef_cone(pp)
A polyhedral cone in ambient dimension 1

julia> dim(nef)
1
```
"""
nef_cone(v::NormalToricVariety) = Cone(pm_object(v).NEF_CONE)
export nef_cone


"""
    mori_cone(v::NormalToricVariety)

Computes the mori cone of the normal toric variety `v`.

# Examples
```jldoctest
julia> pp = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> mori = mori_cone(pp)
A polyhedral cone in ambient dimension 1

julia> dim(mori)
1
```
"""
mori_cone(v::NormalToricVariety) = Cone(pm_object(v).MORI_CONE)
export mori_cone


@doc Markdown.doc"""
    fan(v::AbstractNormalToricVariety)

Computes the fan of an abstract normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> fan(p2)
A polyhedral fan in ambient dimension 2
```
"""
function fan(v::AbstractNormalToricVariety)
    return PolyhedralFan(pm_object(v))
end
export fan


@doc Markdown.doc"""
    cone(v::AffineNormalToricVariety)

Returns the cone of the affine normal toric variety `v`.

# Examples
```jdoctest
julia> cone(AffineNormalToricVariety(Oscar.positive_hull([1 1; -1 1])))
A polyhedral cone in ambient dimension 2
"""
function cone(v::AffineNormalToricVariety)
    return maximal_cones(fan(v))[1]
end
export cone


############################
# Affine covering
############################


@doc Markdown.doc"""
    affine_open_covering(v::AbstractNormalToricVariety)

Computes an affine open cover of the normal toric variety `v`, i.e. returns a list of affine toric varieties.

# Examples
```jldoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> affine_open_covering(p2)
3-element Vector{AffineNormalToricVariety}:
 A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
 A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
 A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function affine_open_covering(v::AbstractNormalToricVariety)
    charts = Vector{AffineNormalToricVariety}(undef, pm_object(v).N_MAXIMAL_CONES)
    for i in 1:pm_object(v).N_MAXIMAL_CONES
        charts[i] = AffineNormalToricVariety(Cone(Polymake.fan.cone(pm_object(v), i-1)))
    end
    return charts
end
export affine_open_covering
