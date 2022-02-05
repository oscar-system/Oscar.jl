############################
# Dimensions
############################


@doc Markdown.doc"""
    dim(v::AbstractNormalToricVariety)

Return the dimension of the normal toric variety `v`.
"""
function dim(v::AbstractNormalToricVariety)
    return get_attribute!(v, :dim) do
        return pm_object(v).FAN_DIM
    end::Int
end
export dim


@doc Markdown.doc"""
    dim_of_torusfactor(v::AbstractNormalToricVariety)

Return the dimension of the torus factor of the normal toric variety `v`.
"""
function dim_of_torusfactor(v::AbstractNormalToricVariety)
    return get_attribute!(v, :dim_of_torusfactor) do
        if hastorusfactor(v) == false
            return 0
        end
        dimension_of_fan = pm_object(v).FAN_DIM::Int
        ambient_dimension = pm_object(v).FAN_AMBIENT_DIM::Int
        return ambient_dimension - dimension_of_fan
    end::Int
end
export dim_of_torusfactor


@doc Markdown.doc"""
    euler_characteristic(v::AbstractNormalToricVariety)

Return the Euler characteristic of the normal toric variety `v`.
"""
function euler_characteristic(v::AbstractNormalToricVariety)
    return get_attribute!(v, :euler_characteristic) do
        f_vector = Vector{Int}(pm_object(v).F_VECTOR)
        return f_vector[dim(v)]
    end
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
    if has_attribute(v, :cox_ring)
        error("Cox ring already constructed. Coefficient ring must not be changed.")
    end
    if has_attribute(v, :toric_ideal)
        error("Toric ideal already constructed. Coefficient ring must not be changed.")
    end
    set_attribute!(v, :coefficient_ring, coefficient_ring)
end
export set_coefficient_ring


@doc Markdown.doc"""
    coefficient_ring(v::AbstractNormalToricVariety)

This method returns the coefficient_ring of the normal toric variety `v`.
An error is triggered if it is not yet set.
"""
function coefficient_ring(v::AbstractNormalToricVariety)
    if !has_attribute(v, :coefficient_ring)
        error("Coefficient ring not yet set.")
    end
    return get_attribute(v, :coefficient_ring)
end
export coefficient_ring


@doc Markdown.doc"""
    set_coordinate_names(v::AbstractNormalToricVariety, coordinate_names::Vector{String})

Allows to set the names of the homogeneous coordinates. If
the Cox ring of the variety has already been computed, we do
not allow changes of the coordinate names. In this case, an error
is triggered.
"""
function set_coordinate_names(v::AbstractNormalToricVariety, coordinate_names::Vector{String})
    if has_attribute(v, :cox_ring)
        error("Cox ring already constructed. Coordinate names cannot be changed anymore.")
    end
    if length(coordinate_names) != nrays(fan(v))
        throw(ArgumentError("The provided list of coordinate names must match the number of rays in the fan."))
    end
    set_attribute!(v, :coordinate_names, coordinate_names)
end
export set_coordinate_names


@doc Markdown.doc"""
    coordinate_names(v::AbstractNormalToricVariety)

This method returns the names of the homogeneous coordinates of 
the normal toric variety `v`. If they are not yet set an error is returned.
"""
function coordinate_names(v::AbstractNormalToricVariety)
    if !has_attribute(v, :coordinate_names)
        error("Coordinate names not yet set.")
    end
    return get_attribute(v, :coordinate_names)
end
export coordinate_names


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
    return get_attribute!(v, :cox_ring) do
        
        # are coordinates names set? If not, set default values
        if !has_attribute(v, :coordinate_names)
            set_attribute!(v, :coordinate_names, ["x" * string(i) for i in 1:rank(torusinvariant_divisor_group(v))])
        end
        
        # is the coefficient_ring set? If not, set default value
        if !has_attribute(v, :coefficient_ring)
            set_attribute!(v, :coefficient_ring, QQ)
        end
        
        # construct the cox ring
        S, _ = PolynomialRing(coefficient_ring(v), coordinate_names(v))
        weights = [map_from_weil_divisors_to_class_group(v)(x) for x in gens(torusinvariant_divisor_group(v))]        
        return grade(S,weights)[1]
        
    end
end
export cox_ring


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
function stanley_reisner_ideal(v::AbstractNormalToricVariety)
    return get_attribute!(v, :stanley_reisner_ideal) do
        if has_attribute(v, :polyhedron)
            # compute via simplicial complex
            I = pm_object(get_attribute(v, :polyhedron)).FACETS_THRU_VERTICES
            K = SimplicialComplex(I)
            return stanley_reisner_ideal(cox_ring(v), K)
        else
            # compute via primitive collections
            collections = primitive_collections(fan(v))
            vars = Hecke.gens(cox_ring(v))
            SR_generators = [prod(vars[I]) for I in collections]
            return ideal(SR_generators)
        end
    end
end
export stanley_reisner_ideal


@doc Markdown.doc"""
    irrelevant_ideal(v::AbstractNormalToricVariety)

Return the irrelevant ideal of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> length(irrelevant_ideal(p2).gens)
3
```
"""
function irrelevant_ideal(v::AbstractNormalToricVariety)
    return get_attribute!(v, :irrelevant_ideal) do
        # prepare maximal cones
        max_cones = [findall(x->x!=0, l) for l in eachrow(pm_object(v).MAXIMAL_CONES)]
        maximal_cones = Vector{Int}[]
        for c in max_cones
            buffer = zeros(Int, nrays(fan(v)))
            for k in c
                buffer[k] = 1
            end
            push!(maximal_cones, buffer)
        end
        
        # compute generators
        indeterminates = Hecke.gens(cox_ring(v))
        gens = typeof(indeterminates[1])[]
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
    return get_attribute!(antv, :toric_ideal) do
        # is the coefficient_ring set? If not, set default value
        if !has_attribute(antv, :coefficient_ring)
            set_attribute!(antv, :coefficient_ring, QQ)
        end
        
        # construct the toric ideal
        cone = Cone(pm_object(antv).WEIGHT_CONE)
        return toric_ideal(hilbert_basis(cone), coefficient_ring(antv))
    end
end
export toric_ideal

function toric_ideal(ntv::NormalToricVariety)
    isaffine(ntv) || error("Cannot construct affine toric variety from non-affine input")    
    return get_attribute!(() -> toric_ideal(AffineNormalToricVariety(ntv)), ntv, :toric_ideal)
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
function character_lattice(v::AbstractNormalToricVariety)
    return get_attribute!(v, :character_lattice) do
        return free_abelian_group(ambient_dim(fan(v)))
    end
end
export character_lattice


@doc Markdown.doc"""
    torusinvariant_divisor_group(v::AbstractNormalToricVariety)

Return the torusinvariant divisor group of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> torusinvariant_divisor_group(p2)
GrpAb: Z^3
```
"""
function torusinvariant_divisor_group(v::AbstractNormalToricVariety)
    return get_attribute!(v, :torusinvariant_divisor_group) do
        return free_abelian_group(nrays(fan(v)))
    end
end
export torusinvariant_divisor_group


@doc Markdown.doc"""
    map_from_character_to_principal_divisors(v::AbstractNormalToricVariety)

Return the map from the character lattice to the group of principal divisors of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

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
    return get_attribute!(v, :map_from_character_to_principal_divisors) do
        mat = transpose(matrix(ZZ, rays(fan(v))))
        return hom(character_lattice(v), torusinvariant_divisor_group(v), mat)
    end
end
export map_from_character_to_principal_divisors


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
function torusinvariant_prime_divisors(v::AbstractNormalToricVariety)
    return get_attribute!(v, :torusinvariant_prime_divisors) do
        ti_divisors = torusinvariant_divisor_group(v)
        prime_divisors = ToricDivisor[]
        for i in 1:rank(ti_divisors)
            coeffs = zeros(Int, rank(ti_divisors))
            coeffs[i] = 1
            push!(prime_divisors, ToricDivisor(v,coeffs))
        end
        return prime_divisors
    end
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
function class_group(v::AbstractNormalToricVariety)
    return get_attribute!(v, :class_group) do
        return codomain(map_from_weil_divisors_to_class_group(v))
    end
end
export class_group


@doc Markdown.doc"""
    map_from_weil_divisors_to_class_group(v::AbstractNormalToricVariety)

Return the map from the group of Weil divisors to the class of group of a normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> map_from_weil_divisors_to_class_group(p2)
Map with following data
Domain:
=======
Abelian group with structure: Z^3
Codomain:
=========
Abelian group with structure: Z
```
"""
function map_from_weil_divisors_to_class_group(v::AbstractNormalToricVariety)
    return get_attribute!(v, :map_from_weil_divisors_to_class_group) do
        map1 = cokernel(map_from_character_to_principal_divisors(v))[2]
        map2 = inv(snf(codomain(map1))[2])
        return map1*map2
    end        
end
export map_from_weil_divisors_to_class_group


@doc Markdown.doc"""
    map_from_cartier_divisor_group_to_torus_invariant_divisor_group(v::AbstractNormalToricVariety)

Return the embedding of the group of Cartier divisors into the group of
torus-invariant Weil divisors of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> map_from_cartier_divisor_group_to_torus_invariant_divisor_group(p2)
Identity map with

Domain:
=======
GrpAb: Z^3
```
"""
function map_from_cartier_divisor_group_to_torus_invariant_divisor_group(v::AbstractNormalToricVariety)

    return get_attribute!(v, :map_from_cartier_divisor_group_to_torus_invariant_divisor_group) do
        # check input
        if hastorusfactor(v)
            throw(ArgumentError("Group of the torus-invariant Cartier divisors can only be computed if the variety has no torus factor."))
        end
        
        # identify fan_rays and cones
        fan_rays = transpose(matrix(ZZ, rays(fan(v))))
        max_cones = ray_indices(maximal_cones(fan(v)))
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
        map_to_weil_divisors = zero_matrix(ZZ, number_of_cones * rc, rank(torusinvariant_divisor_group(v)))
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
        embedding = snf(ker[1])[2] * ker[2] * hom(codomain(ker[2]), torusinvariant_divisor_group(v), map_to_weil_divisors)
        
        # return the image of this embedding
        return image(embedding)[2]
    end
end
export map_from_cartier_divisor_group_to_torus_invariant_divisor_group

@doc Markdown.doc"""
    cartier_divisor_group(v::AbstractNormalToricVariety)

Return the Cartier divisor group of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> cartier_divisor_group(p2)
GrpAb: Z^3
```
"""
function cartier_divisor_group(v::AbstractNormalToricVariety)
    return get_attribute!(v, :cartier_divisor_group) do
        return domain(map_from_cartier_divisor_group_to_torus_invariant_divisor_group(v))
    end
end
export cartier_divisor_group


@doc Markdown.doc"""
    map_from_cartier_divisor_group_to_picard_group(v::AbstractNormalToricVariety)

Return the map from the Cartier divisors to the Picard group 
of an abstract normal toric variety `v`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> map_from_cartier_divisor_group_to_picard_group(p2)
Map with following data
Domain:
=======
Abelian group with structure: Z^3
Codomain:
=========
Abelian group with structure: Z
```
"""
function map_from_cartier_divisor_group_to_picard_group(v::AbstractNormalToricVariety)
    return get_attribute!(v, :map_from_cartier_divisor_group_to_picard_group) do
        # check input
        if hastorusfactor(v)
            throw(ArgumentError("Group of the torus-invariant Cartier divisors can only be computed if the variety has no torus factor."))
        end
        
        # compute mapping
        map1 = map_from_cartier_divisor_group_to_torus_invariant_divisor_group(v)
        map2 = map_from_weil_divisors_to_class_group(v)
        return restrict_codomain(map1*map2)
        
    end
end
export map_from_cartier_divisor_group_to_picard_group


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
function picard_group(v::AbstractNormalToricVariety)
    return get_attribute!(v, :picard_group) do
        return codomain(map_from_cartier_divisor_group_to_picard_group(v))
    end
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
julia> pp = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> nef = nef_cone(pp)
A polyhedral cone in ambient dimension 1

julia> dim(nef)
1
```
"""
function nef_cone(v::NormalToricVariety)
    return get_attribute!(v, :nef_cone) do
        return Cone(pm_object(v).NEF_CONE)
    end
end
export nef_cone


"""
    mori_cone(v::NormalToricVariety)

Return the mori cone of the normal toric variety `v`.

# Examples
```jldoctest
julia> pp = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> mori = mori_cone(pp)
A polyhedral cone in ambient dimension 1

julia> dim(mori)
1
```
"""
function mori_cone(v::NormalToricVariety)
    return get_attribute!(v, :mori_cone) do
        return Cone(pm_object(v).MORI_CONE)
    end
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
function fan(v::AbstractNormalToricVariety)
    return get_attribute!(v, :fan) do
        return PolyhedralFan(pm_object(v))
    end
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
function cone(v::AffineNormalToricVariety)
    return get_attribute!(v, :cone) do
        return maximal_cones(fan(v))[1]
    end
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
function affine_open_covering(v::AbstractNormalToricVariety)
    return get_attribute!(v, :affine_open_covering) do
        charts = Vector{AffineNormalToricVariety}(undef, pm_object(v).N_MAXIMAL_CONES)
        for i in 1:pm_object(v).N_MAXIMAL_CONES
            charts[i] = AffineNormalToricVariety(Cone(Polymake.fan.cone(pm_object(v), i-1)))
        end
        return charts
    end    
end
export affine_open_covering
