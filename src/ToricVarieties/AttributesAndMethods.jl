######################
# 1: Attributes of ToricVarieties
######################

@doc Markdown.doc"""
    dim(v::AbstractNormalToricVariety)

Computes the dimension of the normal toric variety `v`.
"""
function dim(v::AbstractNormalToricVariety)
    return pm_ntv(v).FAN_DIM::Int
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
    
    dimension_of_fan = pm_ntv(v).FAN_DIM::Int
    ambient_dimension = pm_ntv(v).FAN_AMBIENT_DIM::Int
    return ambient_dimension - dimension_of_fan
end
export dim_of_torusfactor


@doc Markdown.doc"""
    euler_characteristic(v::AbstractNormalToricVariety)

Computes the Euler characteristic of the normal toric variety `v`.
"""
function euler_characteristic(v::AbstractNormalToricVariety)
    f_vector = Vector{Int}(pm_ntv(v).F_VECTOR)
    return f_vector[dim(v)]
end
export euler_characteristic


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
nef_cone(v::NormalToricVariety) = Cone(pm_ntv(v).NEF_CONE)
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
mori_cone(v::NormalToricVariety) = Cone(pm_ntv(v).MORI_CONE)
export mori_cone


@doc Markdown.doc"""
    affine_open_covering(v::NormalToricVariety)

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
function affine_open_covering(v::NormalToricVariety)
    charts = Vector{AffineNormalToricVariety}(undef, pm_ntv(v).N_MAXIMAL_CONES)
    for i in 1:pm_ntv(v).N_MAXIMAL_CONES
        charts[i] = AffineNormalToricVariety(Cone(Polymake.fan.cone(pm_ntv(v), i-1)))
    end
    return charts
end


@doc Markdown.doc"""
    affine_open_covering(v::AffineNormalToricVariety)

Computes an affine open cover of the affine normal toric variety `v`, i.e. returns this very variety.

# Examples
```jdoctest
julia> C = Oscar.positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> affine_open_covering(antv)
1-element Vector{AffineNormalToricVariety}:
 A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
```
"""
function affine_open_covering(v::AffineNormalToricVariety)
    return [v]
end
export affine_open_covering


@doc Markdown.doc"""
    fan_of_variety(v::NormalToricVariety)

Computes the fan of a normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> fan_of_variety(p2)
A polyhedral fan in ambient dimension 2
```
"""
function fan_of_variety(v::NormalToricVariety)
    return PolyhedralFan(pm_ntv(v))
end
export fan_of_variety


@doc Markdown.doc"""
    fan(v::NormalToricVariety)

A convenience method for fan_of_variety of a normal toric variety `v`.
"""
function fan(v::NormalToricVariety)
    return fan_of_variety(v)
end
export fan


@doc Markdown.doc"""
    torusinvariant_divisor_group(v::NormalToricVariety)

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
    return abelian_group([0 for i in 1:pm_ntv(v).N_RAYS])
end
export torusinvariant_divisor_group


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
    return abelian_group([0 for i in 1:pm_ntv(v).FAN_DIM])
end
export character_lattice


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
    matrix = Matrix{Int}(Oscar.Polymake.common.primitive(pm_ntv(v).RAYS))
    abstract_matrix = AbstractAlgebra.matrix(ZZ, size(matrix,2), size(matrix,1), vec(matrix))
    return hom(character_lattice(v), torusinvariant_divisor_group(v), abstract_matrix)
end
export map_from_character_to_principal_divisors


@doc Markdown.doc"""
    class_group(v::AbstractNormalToricVariety)

Computes the class group of the normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> class_group(p2)
(General) abelian group with relation matrix
[1 0 -1; 0 1 -1]
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
[0 0 0; 0 0 0; 0 0 0; 1 0 -1; 0 1 -1]
with structure of Abelian group with structure: Z
```
"""
function map_from_weil_divisors_to_class_group(v::AbstractNormalToricVariety)
    return cokernel(map_from_character_to_principal_divisors(v))[2]
end
export map_from_weil_divisors_to_class_group


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
    Qx, x = PolynomialRing(QQ, :x=>1:pm_ntv(v).N_RAYS)
    return grade(Qx,gens(class_group(v)))[1]
end
export cox_ring


@doc Markdown.doc"""
    list_of_variables_of_cox_ring(v::NormalToricVariety)

Lists the names of the indeterminates of the Cox ring of a normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> list_of_variables_of_cox_ring(p2)
3-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x1
 x2
 x3
```
"""
function list_of_variables_of_cox_ring(v::AbstractNormalToricVariety)
    return gens(cox_ring(v))
end
export list_of_variables_of_cox_ring


@doc Markdown.doc"""
    stanley_reisner_ideal(v::NormalToricVariety)
Computes the Stanley-Reisner ideal of a normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> length(stanley_reisner_ideal(P2).gens)
1
```
"""
function stanley_reisner_ideal(v::AbstractNormalToricVariety)
    collections = primitive_collections(fan_of_variety(v))
    SR_generators = []
    vars = list_of_variables_of_cox_ring(v)
    for I in collections
        buffer = vars[I[1]]
        for k in 2 : length(I)
           buffer = buffer * vars[I[k]]
        end
        push!(SR_generators, buffer)
    end
    return ideal(SR_generators)
end
export stanley_reisner_ideal


@doc Markdown.doc"""
    irrelevant_ideal(v::NormalToricVariety)

Computes the irrelevant ideal of a normal toric variety `v`.

# Examples
```jdoctest
julia> p2 = toric_projective_space(2)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> length(irrelevant_ideal(p2).gens)
1
```
"""
function irrelevant_ideal(v::NormalToricVariety)
    # prepare maximal cone presentation
    max_cones = [findall(x->x!=0, l) for l in eachrow(pm_ntv(v).MAXIMAL_CONES)]
    n_ray = size(pm_ntv(v).RAYS, 1)
    maximal_cones = []
    for c in max_cones
        buffer = Vector{Int}(undef, n_ray)
        buffer = fill(0, n_ray)
        for k in c
            buffer[k] = 1
        end
        push!(maximal_cones, buffer)
    end
    # compute generators
    indeterminates = list_of_variables_of_cox_ring(v)
    gens = [];
    for i in 1 : length(maximal_cones)
        monom = 1
        for j in 1 : length(maximal_cones[i])
            monom = monom * indeterminates[j]^(1 - maximal_cones[i][j])
        end
        push!(gens, monom)
    end
    # return the ideal
    return ideal(gens)
end
export irrelevant_ideal


######################
# 2: Methods of ToricVarieties
######################


@doc Markdown.doc"""
    ith_betti_number(v::AbstractNormalToricVariety, i::Int)

Compute the i-th Betti number of the normal toric variety `v`.
"""
function ith_betti_number(v::AbstractNormalToricVariety, i::Int)
    if isodd(i)
        return 0
    end
    k = div(i, 2)
    f_vector = Vector{Int}(pm_ntv(v).F_VECTOR)
    pushfirst!(f_vector, 1)
    betti_number = sum((-1)^(i-k) * binomial(i,k) * f_vector[dim(v) - i + 1] for i=k:dim(v))
    return betti_number
    
end
export ith_betti_number


# @doc Markdown.doc"""
#     toric_ideal_binomial_generators(antv::AffineNormalToricVariety)
# 
# Get the exponent vectors corresponding to the generators of the toric ideal
# associated to the affine normal toric variety `antv`.
# 
# # Examples
# Take the cyclic quotient singularity corresponding to the pair of integers
# `(2,5)`.
# ```jldoctest
# julia> C = Oscar.positive_hull([-2 5; 1 0])
# A polyhedral cone in ambient dimension 2
# 
# julia> antv = AffineNormalToricVariety(C)
# A normal toric variety corresponding to a polyhedral fan in ambient dimension 2
# 
# julia> toric_ideal_binomial_generators(antv)
# pm::Matrix<long>
# -1 -1 2 1
# -1 0 3 -1
# 0 -1 -1 2
# ```
# """
# function toric_ideal_binomial_generators(antv::AffineNormalToricVariety)
#     pmntv = pm_ntv(antv)
#     result = pmntv.TORIC_IDEAL.BINOMIAL_GENERATORS
#     return result
# end
# export toric_ideal_binomial_generators
