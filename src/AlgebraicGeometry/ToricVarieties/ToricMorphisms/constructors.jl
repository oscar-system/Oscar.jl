####################################################
# 1: The Julia type for ToricMorphisms
####################################################

@attributes mutable struct ToricMorphism
    domain::AbstractNormalToricVariety
    grid_morphism::GrpAbFinGenMap
    image::AbstractNormalToricVariety
    codomain::AbstractNormalToricVariety
    ToricMorphism(domain, grid_morphism, image, codomain) = new(domain, grid_morphism, image, codomain)
end


####################################################
# 2: Generic constructors
####################################################



@doc raw"""
    toric_morphism(domain::AbstractNormalToricVariety, mapping_matrix::ZZMatrix, codomain::AbstractNormalToricVariety)

Construct the toric morphism with given domain and associated to the lattice
morphism given by the `mapping_matrix`.

If the codomain is left out, it will be determined whether the image of the
domain fan is itself a polyhedral fan. In that case the codomain is assumed to
be the associated toric variety.

All checks can be disabled with `check=false`.

# Examples
```jldoctest
julia> domain = projective_space(NormalToricVariety, 1)
Normal, non-affine, smooth, projective, gorenstein, fano, 1-dimensional toric variety without torusfactor

julia> codomain = hirzebruch_surface(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> mapping_matrix = matrix(ZZ, [0 1])
[0   1]

julia> toric_morphism(domain, mapping_matrix, codomain)
A toric morphism
```
"""
function toric_morphism(domain::AbstractNormalToricVariety, mapping_matrix::ZZMatrix, codomain::AbstractNormalToricVariety; check=true) 
  @req (nrows(mapping_matrix) > 0 && ncols(mapping_matrix) > 0) "The mapping matrix must not be empty"
  return toric_morphism(domain, hom(character_lattice(domain), character_lattice(codomain), mapping_matrix), codomain; check=check)
end
toric_morphism(domain::AbstractNormalToricVariety, mapping_matrix::Matrix{T}, codomain::AbstractNormalToricVariety; check=true) where {T <: IntegerUnion} =
  toric_morphism(domain, matrix(ZZ, mapping_matrix), codomain; check=check)
toric_morphism(domain::AbstractNormalToricVariety, mapping_matrix::Vector{Vector{T}}, codomain::AbstractNormalToricVariety; check=true) where {T <: IntegerUnion} =
  toric_morphism(domain, matrix(ZZ, mapping_matrix), codomain; check=check)
function toric_morphism(domain::AbstractNormalToricVariety, mapping_matrix::ZZMatrix) 
  @req (nrows(mapping_matrix) > 0 && ncols(mapping_matrix) > 0) "The mapping matrix must not be empty"
  return toric_morphism(domain, hom(character_lattice(domain), free_abelian_group(ncols(mapping_matrix)), mapping_matrix))
end
toric_morphism(domain::AbstractNormalToricVariety, mapping_matrix::Matrix{T}) where {T <: IntegerUnion} =
  toric_morphism(domain, matrix(ZZ, mapping_matrix))
toric_morphism(domain::AbstractNormalToricVariety, mapping_matrix::Vector{Vector{T}}) where {T <: IntegerUnion} =
  toric_morphism(domain, matrix(ZZ, mapping_matrix))


@doc raw"""
    function toric_morphism(domain::AbstractNormalToricVariety, grid_morphism::GrpAbFinGenMap, codomain::T=nothing) where {T <: Union{AbstractNormalToricVariety, Nothing}}

Construct the toric morphism from the `domain` to the `codomain` with map given
by the `grid_morphism`.

If the codomain is left out, it will be determined whether the image of the
domain fan is itself a polyhedral fan. In that case the codomain is assumed to
be the associated toric variety.

All checks can be disabled with `check=false`.

# Examples
```jldoctest
julia> domain = projective_space(NormalToricVariety, 1)
Normal, non-affine, smooth, projective, gorenstein, fano, 1-dimensional toric variety without torusfactor

julia> codomain = hirzebruch_surface(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> mapping_matrix = matrix(ZZ, [[0, 1]])
[0   1]

julia> grid_morphism = hom(character_lattice(domain), character_lattice(codomain), mapping_matrix)
Map with following data
Domain:
=======
Abelian group with structure: Z
Codomain:
=========
Abelian group with structure: Z^2

julia> toric_morphism(domain, grid_morphism, codomain)
A toric morphism
```
"""
function toric_morphism(domain::AbstractNormalToricVariety, grid_morphism::GrpAbFinGenMap, codomain::AbstractNormalToricVariety; check=true)
    # avoid empty mapping
    @req (nrows(matrix(grid_morphism)) > 0 && ncols(matrix(grid_morphism)) > 0) "The mapping matrix must not be empty"

    # check for a well-defined map
    @req nrows(matrix(grid_morphism)) == rank(character_lattice(domain)) "The number of rows of the mapping matrix must match the rank of the character lattice of the domain toric variety"

    # compute the image
    image_rays = matrix(ZZ, rays(domain)) * matrix(grid_morphism)
    image_rays = hcat([[Int(image_rays[i, j]) for i in 1:nrows(image_rays)] for j in 1:ncols(image_rays)]...)
    image = normal_toric_variety(polyhedral_fan(image_rays, ray_indices(maximal_cones(domain))))

    # compute the morphism
    @req ncols(matrix(grid_morphism)) == rank(character_lattice(codomain)) "The number of columns of the mapping matrix must match the rank of the character lattice of the codomain toric variety"
    if check
      codomain_cones = maximal_cones(codomain)
      image_cones = [positive_hull(matrix(ZZ, rays(c)) * matrix(grid_morphism)) for c in maximal_cones(domain)]
      for c in image_cones
        @req any(cc -> issubset(c, cc), codomain_cones) "Toric morphism not well-defined"
      end
    end
    return ToricMorphism(domain, grid_morphism, image, codomain)
end
function toric_morphism(domain::AbstractNormalToricVariety, grid_morphism::GrpAbFinGenMap)
    image_rays = matrix(ZZ, rays(domain)) * matrix(grid_morphism)
    image_rays = hcat([[Int(image_rays[i, j]) for i in 1:nrows(image_rays)] for j in 1:ncols(image_rays)]...)
    image = normal_toric_variety(polyhedral_fan(image_rays, ray_indices(maximal_cones(domain))))
    return ToricMorphism(domain, grid_morphism, image, image)
end


####################################################
# 3: Special constructors
####################################################

@doc raw"""
    toric_identity_morphism(variety::AbstractNormalToricVariety)

Construct the toric identity morphism from `variety` to `variety`.

# Examples
```jldoctest
julia> toric_identity_morphism(hirzebruch_surface(NormalToricVariety, 2))
A toric morphism
```
"""
function toric_identity_morphism(variety::AbstractNormalToricVariety)
    r = rank(character_lattice(variety))
    identity_matrix = matrix(ZZ, [[if i==j 1 else 0 end for j in 1:r] for i in 1:r])
    grid_morphism = hom(character_lattice(variety), character_lattice(variety), identity_matrix)
    return ToricMorphism(variety, grid_morphism, variety, variety)
end


####################################################
# 4: Addition and scalar multiplication of morphisms
####################################################

function Base.:+(tm1::ToricMorphism, tm2::ToricMorphism)
    @req domain(tm1) === domain(tm2) "The toric morphisms must have identical domains"
    @req codomain(tm1) === codomain(tm2) "The toric morphisms must have identical codomains"
    return toric_morphism(domain(tm1), grid_morphism(tm1) + grid_morphism(tm2), codomain(tm1))
end


function Base.:-(tm1::ToricMorphism, tm2::ToricMorphism)
    @req domain(tm1) === domain(tm2) "The toric morphisms must have identical domains"
    @req codomain(tm1) === codomain(tm2) "The toric morphisms must have identical codomains"
    return toric_morphism(domain(tm1), grid_morphism(tm1) - grid_morphism(tm2), codomain(tm1))
end


function Base.:*(c::T, tm::ToricMorphism) where T <: IntegerUnion
  new_grid_morphism = hom(domain(grid_morphism(tm)), codomain(grid_morphism(tm)), c * matrix(grid_morphism(tm)))
  return toric_morphism(domain(tm), new_grid_morphism, codomain(tm))
end


####################################################
# 5: Composition of toric morphisms
####################################################

function Base.:*(tm1::ToricMorphism, tm2::ToricMorphism)
    @req codomain(tm1) === domain(tm2) "The codomain of the first toric morphism must be identically the same as the domain of the second morphism"
    return toric_morphism(domain(tm1), grid_morphism(tm1) * grid_morphism(tm2), codomain(tm2))
end


####################################################
# 6: Equality and hash of toric morphisms
####################################################

function Base.:(==)(tm1::ToricMorphism, tm2::ToricMorphism)
    return domain(tm1) == domain(tm2) && codomain(tm1) == codomain(tm2) && grid_morphism(tm1) == grid_morphism(tm2)
end

function Base.hash(tm::ToricMorphism, h::UInt)
    b = 0x1a66f927cae2d409 % UInt
    h = hash(domain(tm), h)
    h = hash(codomain(tm), h)
    h = hash(grid_morphism(tm), h)
    return xor(h, b)
end


######################
# 7: Display
######################

function Base.show(io::IO, tm::ToricMorphism)
    join(io, "A toric morphism")
end
