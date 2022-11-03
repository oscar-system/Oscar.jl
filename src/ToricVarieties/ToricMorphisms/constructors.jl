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
export ToricMorphism


####################################################
# 2: Generic constructors
####################################################


@doc Markdown.doc"""
    ToricMorphism(domain::AbstractNormalToricVariety, mapping_matrix::Vector{Vector{T}}, codomain::T2=nothing) where {T <: IntegerUnion, T2 <: Union{AbstractNormalToricVariety, Nothing}}

Construct the toric morphism with given domain and associated to the lattice morphism given by the `mapping_matrix`.
As optional argument, the codomain of the morphism can be specified.

# Examples
```jldoctest
julia> domain = projective_space(NormalToricVariety, 1)
A normal, non-affine, smooth, projective, gorenstein, fano, 1-dimensional toric variety without torusfactor

julia> mapping_matrix = [[0, 1]]
1-element Vector{Vector{Int64}}:
 [0, 1]

julia> ToricMorphism(domain, mapping_matrix)
A toric morphism
```
"""
function ToricMorphism(domain::AbstractNormalToricVariety, mapping_matrix::Vector{Vector{T}}, codomain::T2=nothing) where {T <: IntegerUnion, T2 <: Union{AbstractNormalToricVariety, Nothing}}
    (length(mapping_matrix) > 0 && length(mapping_matrix[1]) > 0) || throw(ArgumentError("The mapping matrix must not be empty"))
    if codomain == nothing
      return ToricMorphism(domain, hom(character_lattice(domain), free_abelian_group(length(mapping_matrix[1])), matrix(ZZ, mapping_matrix)), codomain)
    else
      return ToricMorphism(domain, hom(character_lattice(domain), character_lattice(codomain), matrix(ZZ, mapping_matrix)), codomain)
    end
end
export ToricMorphism


@doc Markdown.doc"""
    ToricMorphism(domain::AbstractNormalToricVariety, mapping_matrix::Matrix{T}, codomain::T2=nothing) where {T <: IntegerUnion, T2 <: Union{AbstractNormalToricVariety, Nothing}}

Construct the toric morphism with given domain and associated to the lattice morphism given by the `mapping_matrix`.

# Examples
```jldoctest
julia> domain = projective_space(NormalToricVariety, 1)
A normal, non-affine, smooth, projective, gorenstein, fano, 1-dimensional toric variety without torusfactor

julia> mapping_matrix = [0 1]
1×2 Matrix{Int64}:
 0  1

julia> ToricMorphism(domain, mapping_matrix)
A toric morphism
```
"""
function ToricMorphism(domain::AbstractNormalToricVariety, mapping_matrix::Matrix{T}, codomain::T2=nothing) where {T <: IntegerUnion, T2 <: Union{AbstractNormalToricVariety, Nothing}}
    (nrows(mapping_matrix) > 0 && ncols(mapping_matrix) > 0) || throw(ArgumentError("The mapping matrix must not be empty"))
    if codomain == nothing
      return ToricMorphism(domain, hom(character_lattice(domain), free_abelian_group(ncols(mapping_matrix)), matrix(ZZ, mapping_matrix)), codomain)
    else
      return ToricMorphism(domain, hom(character_lattice(domain), character_lattice(codomain), matrix(ZZ, mapping_matrix)), codomain)
    end
end
export ToricMorphism


@doc Markdown.doc"""
    ToricMorphism(domain::AbstractNormalToricVariety, mapping_matrix::fmpz_mat, codomain::T=nothing) where {T <: Union{AbstractNormalToricVariety, Nothing}}

Construct the toric morphism with given domain and associated to the lattice morphism given by the `mapping_matrix`.

# Examples
```jldoctest
julia> domain = projective_space(NormalToricVariety, 1)
A normal, non-affine, smooth, projective, gorenstein, fano, 1-dimensional toric variety without torusfactor

julia> codomain = hirzebruch_surface(2)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> mapping_matrix = matrix(ZZ, [0 1])
[0   1]

julia> ToricMorphism(domain, mapping_matrix, codomain)
A toric morphism
```
"""
function ToricMorphism(domain::AbstractNormalToricVariety, mapping_matrix::fmpz_mat, codomain::T=nothing) where {T <: Union{AbstractNormalToricVariety, Nothing}}
    (nrows(mapping_matrix) > 0 && ncols(mapping_matrix) > 0) || throw(ArgumentError("The mapping matrix must not be empty"))
    if codomain == nothing
      return ToricMorphism(domain, hom(character_lattice(domain), free_abelian_group(ncols(mapping_matrix)), mapping_matrix), codomain)
    else
      return ToricMorphism(domain, hom(character_lattice(domain), character_lattice(codomain), mapping_matrix), codomain)
    end
end
export ToricMorphism


@doc Markdown.doc"""
    function ToricMorphism(domain::AbstractNormalToricVariety, grid_morphism::GrpAbFinGenMap, codomain::T=nothing) where {T <: Union{AbstractNormalToricVariety, Nothing}}

Construct the toric morphism from the `domain` to the `codomain` with map given by the `grid_morphism`.

# Examples
```jldoctest
julia> domain = projective_space(NormalToricVariety, 1)
A normal, non-affine, smooth, projective, gorenstein, fano, 1-dimensional toric variety without torusfactor

julia> codomain = hirzebruch_surface(2)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

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

julia> ToricMorphism(domain, grid_morphism, codomain)
A toric morphism
```
"""
function ToricMorphism(domain::AbstractNormalToricVariety, grid_morphism::GrpAbFinGenMap, codomain::T=nothing) where {T <: Union{AbstractNormalToricVariety, Nothing}}
    # avoid empty mapping
    (nrows(matrix(grid_morphism)) > 0 && ncols(matrix(grid_morphism)) > 0) || throw(ArgumentError("The mapping matrix must not be empty"))

    # check for a well-defined map
    if nrows(matrix(grid_morphism)) !== rank(character_lattice(domain))
      throw(ArgumentError("The number of rows of the mapping matrix must match the rank of the character lattice of the domain toric variety"))
    end

    # compute the image
    image_rays = matrix(ZZ, rays(domain)) * matrix(grid_morphism)
    image_rays = hcat([[Int(image_rays[i, j]) for i in 1:nrows(image_rays)] for j in 1:ncols(image_rays)]...)
    image = NormalToricVariety(PolyhedralFan(image_rays, ray_indices(maximal_cones(domain))))

    # compute the morphism
    if codomain == nothing
      return ToricMorphism(domain, grid_morphism, image, image)
    else
      if ncols(matrix(grid_morphism)) !== rank(character_lattice(codomain))
        throw(ArgumentError("The number of columns of the mapping matrix must match the rank of the character lattice of the codomain toric variety"))
      end
      codomain_cones = maximal_cones(codomain)
      image_cones = [positive_hull(matrix(ZZ, rays(c)) * matrix(grid_morphism)) for c in maximal_cones(domain)]
      for c in image_cones
        any([intersect(c, ic) == c for ic in image_cones]) || throw(ArgumentError("Toric morphism not well-defined"))
      end
      return ToricMorphism(domain, grid_morphism, image, codomain)
    end
end
export ToricMorphism


####################################################
# 3: Special constructors
####################################################

@doc Markdown.doc"""
    ToricIdentityMorphism(variety::AbstractNormalToricVariety)

Construct the toric identity morphism from `variety` to `variety`.

# Examples
```jldoctest
julia> ToricIdentityMorphism(hirzebruch_surface(2))
A toric morphism
```
"""
function ToricIdentityMorphism(variety::AbstractNormalToricVariety)
    r = rank(character_lattice(variety))
    identity_matrix = matrix(ZZ, [[if i==j 1 else 0 end for j in 1:r] for i in 1:r])
    grid_morphism = hom(character_lattice(variety), character_lattice(variety), identity_matrix)
    return ToricMorphism(variety, grid_morphism, variety, variety)
end
export ToricIdentityMorphism


####################################################
# 4: Addition and scalar multiplication of morphisms
####################################################

function Base.:+(tm1::ToricMorphism, tm2::ToricMorphism)
    if domain(tm1) !== domain(tm2)
        throw(ArgumentError("The toric morphism must be defined with identically the same domain"))
    end
    if codomain(tm1) !== codomain(tm2)
        throw(ArgumentError("The toric morphism must be defined with identically the same codomain"))
    end
    return ToricMorphism(domain(tm1), grid_morphism(tm1) + grid_morphism(tm2), codomain(tm1))
end


function Base.:-(tm1::ToricMorphism, tm2::ToricMorphism)
    if domain(tm1) !== domain(tm2)
        throw(ArgumentError("The toric morphism must be defined with identically the same domain"))
    end
    if codomain(tm1) !== codomain(tm2)
        throw(ArgumentError("The toric morphism must be defined with identically the same codomain"))
    end
    return ToricMorphism(domain(tm1), grid_morphism(tm1) - grid_morphism(tm2), codomain(tm1))
end


function Base.:*(c::T, tm::ToricMorphism) where {T <: IntegerUnion, T2 <: Union{AbstractNormalToricVariety, Nothing}}
  new_grid_morphism = hom(domain(grid_morphism(tm)), codomain(grid_morphism(tm)), c * matrix(grid_morphism(tm)))
  return ToricMorphism(domain(tm), new_grid_morphism, codomain(tm))
end


####################################################
# 5: Composition of toric morphisms
####################################################

function Base.:*(tm1::ToricMorphism, tm2::ToricMorphism)
    if codomain(tm1) !== domain(tm2)
        throw(ArgumentError("The codomain of the first toric morphism must be identically the same as the domain of the second morphism"))
    end
    return ToricMorphism(domain(tm1), grid_morphism(tm1) * grid_morphism(tm2), codomain(tm2))
end


####################################################
# 6: Equality of toric morphisms
####################################################

function Base.:(==)(tm1::ToricMorphism, tm2::ToricMorphism)
    return domain(tm1) == domain(tm2) && codomain(tm1) == codomain(tm2) && grid_morphism(tm1) == grid_morphism(tm2)
end


######################
# 6: Display
######################

function Base.show(io::IO, tm::ToricMorphism)
    join(io, "A toric morphism")
end
