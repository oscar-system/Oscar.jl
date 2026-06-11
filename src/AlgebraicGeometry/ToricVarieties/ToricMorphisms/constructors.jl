####################################################
# 1: The Julia type for ToricMorphisms
####################################################

@attributes mutable struct ToricMorphism{
  DomainType<:AbsCoveredScheme,
  CodomainType<:AbsCoveredScheme,
  BaseMorphismType,
} <: AbsCoveredSchemeMorphism{
  DomainType,
  CodomainType,
  BaseMorphismType,
  ToricMorphism,
}
  domain::NormalToricVarietyType
  lattice_homomorphism::FinGenAbGroupHom
  codomain::NormalToricVarietyType
  function ToricMorphism(domain, lattice_homomorphism, codomain)
    result = new{typeof(domain),typeof(codomain),Nothing}(
      domain, lattice_homomorphism, codomain
    )
  end
end

####################################################
# 2: Generic constructors
####################################################

@doc raw"""
    toric_morphism(domain::NormalToricVarietyType, mapping_matrix::ZZMatrix, codomain::NormalToricVarietyType; check=true)

Construct the toric morphism with given domain and associated to the lattice
morphism given by the `mapping_matrix`.

If the codomain is left out, it will be determined whether the image of the
domain fan is itself a polyhedral fan. In that case the codomain is assumed to
be the associated toric variety.

All checks can be disabled with `check=false`.

# Examples
```jldoctest
julia> domain = projective_space(NormalToricVariety, 1)
Normal toric variety

julia> codomain = hirzebruch_surface(NormalToricVariety, 2)
Normal toric variety

julia> mapping_matrix = matrix(ZZ, [0 1])
[0   1]

julia> toric_morphism(domain, mapping_matrix, codomain)
Toric morphism
```
"""
function toric_morphism(
  domain::NormalToricVarietyType,
  mapping_matrix::ZZMatrix,
  codomain::NormalToricVarietyType;
  check=true,
)
  @req (nrows(mapping_matrix) > 0 && ncols(mapping_matrix) > 0) "The mapping matrix must not be empty"
  return toric_morphism(
    domain,
    hom(
      lattice_of_one_parameter_subgroups(domain),
      lattice_of_one_parameter_subgroups(codomain),
      mapping_matrix,
    ),
    codomain;
    check=check,
  )
end
toric_morphism(
  domain::NormalToricVarietyType,
  mapping_matrix::Matrix{T},
  codomain::NormalToricVarietyType;
  check=true,
) where {T<:IntegerUnion} = toric_morphism(
  domain, matrix(ZZ, mapping_matrix), codomain; check=check
)
toric_morphism(
  domain::NormalToricVarietyType,
  mapping_matrix::Vector{Vector{T}},
  codomain::NormalToricVarietyType;
  check=true,
) where {T<:IntegerUnion} = toric_morphism(
  domain, matrix(ZZ, mapping_matrix), codomain; check=check
)
function toric_morphism(domain::NormalToricVarietyType, mapping_matrix::ZZMatrix)
  @req (nrows(mapping_matrix) > 0 && ncols(mapping_matrix) > 0) "The mapping matrix must not be empty"
  return toric_morphism(
    domain,
    hom(
      lattice_of_one_parameter_subgroups(domain),
      free_abelian_group(ncols(mapping_matrix)),
      mapping_matrix,
    ),
  )
end
toric_morphism(
  domain::NormalToricVarietyType, mapping_matrix::Matrix{T}
) where {T<:IntegerUnion} = toric_morphism(
  domain, matrix(ZZ, mapping_matrix)
)
toric_morphism(
  domain::NormalToricVarietyType, mapping_matrix::Vector{Vector{T}}
) where {T<:IntegerUnion} = toric_morphism(
  domain, matrix(ZZ, mapping_matrix)
)

@doc raw"""
    toric_morphism(domain::NormalToricVarietyType, lattice_homomorphism::FinGenAbGroupHom, codomain::NormalToricVarietyType; check=true)

Construct the toric morphism from `domain` to `codomain` induced by the
homomomorphism of lattices of one-parameter subgroups
`lattice_homomorphism` as in Definition 3.3.1 of [CLS11](@cite).

If the codomain is left out, it will be determined whether the image of the
domain fan is itself a polyhedral fan. In that case the codomain is assumed to
be the associated toric variety.

All checks can be disabled with `check=false`.

# Examples
```jldoctest
julia> domain = projective_space(NormalToricVariety, 1)
Normal toric variety

julia> codomain = hirzebruch_surface(NormalToricVariety, 2)
Normal toric variety

julia> mapping_matrix = matrix(ZZ, [[0, 1]])
[0   1]

julia> lattice_homomorphism = hom(lattice_of_one_parameter_subgroups(domain), lattice_of_one_parameter_subgroups(codomain), mapping_matrix)
Map
  from Z
  to Z^2

julia> toric_morphism(domain, lattice_homomorphism, codomain)
Toric morphism
```
"""
function toric_morphism(
  domain::NormalToricVarietyType,
  lattice_homomorphism::FinGenAbGroupHom,
  codomain::NormalToricVarietyType;
  check=true,
)
  # avoid empty mapping
  @req (nrows(matrix(lattice_homomorphism)) > 0 && ncols(matrix(lattice_homomorphism)) > 0) "The mapping matrix must not be empty"

  # check for a well-defined map
  @req nrows(matrix(lattice_homomorphism)) ==
    torsion_free_rank(lattice_of_one_parameter_subgroups(domain)) "The number of rows of the mapping matrix must match the rank of the lattice of one-parameter subgroups of the domain toric variety"

  # compute the morphism
  @req ncols(matrix(lattice_homomorphism)) ==
    torsion_free_rank(lattice_of_one_parameter_subgroups(codomain)) "The number of columns of the mapping matrix must match the rank of the lattice of one-parameter subgroups of the codomain toric variety"
  if check
    codomain_cones = maximal_cones(codomain)
    image_cones = [
      positive_hull(matrix(ZZ, rays(c)) * matrix(lattice_homomorphism)) for
      c in maximal_cones(domain)
    ]
    for c in image_cones
      @req any(cc -> issubset(c, cc), codomain_cones) "Toric morphism not well-defined"
    end
  end
  return ToricMorphism(domain, lattice_homomorphism, codomain)
end
function toric_morphism(
  domain::NormalToricVarietyType, lattice_homomorphism::FinGenAbGroupHom; check=true
)
  image = transform(domain, matrix(lattice_homomorphism); check=check)
  return ToricMorphism(domain, lattice_homomorphism, image)
end

####################################################
# 3: Special constructors
####################################################

@doc raw"""
    toric_identity_morphism(variety::NormalToricVarietyType)

Construct the toric identity morphism from `variety` to `variety`.

# Examples
```jldoctest
julia> toric_identity_morphism(hirzebruch_surface(NormalToricVariety, 2))
Toric morphism
```
"""
function toric_identity_morphism(variety::NormalToricVarietyType)
  r = torsion_free_rank(lattice_of_one_parameter_subgroups(variety))
  identity_matrix = matrix(ZZ, [[
    if i == j
      1
    else
      0
    end for j in 1:r
  ] for i in 1:r])

  # Despite the notation, the homomorphism is between lattices of
  # one-parameter subgroups
  lattice_homomorphism = hom(
    lattice_of_one_parameter_subgroups(variety),
    lattice_of_one_parameter_subgroups(variety),
    identity_matrix,
  )
  return ToricMorphism(variety, lattice_homomorphism, variety)
end

####################################################
# 4: Addition and scalar multiplication of morphisms
####################################################

function Base.:+(tm1::ToricMorphism, tm2::ToricMorphism)
  @req domain(tm1) === domain(tm2) "The toric morphisms must have identical domains"
  @req codomain(tm1) === codomain(tm2) "The toric morphisms must have identical codomains"
  return toric_morphism(
    domain(tm1), lattice_homomorphism(tm1) + lattice_homomorphism(tm2), codomain(tm1)
  )
end

function Base.:-(tm1::ToricMorphism, tm2::ToricMorphism)
  @req domain(tm1) === domain(tm2) "The toric morphisms must have identical domains"
  @req codomain(tm1) === codomain(tm2) "The toric morphisms must have identical codomains"
  return toric_morphism(
    domain(tm1), lattice_homomorphism(tm1) - lattice_homomorphism(tm2), codomain(tm1)
  )
end

function Base.:*(c::T, tm::ToricMorphism) where {T<:IntegerUnion}
  new_lattice_homomorphism = hom(
    domain(lattice_homomorphism(tm)),
    codomain(lattice_homomorphism(tm)),
    c * matrix(lattice_homomorphism(tm)),
  )
  return toric_morphism(domain(tm), new_lattice_homomorphism, codomain(tm))
end

####################################################
# 5: Composition of toric morphisms
####################################################

function Base.:*(tm1::ToricMorphism, tm2::ToricMorphism)
  @req codomain(tm1) === domain(tm2) "The codomain of the first toric morphism must be identically the same as the domain of the second morphism"
  return toric_morphism(
    domain(tm1), lattice_homomorphism(tm1) * lattice_homomorphism(tm2), codomain(tm2)
  )
end

####################################################
# 6: Equality and hash of toric morphisms
####################################################

function Base.:(==)(tm1::ToricMorphism, tm2::ToricMorphism)
  return domain(tm1) === domain(tm2) && codomain(tm1) === codomain(tm2) &&
         lattice_homomorphism(tm1) == lattice_homomorphism(tm2)
end

function Base.hash(tm::ToricMorphism, h::UInt)
  b = 0x1a66f927cae2d409 % UInt
  h = hash(domain(tm), h)
  h = hash(codomain(tm), h)
  h = hash(lattice_homomorphism(tm), h)
  return xor(h, b)
end

######################
# 7: Display
######################

function Base.show(io::IO, tm::ToricMorphism)
  join(io, "Toric morphism")
end

Base.show(io::IO, ::MIME"text/plain", tm::ToricMorphism) = Base.show(pretty(io), tm)
