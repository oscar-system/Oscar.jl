@doc raw"""
    domain(tm::ToricMorphism)

Return the domain of the toric morphism `tm`.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> domain(toric_identity_morphism(F4))
Normal toric variety
```
"""
domain(tm::ToricMorphism) = tm.domain


@doc raw"""
    codomain(tm::ToricMorphism)

Return the codomain of the toric morphism `tm`.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> codomain(toric_identity_morphism(F4))
Normal toric variety
```
"""
codomain(tm::ToricMorphism) = tm.codomain


@doc raw"""
    lattice_homomorphism(tm::ToricMorphism)

Return the underlying homomorphism of lattices of one parameter subgroups of the toric morphism `tm` as in Definition 3.3.1 of [CLS11](@cite).

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> lattice_homomorphism(toric_identity_morphism(F4))
Map
  from Z^2
  to Z^2
```
"""
lattice_homomorphism(tm::ToricMorphism) = tm.lattice_homomorphism


@doc raw"""
    morphism_on_torusinvariant_weil_divisor_group(tm::ToricMorphism)

For a given toric morphism `tm`, this method computes the corresponding
map of the torusinvariant Weil divisors.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> morphism_on_torusinvariant_weil_divisor_group(toric_identity_morphism(F4))
Map
  from Z^4
  to Z^4
```
"""
@attr FinGenAbGroupHom function morphism_on_torusinvariant_weil_divisor_group(tm::ToricMorphism)
    d = domain(tm)
    cod = codomain(tm)
    cod_rays = matrix(ZZ, rays(cod))
    images = matrix(ZZ, rays(d)) * matrix(lattice_homomorphism(tm))
    mapping_matrix = matrix(ZZ, zeros(ZZ, torsion_free_rank(torusinvariant_weil_divisor_group(cod)), 0))
    for i in 1:nrows(images)
      v = [images[i,k] for k in 1:ncols(images)]
      j = findfirst(j -> v in maximal_cones(cod)[j], 1:n_maximal_cones(cod))
      m = reduce(vcat, [Int(ray_indices(maximal_cones(cod))[j, k]) * cod_rays[k:k, :] for k in 1:n_rays(cod)])
      mapping_matrix = hcat(mapping_matrix, solve(transpose(m), transpose(images[i:i, :]); side = :right))
    end
    return hom(torusinvariant_weil_divisor_group(d), torusinvariant_weil_divisor_group(cod), transpose(mapping_matrix))
end


@doc raw"""
    morphism_on_torusinvariant_cartier_divisor_group(tm::ToricMorphism)

For a given toric morphism `tm`, this method computes the corresponding
map of the Cartier divisors.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> morphism_on_torusinvariant_cartier_divisor_group(toric_identity_morphism(F4))
Map
  from Z^4
  to Z^4
```
"""
@attr FinGenAbGroupHom function morphism_on_torusinvariant_cartier_divisor_group(tm::ToricMorphism)
    domain_variety = domain(tm)
    codomain_variety = codomain(tm)
    domain_embedding = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(domain_variety)
    morphism_on_weil_divisors = morphism_on_torusinvariant_weil_divisor_group(tm)
    codomain_post_inverse = postinverse(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(codomain_variety))
    return domain_embedding * morphism_on_weil_divisors * codomain_post_inverse
end


@doc raw"""
    morphism_on_class_group(tm::ToricMorphism)

For a given toric morphism `tm`, this method computes the corresponding
map of the Class groups.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> morphism_on_class_group(toric_identity_morphism(F4))
Map
  from Z^2
  to Z^2
```
"""
@attr FinGenAbGroupHom function morphism_on_class_group(tm::ToricMorphism)
    domain_variety = domain(tm)
    codomain_variety = codomain(tm)
    domain_preinverse = preinverse(map_from_torusinvariant_weil_divisor_group_to_class_group(domain_variety))
    morphism_on_weil_divisors = morphism_on_torusinvariant_weil_divisor_group(tm)
    codomain_projection = map_from_torusinvariant_weil_divisor_group_to_class_group(codomain_variety)
    return domain_preinverse * morphism_on_weil_divisors * codomain_projection
end


@doc raw"""
    morphism_on_picard_group(tm::ToricMorphism)

For a given toric morphism `tm`, this method computes the corresponding
map of the Picard groups.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> morphism_on_picard_group(toric_identity_morphism(F4))
Map
  from Z^2
  to Z^2
```
"""
@attr FinGenAbGroupHom function morphism_on_picard_group(tm::ToricMorphism)
    domain_variety = domain(tm)
    codomain_variety = codomain(tm)
    domain_preinverse = preinverse(map_from_torusinvariant_cartier_divisor_group_to_picard_group(domain_variety))
    morphism_on_cartier_divisors = morphism_on_torusinvariant_cartier_divisor_group(tm)
    codomain_projection = map_from_torusinvariant_cartier_divisor_group_to_picard_group(codomain_variety)
    return domain_preinverse * morphism_on_cartier_divisors * codomain_projection
end


########################################################################
# Functionality for the associated CoveredSchemeMorphism               #
########################################################################

@doc raw"""
    covering_morphism(f::ToricMorphism)

For a given toric morphism `tm`, we can compute the corresponding
morphism of covered schemes. The following demonstrates this for the
blow-up morphism of a blow-up of the projective space.

# Examples
```jldoctest
julia> IP2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> bl = blow_up(IP2, [1, 1]);

julia> cov_bl = covering_morphism(bl);

julia> domain(cov_bl)
Covering
  described by patches
    1: normal toric variety
    2: normal toric variety
    3: normal toric variety
    4: normal toric variety
  in the coordinate(s)
    1: [x_1_1, x_2_1]
    2: [x_1_2, x_2_2]
    3: [x_1_3, x_2_3]
    4: [x_1_4, x_2_4]

julia> codomain(cov_bl)
Covering
  described by patches
    1: normal toric variety
    2: normal toric variety
    3: normal toric variety
  in the coordinate(s)
    1: [x_1_1, x_2_1]
    2: [x_1_2, x_2_2]
    3: [x_1_3, x_2_3]
```
"""
@attr CoveringMorphism function covering_morphism(f::ToricMorphism)
  # TODO: If f is a blowup morphism, we can simplify 
  # the matchings of cones below.
  X = domain(f)
  Y = codomain(f)

  # Find the image cones
  codomain_cones = maximal_cones(Y)
  domain_cones = maximal_cones(X)
  A = matrix(lattice_homomorphism(f))
  image_cones = [positive_hull(matrix(ZZ, rays(c)) * A) for c in domain_cones]

  # construct the corresponding morphism of rings
  morphism_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  domain_cov = default_covering(X) # ordering of the patches must be the same as `maximal_cones`
  codomain_cov = default_covering(Y)
  for i in 1:n_maximal_cones(X)
    U = domain_cov[i] # The corresponding chart in the domain
    k = findfirst(x-> is_subset(image_cones[i], x), codomain_cones)
    V = codomain_cov[k] # The chart in the codomain whose cone contains the image of the cone of U

    wc1 = weight_cone(U)
    hb_U = hilbert_basis(wc1) # corresponds to the variables of OO(U)
    wc2 = weight_cone(V)
    hb_V = hilbert_basis(wc2)

    At = transpose(A) # the matrix for the dual of the lattice_map
    hb_V_mat = matrix(ZZ, hb_V)
    hb_V_img = hb_V_mat * At
    hb_U_mat = matrix(ZZ, hb_U)
    #sol = solve(hb_U_mat, hb_V_img; side = :left) # deprecated! It was giving vectors with negative entries.
    id_mat = identity_matrix(ZZ, nrows(hb_U_mat))
    offset_mat = zero_matrix(ZZ, nrows(hb_U_mat), 1)
    sol_list = Vector{ZZMatrix}()
    for k in 1:nrows(hb_V_img)
      b = transpose(hb_V_img[k:k, :])
      push!(sol_list, solve_mixed(ZZMatrix, transpose(hb_U_mat), b, id_mat, offset_mat))
    end
    # For some weird reason `solve_mixed` solves A*x = b, but returns x as a 1xn-matrix (not nx1).
    sol = vcat(sol_list...)

    @assert sol*hb_U_mat == hb_V_img
    @assert all(>=(0), sol)

    # assemble the monomials where the variables of OO(V) are mapped
    imgs = [prod(gens(OO(U))[k]^sol[i, k] for k in 1:ngens(OO(U)); init=one(OO(U))) for i in 1:nrows(sol)]
    morphism_dict[U] = morphism(U, V, imgs)
  end
  return CoveringMorphism(domain_cov, codomain_cov, morphism_dict, check=false)
end

# Some helper functions for conversion. 
# These are here, because we didn't know better and documentation on how to convert polymake 
# objects is too poor to be used by us. Please improve if you know how to!
function _my_mult(u::PointVector{ZZRingElem}, A::ZZMatrix)
  m = length(u)
  m == nrows(A) || error("sizes incompatible")
  n = ncols(A)
  result = zero_matrix(ZZ, 1, n)
  for k in 1:n
    result[1, k] = sum(u[i]*A[i, k] for i in 1:m; init=zero(ZZ))
  end
  return result
end

function _to_ZZ_matrix(u::PointVector{ZZRingElem})
  n = length(u)
  result = zero_matrix(ZZ, 1, n)
  for i in 1:n
    result[1, i] = u[i]
  end
  return result
end

@attr CoveredSchemeMorphism function underlying_morphism(f::ToricMorphism)
  return CoveredSchemeMorphism(domain(f), codomain(f), covering_morphism(f))
end
