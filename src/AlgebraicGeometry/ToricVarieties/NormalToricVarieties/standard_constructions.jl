@doc raw"""
    affine_space(::Type{NormalToricVariety}, d::Int)

Construct the (toric) affine space of dimension `d`.

# Examples
```jldoctest
julia> affine_space(NormalToricVariety, 2)
Normal toric variety
```
"""
function affine_space(::Type{NormalToricVariety}, d::Int)
  @req d >= 0 "Dimension must be non-negative"
  C = positive_hull(identity_matrix(ZZ, d))
  variety = normal_toric_variety(C)
  return variety
end


@doc raw"""
    affine_space(::Type{NormalToricVariety}, var_symbols::AbstractVector{<:VarName})

Construct a (toric) affine space, specifying the variable names.

# Examples
```jldoctest
julia> affine_space(NormalToricVariety, [:x, :y, :z])
Normal toric variety

julia> affine_space(NormalToricVariety, ["x", "y", "z"])
Normal toric variety

julia> affine_space(NormalToricVariety, ['x', 'y', 'z'])
Normal toric variety
```
"""
function affine_space(::Type{NormalToricVariety}, var_names::AbstractVector{<:VarName})
  variety = affine_space(NormalToricVariety, length(var_names))
  set_coordinate_names(variety, var_names)
  return variety
end


@doc raw"""
    projective_space(::Type{NormalToricVariety}, d::Int)

Construct the projective space of dimension `d`.

# Examples
```jldoctest
julia> projective_space(NormalToricVariety, 2)
Normal toric variety
```
"""
function projective_space(::Type{NormalToricVariety}, d::Int)
  @req d >= 0 "Dimension must be non-negative"
  f = normal_fan(Oscar.simplex(d))
  pm_object = Polymake.fulton.NormalToricVariety(Oscar.pm_object(f))
  variety = NormalToricVariety(pm_object)
  weights = matrix(ZZ,ones(Int,d+1,1))
  tiwd_group = free_abelian_group(d+1)
  cl_group = free_abelian_group(1)
  char_lattice = free_abelian_group(d)
  mat = transpose(matrix(ZZ, rays(variety)))
  set_attribute!(variety, :map_from_character_lattice_to_torusinvariant_weil_divisor_group, hom(char_lattice, tiwd_group,  mat))
  set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(tiwd_group, cl_group, weights))
  return variety
end


@doc raw"""
    weighted_projective_space(::Type{NormalToricVariety}, w::Vector{T}) where {T <: IntegerUnion}

Construct the weighted projective space corresponding to the weights `w`.

# Examples
```jldoctest
julia> weighted_projective_space(NormalToricVariety, [2, 3, 1])
Normal toric variety
```
"""
function weighted_projective_space(::Type{NormalToricVariety}, w::Vector{T}) where {T <: IntegerUnion}
  @req length(w) > 0 "At least one weight must be specified"
  if all(is_one, w)
    return projective_space(NormalToricVariety, length(w)-1)
  end
  pmntv = Polymake.fulton.weighted_projective_space(w)
  variety = NormalToricVariety(pmntv)
  weights = matrix(ZZ,hcat(w))
  tiwd_group = free_abelian_group(length(w))
  cl_group = free_abelian_group(1)
  char_lattice = free_abelian_group(length(w)-1)
  mat = transpose(matrix(ZZ, rays(variety)))
  set_attribute!(variety, :map_from_character_lattice_to_torusinvariant_weil_divisor_group, hom(char_lattice, tiwd_group,  mat))
  set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(tiwd_group, cl_group, weights))
  return variety
end


@doc raw"""
    hirzebruch_surface(::Type{NormalToricVariety}, r::Int)

Construct the `r`-th Hirzebruch surface.

# Examples
```jldoctest
julia> hirzebruch_surface(NormalToricVariety, 5)
Normal toric variety
```
"""
function hirzebruch_surface(::Type{NormalToricVariety}, r::Int)
  fan_rays = [1 0; 0 1; -1 r; 0 -1]
  cones = IncidenceMatrix([[1, 2], [2, 3], [3, 4], [4, 1]])
  variety = normal_toric_variety(cones, fan_rays; non_redundant = true)
  weights = matrix(ZZ, [1 0; 0 1; 1 0; r 1])
  vars = ["t1", "x1", "t2", "x2"]
  set_coordinate_names(variety, vars)
  tiwd_group = free_abelian_group(4)
  cl_group = free_abelian_group(2)
  char_lattice = free_abelian_group(2)
  mat = transpose(matrix(ZZ, rays(variety)))
  set_attribute!(variety, :map_from_character_lattice_to_torusinvariant_weil_divisor_group, hom(char_lattice, tiwd_group,  mat))
  set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(tiwd_group, cl_group, weights))
  return variety
end


@doc raw"""
    del_pezzo_surface(::Type{NormalToricVariety}, b::Int)

Construct the del Pezzo surface with `b` blowups for `b` at most 3.

# Examples
```jldoctest
julia> del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety
```
"""
function del_pezzo_surface(::Type{NormalToricVariety}, b::Int)
  # check for valid input
  @req b >= 0 "Number of blowups for construction of del Pezzo surfaces must be non-negative"
  @req b <= 3 "Del Pezzo surfaces with more than three blowups are realized as subvarieties of toric ambient spaces. This is currently not supported"

  # special case of projective space
  if b == 0
    return projective_space(NormalToricVariety, 2)
  end

  # construct the "true" toric del Pezzo surfaces
  if b == 1
    fan_rays = [1 0; 0 1; -1 -1; 1 1]
    cones = IncidenceMatrix([[1, 4], [2, 4], [1, 3], [2, 3]])
  end
  if b == 2
    fan_rays = [1 0; 0 1; -1 -1; 1 1; 0 -1]
    cones = IncidenceMatrix([[1, 4], [2, 4], [1, 5], [5, 3], [2, 3]])
  end
  if b == 3
    fan_rays = [1 0; 0 1; -1 -1; 1 1; 0 -1; -1 0]
    cones = IncidenceMatrix([[1, 4], [2, 4], [1, 5], [5, 3], [2, 6], [6, 3]])
  end
  variety = normal_toric_variety(cones, fan_rays; non_redundant = true)
  tiwd_group = free_abelian_group(b+3)
  cl_group = free_abelian_group(b+1)
  char_lattice = free_abelian_group(2)
  mat = transpose(matrix(ZZ, rays(variety)))
  set_attribute!(variety, :map_from_character_lattice_to_torusinvariant_weil_divisor_group, hom(char_lattice, tiwd_group, mat))
  if b == 1
    weights = matrix(ZZ, [1 1; 1 1; 1 0; 0 -1])
    set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(tiwd_group, cl_group, weights))
  end
  if b == 2
    weights = matrix(ZZ, [1 1 1; 1 1 0; 1 0 1; 0 -1 0; 0 0 -1])
    set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(tiwd_group, cl_group, weights))
  end
  if b == 3
    weights = matrix(ZZ, [1 1 1 0; 1 1 0 1; 1 0 1 1; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1])
    set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, hom(tiwd_group, cl_group, weights))
  end
  vars = ["x1", "x2", "x3", "e1", "e2", "e3"]
  set_coordinate_names(variety, vars[1:(3 + b)])
  return variety
end

@doc raw"""
    Base.:*(v::NormalToricVarietyType, w::NormalToricVarietyType)

Return the Cartesian/direct product of two normal toric varieties `v` and `w`.

By default, we prepend an "x" to all homogeneous coordinate names of the first factor
`v` and a "y" to all homogeneous coordinate names of the second factor `w`. This default
can be overwritten by invoking `set_coordinate_names` after creating the variety
(cf. [`set_coordinate_names(v::NormalToricVarietyType, coordinate_names::Vector{String})`](@ref)).

*Important*: Recall that the coordinate names can only be changed as long as the toric
variety in question is not finalized (cf. [`is_finalized(v::NormalToricVarietyType)`](@ref)).

Crucially, the order of the homogeneous coordinates is not shuffled. To be more
specific, assume that `v` has ``n_1`` and `w` has ``n_2`` homogeneous coordinates. Then
`v * w` has ``n_1 + n_2`` homogeneous coordinates. The first ``n_1`` of these coordinates
are those of `v` and appear in the very same order as they do for `v`. The remaining
``n_2`` homogeneous coordinates are those of `w` and appear in the very same order as they
do for `w`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> v1 = P2 * P2
Normal toric variety

julia> cox_ring(v1)
Multivariate polynomial ring in 6 variables over QQ graded by
  xx1 -> [1 0]
  xx2 -> [1 0]
  xx3 -> [1 0]
  yx1 -> [0 1]
  yx2 -> [0 1]
  yx3 -> [0 1]

julia> v2 = P2 * P2
Normal toric variety

julia> set_coordinate_names(v2, ["x1", "x2", "x3", "y1", "y2", "y3"])


julia> cox_ring(v2)
Multivariate polynomial ring in 6 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [1 0]
  x3 -> [1 0]
  y1 -> [0 1]
  y2 -> [0 1]
  y3 -> [0 1]
```
"""
function Base.:*(v::NormalToricVarietyType, w::NormalToricVarietyType)
    product = normal_toric_variety(polyhedral_fan(v)*polyhedral_fan(w))
    set_coordinate_names(product, vcat(["x$(i)" for i in coordinate_names(v)], ["y$(i)" for i in coordinate_names(w)]))
    return product
end



@doc raw"""
    normal_toric_variety_from_star_triangulation(P::Polyhedron)

Return a toric variety that was obtained from a fine regular
star triangulation of the lattice points of the polyhedron P.
This is particularly useful when the lattice points of the
polyhedron in question admit many triangulations.

# Examples
```jldoctest
julia> P = convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1])
Polyhedron in ambient dimension 3

julia> v = normal_toric_variety_from_star_triangulation(P)
Normal toric variety
```
"""
function normal_toric_variety_from_star_triangulation(P::Polyhedron)
  # Find position of origin in the lattices points of the polyhedron P
  pts = matrix(ZZ, lattice_points(P))
  indices = findall(k -> is_zero_row(pts, k), 1:nrows(pts))
  @req length(indices) == 1 "Polyhedron must contain origin (exactly once)"

  # Change order of lattice points s.t. zero is the first point
  tmp = pts[1:1,:]
  pts[1:1,:] = pts[indices[1]:indices[1],:]
  pts[indices[1]:indices[1],:] = tmp

  # Find one triangulation and turn it into the maximal cones of the toric variety in question. Note that:
  # (a) needs to be converted to incidence matrix
  # (b) one has to remove origin from list of indices (as removed above)
  max_cones = IncidenceMatrix([[c[i]-1 for i in 2:length(c)] for c in _find_full_star_triangulation(pts)])

  # Rays are all but the zero vector at the first position of pts
  integral_rays = reduce(vcat, [pts[k:k,:] for k in 2:nrows(pts)])

  # construct the variety
  return normal_toric_variety(max_cones, integral_rays; non_redundant = true)
end



@doc raw"""
    normal_toric_varieties_from_star_triangulations(P::Polyhedron)

Return the list of toric varieties obtained from fine regular
star triangulations of the polyhedron P. With this we can
compute the two phases of the famous conifold transition.

# Examples
```jldoctest
julia> P = convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1])
Polyhedron in ambient dimension 3

julia> (v1, v2) = normal_toric_varieties_from_star_triangulations(P)
2-element Vector{NormalToricVariety}:
 Normal toric variety
 Normal toric variety

julia> stanley_reisner_ideal(v1)
Ideal generated by
  x1*x4

julia> stanley_reisner_ideal(v2)
Ideal generated by
  x2*x3
```
"""
function normal_toric_varieties_from_star_triangulations(P::Polyhedron)
  # find position of origin in the lattices points of the polyhedron P
  pts = matrix(ZZ, lattice_points(P))
  indices = findall(k -> is_zero_row(pts, k), 1:nrows(pts))
  @req length(indices) == 1 "Polyhedron must contain origin (exactly once)"

  # change order of lattice points s.t. zero is the first point
  tmp = pts[1,:]
  pts[1,:] = pts[indices[1],:]
  pts[indices[1],:] = tmp

  # triangulate the points - note that we enforce full, i.e. want all points to be rays
  trias = star_triangulations(pts; full = true)

  # rays are all but the zero vector at the first position of pts
  integral_rays = reduce(vcat, [pts[k:k,:] for k in 2:nrows(pts)])

  # trias contains the max_cones as list of lists
  # (a) needs to be converted to incidence matrix
  # (b) one has to remove origin from list of indices (as removed above)
  max_cones = [IncidenceMatrix([[c[i]-1 for i in 2:length(c)] for c in t]) for t in trias]

  # construct the varieties
  return [normal_toric_variety(cones, integral_rays; non_redundant = true) for cones in max_cones]
end



@doc raw"""
    normal_toric_variety_from_glsm(charges::ZZMatrix)

Return one toric variety with the desired
GLSM charges. This can be particularly useful provided that
there are many such toric varieties.

# Examples
```jldoctest
julia> charges = [[1, 1, 1]]
1-element Vector{Vector{Int64}}:
 [1, 1, 1]

julia> normal_toric_variety_from_glsm(charges)
Normal toric variety
```

For convenience, we also support:
- `normal_toric_variety_from_glsm(charges::Vector{Vector{Int}})`
- `normal_toric_variety_from_glsm(charges::Vector{Vector{ZZRingElem}})`
"""
function normal_toric_variety_from_glsm(charges::ZZMatrix)

  # find the ray generators
  G1 = free_abelian_group(ncols(charges))
  G2 = free_abelian_group(nrows(charges))
  map = hom(G1, G2, transpose(charges))
  ker = kernel(map)
  embedding = snf(ker[1])[2] * ker[2]
  integral_rays = transpose(embedding.map)

  # identify the points to be triangulated
  pts = matrix(ZZ, zeros(nrows(integral_rays)+1, ncols(integral_rays)))
  pts[2:end, :] = integral_rays

  # construct varieties
  triang = _find_full_star_triangulation(pts)
  max_cones = IncidenceMatrix([[c[i]-1 for i in 2:length(c)] for c in triang])
  variety = normal_toric_variety(max_cones, integral_rays; non_redundant = true)

  # set the attributes and return the variety
  char_lattice = free_abelian_group(ambient_dim(variety))
  mat = transpose(matrix(ZZ, rays(variety)))
  set_attribute!(variety, :map_from_character_lattice_to_torusinvariant_weil_divisor_group, hom(char_lattice, G1,  mat))
  set_attribute!(variety, :map_from_torusinvariant_weil_divisor_group_to_class_group, map)
  return variety
end
normal_toric_variety_from_glsm(charges::Vector{Vector{T}}) where {T <: IntegerUnion} = normal_toric_variety_from_glsm(matrix(ZZ, charges))



@doc raw"""
    normal_toric_varieties_from_glsm(charges::ZZMatrix)

Return all toric variety with the desired
GLSM charges. This computation may take a long time if
there are many such toric varieties.

# Examples
```jldoctest
julia> charges = [[1, 1, 1]]
1-element Vector{Vector{Int64}}:
 [1, 1, 1]

julia> normal_toric_varieties_from_glsm(charges)
1-element Vector{NormalToricVariety}:
 Normal toric variety

julia> varieties = normal_toric_varieties_from_glsm(matrix(ZZ, [1 2 3 4 6 0; -1 -1 -2 -2 -3 1]))
1-element Vector{NormalToricVariety}:
 Normal toric variety

julia> cox_ring(varieties[1])
Multivariate polynomial ring in 6 variables over QQ graded by
  x1 -> [1 -1]
  x2 -> [2 -1]
  x3 -> [3 -2]
  x4 -> [4 -2]
  x5 -> [6 -3]
  x6 -> [0 1]
```

For convenience, we also support:
- `normal_toric_varieties_from_glsm(charges::Vector{Vector{Int}})`
- `normal_toric_varieties_from_glsm(charges::Vector{Vector{ZZRingElem}})`
"""
function normal_toric_varieties_from_glsm(charges::ZZMatrix)
  # find the ray generators
  G1 = free_abelian_group(ncols(charges))
  G2 = free_abelian_group(nrows(charges))
  map = hom(G1, G2, transpose(charges))
  ker = kernel(map)
  embedding = snf(ker[1])[2] * ker[2]
  initial_rays = transpose(embedding.map)

  # identify the points to be triangulated
  pts = zeros(QQ, nrows(initial_rays), ncols(charges)-nrows(charges))
  for i in 1:nrows(initial_rays)
    pts[i, :] = [ZZRingElem(c) for c in initial_rays[i, :]]
  end
  zero = [0 for i in 1:ncols(charges)-nrows(charges)]
  pts = vcat(matrix(QQ, transpose(zero)), matrix(QQ, pts))

  # construct varieties
  integral_rays = reduce(vcat, [pts[k:k,:] for k in 2:nrows(pts)])
  max_cones = [IncidenceMatrix([[c[i]-1 for i in 2:length(c)] for c in t]) for t in star_triangulations(pts; full = true)]
  varieties = [normal_toric_variety(cones, integral_rays; non_redundant = true) for cones in max_cones]

  # set the map from Div_T -> Cl to the desired matrix
  for v in varieties
    char_lattice = free_abelian_group(ambient_dim(v))
    mat = transpose(matrix(ZZ, rays(v)))
    set_attribute!(v, :map_from_character_lattice_to_torusinvariant_weil_divisor_group, hom(char_lattice, G1,  mat))
    set_attribute!(v, :map_from_torusinvariant_weil_divisor_group_to_class_group, map)
  end

  # return the varieties
  return varieties
end
normal_toric_varieties_from_glsm(charges::Vector{Vector{T}}) where {T <: IntegerUnion} = normal_toric_varieties_from_glsm(matrix(ZZ, charges))
