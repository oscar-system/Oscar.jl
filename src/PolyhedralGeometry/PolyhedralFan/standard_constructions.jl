#TODO: inward/outward options? via polymake changes?

"""
    normal_fan(P::Polyhedron)

Return the normal fan of `P`. The maximal cones of the normal fan of `P` are
dual to the edge cones at the vertices of `P`.

# Examples
The rays of a normal fan of a cube point in every positive and negative unit
direction.
```jldoctest
julia> C = cube(3);

julia> NF = normal_fan(C)
Polyhedral fan in ambient dimension 3

julia> rays(NF)
6-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
 [-1, 0, 0]
 [0, 1, 0]
 [0, -1, 0]
 [0, 0, 1]
 [0, 0, -1]
```
"""
function normal_fan(P::Polyhedron{T}) where {T<:scalar_types}
  pmp = pm_object(P)
  pmnf = Polymake.fan.normal_fan(pmp)
  return PolyhedralFan{T}(pmnf, coefficient_field(P))
end

"""
    face_fan(P::Polyhedron)

Return the face fan of `P`. The polytope `P` has to contain the origin, then
the maximal cones of the face fan of `P` are the cones over the facets of `P`.

# Examples
By definition, this bounded polyhedron's number of facets equals the amount of
maximal cones of its face fan.
```jldoctest
julia> C = cross_polytope(3);

julia> FF = face_fan(C)
Polyhedral fan in ambient dimension 3

julia> n_maximal_cones(FF) == n_facets(C)
true
```
"""
function face_fan(P::Polyhedron{T}) where {T<:scalar_types}
  pmp = pm_object(P)
  pmff = Polymake.fan.face_fan(pmp)
  return PolyhedralFan{T}(pmff, coefficient_field(P))
end

###############################################################################
## Star subdivision
###############################################################################

@doc raw"""
    star_subdivision(PF::PolyhedralFan, exceptional_ray::AbstractVector{<:IntegerUnion})

Return the star subdivision of a polyhedral fan by a primitive element of
the underlying lattice. We follow the definition at the top of page 515 in
[CLS11](@cite).

# Examples
```jldoctest
julia> fan = normal_fan(simplex(3))
Polyhedral fan in ambient dimension 3

julia> exceptional_ray = [1, 1, 1];

julia> star = star_subdivision(fan, exceptional_ray)
Polyhedral fan in ambient dimension 3

julia> rays(star)
5-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
 [0, 1, 0]
 [0, 0, 1]
 [-1, -1, -1]
 [1, 1, 1]

julia> ray_indices(maximal_cones(star))
6×5 IncidenceMatrix
[2, 3, 5]
[1, 3, 5]
[1, 2, 5]
[2, 3, 4]
[1, 3, 4]
[1, 2, 4]
```
"""
function star_subdivision(
  Sigma::_FanLikeType, exceptional_ray::AbstractVector{<:IntegerUnion}
)

  # Check if exceptional_ray is primitive in ZZ^d, i.e. gcd(exceptional_ray)==1.
  @req ambient_dim(Sigma) == length(exceptional_ray) "New ray is not $(ambient_dim(Sigma))-dimensional"
  @req gcd(exceptional_ray) == 1 "The exceptional ray r is not a primitive element of the lattice Z^d with d = length(r)"
  @req lineality_dim(Sigma) == 0 "star_subdivision does not work for polyhedral fans with lineality."

  old_rays = matrix(ZZ, rays(Sigma))
  # In case the exceptional ray is an old ray.
  exceptional_ray_index = findfirst(
    i -> vec(old_rays[i, :]) == exceptional_ray, 1:nrows(old_rays)
  )
  new_rays = old_rays
  if isnothing(exceptional_ray_index)
    new_rays = vcat(old_rays, matrix(ZZ, [exceptional_ray]))
    exceptional_ray_index = n_rays(Sigma) + 1
  end
  mc_old = maximal_cones(IncidenceMatrix, Sigma)

  facet_normals = matrix(QQ, pm_object(Sigma).FACET_NORMALS)
  refinable_cones = _get_maximal_cones_containing_vector(Sigma, exceptional_ray)
  @req length(refinable_cones) > 0 "$exceptional_ray not contained in support of fan."
  new_cones = _get_refinable_facets(
    Sigma, exceptional_ray, refinable_cones, facet_normals, mc_old
  )
  for nc in new_cones
    push!(nc, exceptional_ray_index)
  end
  for i in 1:n_maximal_cones(Sigma)
    if !(i in refinable_cones)
      push!(new_cones, Vector{Int}(Polymake.row(mc_old, i)))
    end
  end

  return polyhedral_fan(
    coefficient_field(Sigma),
    IncidenceMatrix([nc for nc in new_cones]),
    new_rays;
    non_redundant=true,
  )
end

function _get_refinable_facets(
  Sigma::_FanLikeType,
  exceptional_ray::AbstractVector{<:IntegerUnion},
  refinable_cones::Vector{Int},
  facet_normals::MatElem,
  mc_old::IncidenceMatrix,
)
  new_cones = Vector{Int}[]
  v_facet_signs = _facet_signs(facet_normals, exceptional_ray)
  R = rays(Sigma)
  hd = pm_object(Sigma).HASSE_DIAGRAM
  hd_graph = Graph{Directed}(hd.ADJACENCY)
  hd_maximal_cones = inneighbors(hd_graph, hd.TOP_NODE + 1)
  Sfacets = pm_object(Sigma).MAXIMAL_CONES_FACETS
  for mc_index in refinable_cones
    mc_indices = Polymake.row(mc_old, mc_index)
    mc_facet_indices = Polymake.findnz(Sfacets[mc_index, :])[1]
    mc_hd_index = hd_maximal_cones[findfirst(
      i ->
        Polymake.to_one_based_indexing(Polymake._get_entry(hd.FACES, i - 1)) == mc_indices,
      hd_maximal_cones,
    )]
    refinable_facets = _get_refinable_facets_of_cone(
      mc_hd_index, facet_normals, hd, hd_graph, v_facet_signs, R, mc_facet_indices
    )
    append!(new_cones, refinable_facets)
    # If all facets contain exceptional_ray, then the current maximal cone is just
    # exceptional_ray.
    length(refinable_facets) > 0 || push!(new_cones, Vector{Int}(mc_indices))
  end
  return unique(new_cones)
end

function _get_refinable_facets_of_cone(
  mc_hd_index::Int,
  facet_normals::MatElem,
  hd,
  hd_graph::Graph{Directed},
  v_facet_signs::AbstractVector,
  R::AbstractVector{<:RayVector},
  mcfi::Vector{Int},
)
  refinable_facets = Vector{Int}[]
  mc_indices = Vector{Int}(
    Polymake.to_one_based_indexing(Polymake._get_entry(hd.FACES, mc_hd_index - 1))
  )
  for fc_index in inneighbors(hd_graph, mc_hd_index)
    fc_indices = Polymake.to_one_based_indexing(Polymake._get_entry(hd.FACES, fc_index - 1))
    length(fc_indices) > 0 || return refinable_facets # The only facet was 0
    inner_ray = sum([R[i] for i in fc_indices])
    fc_facet_signs = _facet_signs(facet_normals, inner_ray)
    if (!_check_containment_via_facet_signs(v_facet_signs[mcfi], fc_facet_signs[mcfi]))
      push!(refinable_facets, Vector{Int}(fc_indices))
    end
  end
  return refinable_facets
end

_facet_signs(F::MatElem, v::AbstractVector) = sign.(Int, F * v)[:, 1]

function _check_containment_via_facet_signs(smaller::Vector{Int}, bigger::Vector{Int})
  for a in zip(smaller, bigger)
    p = prod(a)
    if p == 0
      a[1] == 0 || return false
    end
    p >= 0 || return false # Both facet vectors must point in the same direction.
  end
  return true
end

function _get_maximal_cones_containing_vector(
  Sigma::_FanLikeType, v::AbstractVector{<:RationalUnion}
)
  # Make sure these are computed as otherwise performance degrades.
  pm_object(Sigma).FACET_NORMALS
  pm_object(Sigma).MAXIMAL_CONES_FACETS
  return findall(mc -> v in mc, maximal_cones(Sigma))
end

@doc raw"""
    star_subdivision(PF::PolyhedralFan, n::Int)

Return the star subdivision of a polyhedral fan at its n-th torus orbit.
Note that this torus orbit need not be maximal. We follow definition 3.3.17
of [CLS11](@cite).

# Examples
```jldoctest
julia> star = star_subdivision(normal_fan(simplex(3)), 1)
Polyhedral fan in ambient dimension 3

julia> rays(star)
5-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
 [0, 1, 0]
 [0, 0, 1]
 [-1, -1, -1]
 [1, 1, 1]

julia> ray_indices(maximal_cones(star))
6×5 IncidenceMatrix
[2, 3, 5]
[1, 3, 5]
[1, 2, 5]
[2, 3, 4]
[1, 3, 4]
[1, 2, 4]
```
"""
function star_subdivision(Sigma::_FanLikeType, n::Int)
  cones_Sigma = cones(Sigma)
  tau = Polymake.row(cones_Sigma, n)
  R = matrix(ZZ, rays(Sigma))
  exceptional_ray = vec(sum([R[i, :] for i in tau]))
  exceptional_ray = exceptional_ray ./ gcd(exceptional_ray)
  return star_subdivision(Sigma, exceptional_ray)
end

@doc raw"""
    transform(F::_FanLikeType, A::AbstractMatrix; check=true) where

Get the image of a fan `F` under the matrix `A`. The default is to check
whether the images of the maximal cones really form a fan.
"""
function transform(
  F::_FanLikeType,
  A::Union{AbstractMatrix{<:Union{Number,FieldElem}},MatElem{<:FieldElem}};
  check::Bool=true,
)
  @req ncols(A) == ambient_dim(F) "Incompatible dimension of fan and transformation matrix"
  OT = _scalar_type_to_polymake(_get_scalar_type(F))
  return _transform(F, Polymake.Matrix{OT}(A); check)
end

function _transform(F::_FanLikeType, A::Polymake.Matrix; check::Bool)
  OT = _scalar_type_to_polymake(_get_scalar_type(F))
  FT = typeof(F)
  R = pm_object(F).RAYS * transpose(A)
  L = pm_object(F).LINEALITY_SPACE * transpose(A)
  MC = pm_object(F).MAXIMAL_CONES
  opt = Polymake.OptionSet(Dict(["lineality_space" => L]))
  if check
    result = Polymake.fan.check_fan(R, MC, opt)
    return FT(result, coefficient_field(F))
  else
    result = Polymake.fan.PolyhedralFan{OT}(; RAYS=R, LINEALITY_SPACE=L, MAXIMAL_CONES=MC)
    return FT(result, coefficient_field(F))
  end
end

###############################################################################
## Cartesian/Direct product
###############################################################################

@doc raw"""
    *(PF1::PolyhedralFan{QQFieldElem}, PF2::PolyhedralFan{QQFieldElem})

Return the Cartesian/direct product of two polyhedral fans.

# Examples
```jldoctest
julia> normal_fan(simplex(2))*normal_fan(simplex(3))
Polyhedral fan in ambient dimension 5
```
"""
function Base.:*(PF1::PolyhedralFan{QQFieldElem}, PF2::PolyhedralFan{QQFieldElem})
  prod = Polymake.fan.product(pm_object(PF1), pm_object(PF2))
  return PolyhedralFan{QQFieldElem}(prod, QQ)
end

#################################################################################
## Hyperplane arrangements
#################################################################################

@doc raw"""
    arrangement_polynomial([ring::MPolyRing{<: FieldElem},] A::MatElem{<: FieldElem})

Given some $A\in\mathbb{F}^{n\times d},$ return the product of the linear forms
corresponding to the rows.

Let $A$ be a $n\times d$ matrix with entries from a field $\mathbb{F}$.
The rows of $A$  are the normal vectors for a hyperplane arrangement
$$\mathcal{A} = \{H_{1},\dots,H_{n}:H_{i}\subset \mathbb{F}^{d}\}.$$
We have $H_{i} = V(\alpha_{i})$, where $\alpha_{i}\in\mathbb{F}[x_{1},\dots,x_{d}]$
is a linear form whose coefficients are the entries of the $i$th row.

Then we have $$\cup_{H_{i}\in\mathcal{A}}H_{i} = V(\Pi^{n}_{i=1}\alpha_{i}).$$

Optionally one can select to use columns instead of rows in the following way:
```
arrangement_polynomial(...  ; hyperplanes=:in_cols)
```

# Example using standard ring and then custom ring.
```jldoctest
julia> A = matrix(QQ,[1 2 5//2; 0 0 1; 2 3 2; 1//2 3 5; 3 1 2; 7 8 1])
[   1   2   5//2]
[   0   0      1]
[   2   3      2]
[1//2   3      5]
[   3   1      2]
[   7   8      1]

julia> factor(arrangement_polynomial(A))
(1//4) * (2*x1 + 3*x2 + 2*x3) * (7*x1 + 8*x2 + x3) * (x1 + 6*x2 + 10*x3) * (2*x1 + 4*x2 + 5*x3) * x3 * (3*x1 + x2 + 2*x3)

julia> R,_ = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> factor(arrangement_polynomial(R, A))
(1//4) * (2*x + 3*y + 2*z) * (7*x + 8*y + z) * (x + 6*y + 10*z) * (2*x + 4*y + 5*z) * z * (3*x + y + 2*z)
```

To use the columns instead, proceed in the following way:
```jldoctest
julia> A = matrix(QQ,[1 0 2 1//2 3 7;2 0 3 3 1 8;5//2 1 2 5 2 1]);

julia> factor(arrangement_polynomial(A; hyperplanes=:in_cols))
(1//4) * (2*x1 + 3*x2 + 2*x3) * (7*x1 + 8*x2 + x3) * (x1 + 6*x2 + 10*x3) * (2*x1 + 4*x2 + 5*x3) * x3 * (3*x1 + x2 + 2*x3)
```
"""
function arrangement_polynomial(A::MatElem{<:FieldElem}; hyperplanes=:in_rows)
  F = base_ring(A)
  nvars = hyperplanes == :in_rows ? ncols(A) : nrows(A)
  P, _ = polynomial_ring(F, nvars; cached=false)
  return arrangement_polynomial(P, A; hyperplanes)
end
function arrangement_polynomial(
  ring::MPolyRing{<:FieldElem}, A::MatElem{<:FieldElem}; hyperplanes=:in_rows
)
  if hyperplanes == :in_cols
    return arrangement_polynomial(ring, transpose(A))
  else
    @req dim(ring) == ncols(A) "dimension of ring must be number of rows of input matrix"
    @req base_ring(A) == coefficient_ring(ring) "entries of input matrix must be coefficients from ring"
    x = gens(ring)
    return prod(A * x)
  end
end
function arrangement_polynomial(A::AbstractVector{<:AbstractVector{<:FieldElem}})
  @req length(A) > 0 "At least one hyperplane needs to be provided"
  nvars = length(A[1])
  @req all(x -> length(x) == nvars, A) "All hyperplanes need to have the same dimension"
  F = parent(first(first(A)))
  return arrangement_polynomial(matrix(F, A))
end
function arrangement_polynomial(
  ring::MPolyRing{<:FieldElem}, A::AbstractVector{<:AbstractVector{<:FieldElem}}
)
  @req length(A) > 0 "At least one hyperplane needs to be provided"
  nvars = length(A[1])
  @req all(x -> length(x) == nvars, A) "All hyperplanes need to have the same dimension"
  F = parent(first(first(A)))
  return arrangement_polynomial(ring, matrix(F, A))
end
function arrangement_polynomial(A::AbstractMatrix{<:FieldElem}; hyperplanes=:in_rows)
  F = parent(first(A))
  return arrangement_polynomial(matrix(F, A); hyperplanes)
end
function arrangement_polynomial(
  ring::MPolyRing{<:FieldElem}, A::AbstractMatrix{<:FieldElem}; hyperplanes=:in_rows
)
  F = parent(first(A))
  return arrangement_polynomial(ring, matrix(F, A); hyperplanes)
end
