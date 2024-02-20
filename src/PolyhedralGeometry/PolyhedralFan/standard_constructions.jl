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
function normal_fan(P::Polyhedron{T}) where T<:scalar_types
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
function face_fan(P::Polyhedron{T}) where T<:scalar_types
   pmp = pm_object(P)
   pmff = Polymake.fan.face_fan(pmp)
   return PolyhedralFan{T}(pmff, coefficient_field(P))
end


###############################################################################
## Star subdivision
###############################################################################


@doc raw"""
    star_subdivision(PF::PolyhedralFan, new_ray::AbstractVector{<:IntegerUnion})

Return the star subdivision of a polyhedral fan by a primitive element of
the underlying lattice. We follow the definition at the top of page 515 in
[CLS11](@cite).

# Examples
```jldoctest
julia> fan = normal_fan(simplex(3))
Polyhedral fan in ambient dimension 3

julia> new_ray = [1, 1, 1];

julia> star = star_subdivision(fan, new_ray)
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
function star_subdivision(Sigma::_FanLikeType, new_ray::AbstractVector{<:IntegerUnion})
  
  # Check if new_ray is primitive in ZZ^d, i.e. gcd(new_ray)==1.
  @req ambient_dim(Sigma) == length(new_ray) "New ray is not $(ambient_dim(Sigma))-dimensional"
  @req gcd(new_ray) == 1 "The new ray r is not a primitive element of the lattice Z^d with d = length(r)"
  @req lineality_dim(Sigma) == 0 "star_subdivision does not work for polyhedral fans with lineality."
  
  old_rays = matrix(ZZ, rays(Sigma))
  # In case the new ray is an old ray.
  new_ray_index = findfirst(i->vec(old_rays[i,:])==new_ray, 1:nrows(old_rays))
  new_rays = old_rays
  if isnothing(new_ray_index)
    new_rays = vcat(old_rays, matrix(ZZ, [new_ray]))
    new_ray_index = n_rays(Sigma)+1
  end
  mc_old = maximal_cones(IncidenceMatrix, Sigma)
  
  facet_normals = pm_object(Sigma).FACET_NORMALS
  refinable_cones = _get_maximal_cones_containing_vector(Sigma, new_ray)
  @req length(refinable_cones)>0 "$new_ray not contained in support of fan."
  new_cones = _get_refinable_facets(Sigma, new_ray, refinable_cones, facet_normals, mc_old)
  for nc in new_cones
    push!(nc, new_ray_index)
  end
  for i in 1:n_maximal_cones(Sigma)
    if !(i in refinable_cones)
      push!(new_cones, Vector{Int}(Polymake.row(mc_old, i)))
    end
  end
  
  return polyhedral_fan(coefficient_field(Sigma), IncidenceMatrix([nc for nc in new_cones]), new_rays; non_redundant=true)
end

function _get_refinable_facets(Sigma::_FanLikeType, new_ray::AbstractVector{<:IntegerUnion}, refinable_cones::Vector{Int}, facet_normals::AbstractMatrix, mc_old::IncidenceMatrix)
  new_cones = Vector{Int}[]
  v_facet_signs = _facet_signs(facet_normals, new_ray)
  R = pm_object(Sigma).RAYS
  hd = pm_object(Sigma).HASSE_DIAGRAM
  hd_graph = Graph{Directed}(hd.ADJACENCY)
  hd_maximal_cones = inneighbors(hd_graph, hd.TOP_NODE+1)
  Sfacets = pm_object(Sigma).MAXIMAL_CONES_FACETS
  for mc_index in refinable_cones
    mc_indices = Polymake.row(mc_old, mc_index)
    mc_facet_indices = Polymake.findnz(Sfacets[mc_index,:])[1]
    mc_hd_index = hd_maximal_cones[findfirst(i -> Polymake.to_one_based_indexing(Polymake._get_entry(hd.FACES, i-1)) == mc_indices, hd_maximal_cones)]
    refinable_facets = _get_refinable_facets_of_cone(mc_hd_index, facet_normals, hd, hd_graph, v_facet_signs, R, mc_facet_indices)
    append!(new_cones, refinable_facets)
    # If all facets contain new_ray, then the current maximal cone is just
    # new_ray.
    length(refinable_facets) > 0 || push!(new_cones, Vector{Int}(mc_indices))
  end
  return unique(new_cones)
end

function _get_refinable_facets_of_cone(mc_hd_index::Int, facet_normals::AbstractMatrix, hd, hd_graph::Graph{Directed}, v_facet_signs::AbstractVector, R::AbstractMatrix, mcfi::Vector{Int})
  refinable_facets = Vector{Int}[]
  mc_indices = Vector{Int}(Polymake.to_one_based_indexing(Polymake._get_entry(hd.FACES, mc_hd_index-1)))
  for fc_index in inneighbors(hd_graph, mc_hd_index)
    fc_indices = Polymake.to_one_based_indexing(Polymake._get_entry(hd.FACES, fc_index-1))
    length(fc_indices) > 0 || return refinable_facets # The only facet was 0
    inner_ray = sum([R[i,:] for i in fc_indices])
    fc_facet_signs = _facet_signs(facet_normals, inner_ray)
    if(!_check_containment_via_facet_signs(v_facet_signs[mcfi], fc_facet_signs[mcfi]))
      push!(refinable_facets, Vector{Int}(fc_indices))
    end
  end
  return refinable_facets
end

# FIXME: Small workaround, since sign does not work for polymake types.
_int_sign(e) = e>0 ? 1 : (e<0 ? -1 : 0)
_facet_signs(F::AbstractMatrix, v::AbstractVector{<:IntegerUnion}) = [_int_sign(e) for e in (F*Polymake.Vector{Polymake.Integer}(v))]
_facet_signs(F::AbstractMatrix, v::AbstractVector) = [_int_sign(e) for e in (F*v)]
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
function _get_maximal_cones_containing_vector(Sigma::_FanLikeType, v::AbstractVector{<:IntegerUnion})
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
  @req length(tau) > 1 "Cannot subdivide cone $n as it is generated by a single ray"
  R = matrix(ZZ, rays(Sigma))
  newray = vec(sum([R[i,:] for i in tau]))
  newray = newray ./ gcd(newray)
  return star_subdivision(Sigma, newray)
end



@doc raw"""
    transform(F::_FanLikeType, A::AbstractMatrix; check=true) where

Get the image of a fan `F` under the matrix `A`. The default is to check
whether the images of the maximal cones really form a fan.
"""
function transform(F::_FanLikeType, A::Union{AbstractMatrix{<:Union{Number, FieldElem}}, MatElem{<:FieldElem}}; check::Bool=true)
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
    result = Polymake.fan.PolyhedralFan{OT}(RAYS=R, LINEALITY_SPACE=L, MAXIMAL_CONES=MC)
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
