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
   return PolyhedralFan{T}(pmnf)
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

julia> n_maximal_cones(FF) == nfacets(C)
true
```
"""
function face_fan(P::Polyhedron{T}) where T<:scalar_types
   pmp = pm_object(P)
   pmff = Polymake.fan.face_fan(pmp)
   return PolyhedralFan{T}(pmff)
end


###############################################################################
## Star subdivision
###############################################################################


@doc raw"""
    star_subdivision(PF::PolyhedralFan, new_ray::Vector{Int64})

Return the star subdivision of a polyhedral fan by a primitive element of
the underlying lattice. We follow the definition at the top of page 515 in
[CLS11](@cite).

# Examples
```jldoctest
julia> fan = normal_fan(simplex(3))
Polyhedral fan in ambient dimension 3

julia> rays(fan)
4-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
 [0, 1, 0]
 [0, 0, 1]
 [-1, -1, -1]

julia> ray_indices(maximal_cones(fan))
4×4 IncidenceMatrix
[1, 2, 3]
[2, 3, 4]
[1, 3, 4]
[1, 2, 4]

julia> new_ray = [1, 1, 1]
3-element Vector{Int64}:
 1
 1
 1

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
[1, 3, 4]
[1, 3, 5]
[2, 3, 4]
[1, 2, 5]
[1, 2, 4]
```
"""
function star_subdivision(Sigma::_FanLikeType{T}, new_ray::Vector{Int64}) where T<:scalar_types
  
  # Check if new_ray is primitive.
  # Note that new_ray is primitive iff v/n is a lattice point in Z^d (d = length(new_ray)) iff n = 1 or n == -1.
  # For our application, the lattice in question is always Z^d for a suitable integer d.
  # -> We check if gcd(all entries of new_ray) is 1.
  @req ambient_dim(Sigma) == length(new_ray) "New ray cannot be a primitive element"
  @req gcd(new_ray) == 1 "The new ray r is not a primitive element of the lattice Z^d with d = length(r)"
  
  # Construct new rays
  old_rays = matrix(ZZ,rays(Sigma))
  new_rays = vcat(old_rays, matrix(ZZ, [new_ray]))
  
  mc_orig = maximal_cones(IncidenceMatrix, Sigma)

  refinable_cones = _get_maximal_cones_containing_vector(Sigma, new_ray)
  new_cones = Set{Vector{Int}}()
  for rc in refinable_cones
    new_cones = Base.union(new_cones, Set{Vector{Int}}(_refine_maximal_cone(Sigma, new_ray, rc)))
  end
  for nc in new_cones
    push!(nc, nrays(Sigma)+1)
  end
  for i in 1:n_maximal_cones(Sigma)
    if !(i in refinable_cones)
      push!(new_cones, Vector{Int}(Polymake.row(mc_orig)))
    end
  end
  
  return polyhedral_fan(T, new_rays, IncidenceMatrix([nc for nc in new_cones]); non_redundant=true)
  
end

# FIXME: Small workaround, since sign does not work for polymake types.
_int_sign(e) = e>0 ? 1 : (e<0 ? -1 : 0)
_facet_signs(F::AbstractMatrix, v::AbstractVector) = [_int_sign(e) for e in (F*v)]
function _check_containment_via_facet_signs(smaller::Vector{Int}, bigger::Vector{Int}; zero_unused=false)
  # zero_unused indicates whether zero entries in bigger indicate that the cone
  # does not use the facet. Otherwise this is treated like an equation.
  for a in zip(smaller, bigger)
    if !(zero_unused && a[2] == 0)
      p = prod(a)
      if p == 0
        a[1] == 0 || return false
      elseif p < 0
        return false
      end
    end
  end
  return true
end
function _get_maximal_cones_containing_vector(Sigma::_FanLikeType{T}, v::Vector{Int64}) where T<:scalar_types
  maximal_cones_signs = Matrix{Int}(pm_object(Sigma).MAXIMAL_CONES_FACETS)
  v_facet_signs = _facet_signs(pm_object(Sigma).FACET_NORMALS, v)
  return filter(i -> _check_containment_via_facet_signs(v_facet_signs, maximal_cones_signs[i,:]; zero_unused=true), 1:n_maximal_cones(Sigma))
end
function _refine_maximal_cone(Sigma::_FanLikeType{T}, new_ray::Vector, mc_index::Int) where T<:scalar_types
  v_facet_signs = _facet_signs(pm_object(Sigma).FACET_NORMALS, new_ray)
  R = pm_object(Sigma).RAYS
  hd = pm_object(Sigma).HASSE_DIAGRAM
  hdg = Graph{Directed}(hd.ADJACENCY)
  mc_indices = Polymake.row(pm_object(Sigma).MAXIMAL_CONES, mc_index)
  mcs = inneighbors(hdg, hd.TOP_NODE+1)
  mc_hd_index = mcs[findfirst(i -> Polymake.to_one_based_indexing(Polymake._get_entry(hd.FACES, i-1)) == mc_indices, mcs)]
  result = Vector{Int}[]
  for fc_index in inneighbors(hdg, mc_hd_index)
    fc_indices = Polymake.to_one_based_indexing(Polymake._get_entry(hd.FACES, fc_index-1))
    inner_ray = sum([R[i,:] for i in fc_indices])
    fc_facet_signs = _facet_signs(pm_object(Sigma).FACET_NORMALS, inner_ray)
    if(!_check_containment_via_facet_signs(v_facet_signs, fc_facet_signs))
      push!(result, Vector{Int}(fc_indices))
    end
  end
  return result
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
[1, 3, 5]
[2, 3, 5]
[1, 2, 5]
[2, 3, 4]
[1, 3, 4]
[1, 2, 4]
```
"""
function star_subdivision(Sigma::_FanLikeType{T}, n::Int) where T<:scalar_types
  
  # check if n-th cone exist
  @req n <= n_cones(Sigma) "Cannot subdivide cone $n as it does not exist"
  
  # check if n-th cone can be subdivided
  cones_Sigma = cones(Sigma)
  tau = Polymake.row(cones_Sigma, n)
  @req length(tau) > 1 "Cannot subdivide cone $n as it is generated by a single ray"
  
  # check if nth cone is smooth
  @req _cone_is_smooth(Sigma, tau) "Cannot subdivide maximal cone $n as it is not smooth"
  
  # compute new ray
  R = Polymake.common.primitive(Oscar.pm_object(Sigma).RAYS)
  newindex = size(R,1) + 1
  newray = sum([R[i,:] for i in tau])
  
  # compute new rays
  newrays = [R; transpose(newray)]
  
  # compute the new maximal cones
  maximal_cones_Sigma = ray_indices(maximal_cones(Sigma))
  tau_ray_indices = Polymake.row(cones_Sigma, n)
  newmaxcones = (Vector{Int})[]
  for i in 1:n_maximal_cones(Sigma)
    indices_ith_max_cone = Polymake.row(maximal_cones_Sigma, i)
    if issubset(tau_ray_indices, indices_ith_max_cone)
      @req _cone_is_smooth(Sigma, indices_ith_max_cone) "All cones containing sigma need to be smooth"
      for subset in subsets(Vector{Int64}(tau_ray_indices), length(tau_ray_indices)-1)
        tmp = Base.setdiff(indices_ith_max_cone, tau_ray_indices)
        tmp = Base.union(subset, tmp)
        push!(tmp, newindex)
        push!(newmaxcones, tmp)
      end
    else
      push!(newmaxcones, Vector{Int}(indices_ith_max_cone))
    end
  end
  newmaxcones = IncidenceMatrix(newmaxcones)
  
  # return the new fan
  return polyhedral_fan(T, newrays, newmaxcones; non_redundant=true)
  
end


function _cone_is_smooth(PF::_FanLikeType, c::AbstractSet{<:Integer})
  R = matrix(ZZ, rays(PF))
  return _is_unimodular(R[Vector{Int}(c),:])
end

function _is_unimodular(M::ZZMatrix)
  nrows(M) <= ncols(M) || return false
  n = nrows(M)
  return abs(det(snf(M)[:,1:n])) == 1
end


###############################################################################
## Cartesian/Direct product
###############################################################################

@doc raw"""
    *(PF1::PolyhedralFan, PF2::PolyhedralFan)

Return the Cartesian/direct product of two polyhedral fans.

# Examples
```jldoctest
julia> normal_fan(simplex(2))*normal_fan(simplex(3))
Polyhedral fan in ambient dimension 5
```
"""
function Base.:*(PF1::PolyhedralFan, PF2::PolyhedralFan)
    prod = Polymake.fan.product(pm_object(PF1), pm_object(PF2))
    return PolyhedralFan{detect_scalar_type(PolyhedralFan, prod)}(prod)
end
