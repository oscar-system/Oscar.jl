###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

@doc raw"""
    rays(PF::PolyhedralFan)

Return the rays of `PF`.

# Examples
The rays of a normal fan of a cube point in every positive and negative unit
direction.
```jldoctest
julia> C = cube(3);

julia> NF = normal_fan(C);

julia> rays(NF)
6-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
 [-1, 0, 0]
 [0, 1, 0]
 [0, -1, 0]
 [0, 0, 1]
 [0, 0, -1]
```

As for the `Cone`, the rays may be converted to a matrix using the
`matrix(ring, ...)` function.
```jldoctest
julia> C = cube(3);

julia> NF = normal_fan(C);

julia> matrix(QQ, rays(NF))
[ 1    0    0]
[-1    0    0]
[ 0    1    0]
[ 0   -1    0]
[ 0    0    1]
[ 0    0   -1]
```
"""
rays(PF::_FanLikeType{T}) where T<:scalar_types = lineality_dim(PF) == 0 ? _rays(PF) : _empty_subobjectiterator(RayVector{T}, pm_object(PF))
_rays(PF::_FanLikeType{T}) where T<:scalar_types = SubObjectIterator{RayVector{T}}(pm_object(PF), _ray_fan, _nrays(PF))

_ray_fan(::Type{RayVector{T}}, PF::Polymake.BigObject, i::Base.Integer) where T<:scalar_types = RayVector{T}(view(PF.RAYS, i, :))

_vector_matrix(::Val{_ray_fan}, PF::Polymake.BigObject; homogenized=false) = homogenized ? homogenize(PF.RAYS, 0) : PF.RAYS

_matrix_for_polymake(::Val{_ray_fan}) = _vector_matrix

_maximal_cone(::Type{Cone{T}}, PF::Polymake.BigObject, i::Base.Integer) where T<:scalar_types = Cone{T}(Polymake.fan.cone(PF, i - 1))


@doc raw"""                                                 
    rays_modulo_lineality(as, F::PolyhedralFan)
                         
Return the rays of the polyhedral fan `F` up to lineality as a `NamedTuple`
with two iterators. If `F` has lineality `L`, then the iterator
`rays_modulo_lineality` iterates over representatives of the rays of `F/L`.
The iterator `lineality_basis` gives a basis of the lineality space `L`.

# Examples
```jldoctest
julia> P = convex_hull(QQFieldElem, [0 0; 1 0])
Polyhedron in ambient dimension 2

julia> NF = normal_fan(P)
Polyhedral fan in ambient dimension 2

julia> rmlF = rays_modulo_lineality(NF)
(rays_modulo_lineality = RayVector{QQFieldElem}[[1, 0], [-1, 0]], lineality_basis = RayVector{QQFieldElem}[[0, 1]])

julia> rmlF.rays_modulo_lineality
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
 [-1, 0]

julia> rmlF.lineality_basis
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [0, 1]

julia> rays(NF)
0-element SubObjectIterator{RayVector{QQFieldElem}}
```
"""
rays_modulo_lineality(F::_FanLikeType{T}) where T<:scalar_types = rays_modulo_lineality(NamedTuple{(:rays_modulo_lineality, :lineality_basis), Tuple{SubObjectIterator{RayVector{T}}, SubObjectIterator{RayVector{T}}}}, F) 
function rays_modulo_lineality(as::Type{NamedTuple{(:rays_modulo_lineality, :lineality_basis), Tuple{SubObjectIterator{RayVector{T}}, SubObjectIterator{RayVector{T}}}}}, F::_FanLikeType) where T<:scalar_types
    return (
        rays_modulo_lineality = _rays(F),
        lineality_basis = lineality_space(F)
    )
end
rays_modulo_lineality(as::Type{RayVector}, F::_FanLikeType) = _rays(F)
    


@doc raw"""
    maximal_cones(PF::PolyhedralFan)

Return the maximal cones of `PF`.

# Examples
Here we ask for the the number of rays for each maximal cone of the face fan of
the 3-cube and use that `maximal_cones` returns an iterator.
```jldoctest
julia> PF = face_fan(cube(3));

julia> for c in maximal_cones(PF)
       println(nrays(c))
       end
4
4
4
4
4
4
```
"""
maximal_cones(PF::_FanLikeType{T}) where T<:scalar_types = SubObjectIterator{Cone{T}}(pm_object(PF), _maximal_cone, n_maximal_cones(PF))

_ray_indices(::Val{_maximal_cone}, obj::Polymake.BigObject) = obj.MAXIMAL_CONES

@doc raw"""
    cones(PF::PolyhedralFan, cone_dim::Int)

Return an iterator over the cones of `PF` of dimension `cone_dim`.

# Examples
The 12 edges of the 3-cube correspond to the 2-dimensional cones of its face fan:
```jldoctest
julia> PF = face_fan(cube(3));

julia> cones(PF, 2)
12-element SubObjectIterator{Cone{QQFieldElem}}:
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
 Polyhedral cone in ambient dimension 3
```
"""
function cones(PF::_FanLikeType{T}, cone_dim::Int) where T<:scalar_types
    l = cone_dim - length(lineality_space(PF))
    l < 1 && return nothing
    return SubObjectIterator{Cone{T}}(pm_object(PF), _cone_of_dim, size(Polymake.fan.cones_of_dim(pm_object(PF), l), 1), (c_dim = l,))
end

function _cone_of_dim(::Type{Cone{T}}, PF::Polymake.BigObject, i::Base.Integer; c_dim::Int = 0) where T<:scalar_types
    return Cone{T}(Polymake.polytope.Cone{scalar_type_to_polymake[T]}(RAYS = PF.RAYS[collect(Polymake.row(Polymake.fan.cones_of_dim(PF, c_dim), i)),:], LINEALITY_SPACE = PF.LINEALITY_SPACE))
end

_ray_indices(::Val{_cone_of_dim}, PF::Polymake.BigObject; c_dim::Int = 0) = Polymake.fan.cones_of_dim(PF, c_dim)


@doc raw"""
    cones(PF::PolyhedralFan)

Return the ray indices of all non-zero-dimensional
cones in a polyhedral fan.

# Examples
```jldoctest
julia> PF = face_fan(cube(2))
Polyhedral fan in ambient dimension 2

julia> cones(PF)
8×4 IncidenceMatrix
[1, 3]
[2, 4]
[1, 2]
[3, 4]
[1]
[3]
[2]
[4]
```
"""
function cones(PF::_FanLikeType)
  pmo = pm_object(PF)
  ncones = pmo.HASSE_DIAGRAM.N_NODES
  cones = [Polymake._get_entry(pmo.HASSE_DIAGRAM.FACES,i) for i in 0:(ncones-1)];
  cones = filter(x->!(-1 in x) && length(x)>0, cones);
  cones = [Polymake.to_one_based_indexing(x) for x in cones];
  return IncidenceMatrix([Vector{Int}(x) for x in cones])
end


###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################

@doc raw"""
    dim(PF::PolyhedralFan)

Return the dimension of `PF`.

# Examples
This fan in the plane contains a 2-dimensional cone and is thus 2-dimensional
itself.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 -1], IncidenceMatrix([[1, 2], [3]]));

julia> dim(PF)
2
```
"""
dim(PF::_FanLikeType) = pm_object(PF).FAN_DIM::Int

@doc raw"""
    n_maximal_cones(PF::PolyhedralFan)

Return the number of maximal cones of `PF`.

# Examples
The cones given in this construction are non-redundant. Thus there are two
maximal cones.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 -1], IncidenceMatrix([[1, 2], [3]]));

julia> n_maximal_cones(PF)
2
```
"""
n_maximal_cones(PF::_FanLikeType) = pm_object(PF).N_MAXIMAL_CONES::Int

@doc raw"""
    n_cones(PF::PolyhedralFan)

Return the number of cones of `PF`.

# Examples
The cones given in this construction are non-redundant. There are six
cones in this fan.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 -1], IncidenceMatrix([[1, 2], [3]]))
Polyhedral fan in ambient dimension 2

julia> n_cones(PF)
4
```
"""
n_cones(PF::_FanLikeType) = nrows(cones(PF))

@doc raw"""
    ambient_dim(PF::PolyhedralFan)

Return the ambient dimension `PF`, which is the dimension of the embedding
space.

This is equal to the dimension of the fan if and only if the fan is
full-dimensional.

# Examples
The normal fan of the 4-cube is embedded in the same ambient space.
```jldoctest
julia> ambient_dim(normal_fan(cube(4)))
4
```
"""
ambient_dim(PF::_FanLikeType) = pm_object(PF).FAN_AMBIENT_DIM::Int

@doc raw"""
    nrays(PF::PolyhedralFan)

Return the number of rays of `PF`.

# Examples
The 3-cube has 8 vertices. Accordingly, its face fan has 8 rays.
```jldoctest
julia> nrays(face_fan(cube(3)))
8
```
"""
nrays(PF::_FanLikeType) = lineality_dim(PF) == 0 ? _nrays(PF) : 0
_nrays(PF::_FanLikeType) = pm_object(PF).N_RAYS::Int


@doc raw"""
    f_vector(PF::PolyhedralFan)

Compute the vector $(f₁,f₂,...,f_{dim(PF)-1})$` where $f_i$ is the number of
faces of $PF$ of dimension $i$.

# Examples
The f-vector of the normal fan of a polytope is the reverse of the f-vector of
the polytope.
```jldoctest
julia> c = cube(3)
Polyhedron in ambient dimension 3

julia> f_vector(c)
3-element Vector{ZZRingElem}:
 8
 12
 6


julia> nfc = normal_fan(c)
Polyhedral fan in ambient dimension 3

julia> f_vector(nfc)
3-element Vector{ZZRingElem}:
 6
 12
 8
```
"""
function f_vector(PF::_FanLikeType)
    pmf = pm_object(PF)
    ldim = pmf.LINEALITY_DIM
    return Vector{ZZRingElem}(vcat(fill(0,ldim),pmf.F_VECTOR))
end


@doc raw"""
    lineality_dim(PF::PolyhedralFan)

Return the dimension of the lineality space of the polyhedral fan `PF`, i.e.
the dimension of the largest linear subspace.

# Examples
The dimension of the lineality space is zero if and only if the fan is pointed.
```jldoctest
julia> C = convex_hull([0 0; 1 0])
Polyhedron in ambient dimension 2

julia> is_fulldimensional(C)
false

julia> nf = normal_fan(C)
Polyhedral fan in ambient dimension 2

julia> is_pointed(nf)
false

julia> lineality_dim(nf)
1
```
"""
lineality_dim(PF::_FanLikeType) = pm_object(PF).LINEALITY_DIM::Int

###############################################################################
## Points properties
###############################################################################

@doc raw"""
    lineality_space(PF::PolyhedralFan)

Return a non-redundant matrix whose rows are generators of the lineality space
of `PF`.

# Examples
This fan consists of two cones, one containing all the points with $y ≤ 0$ and
one containing all the points with $y ≥ 0$. The fan's lineality is the common
lineality of these two cones, i.e. in $x$-direction.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 0; 0 -1], IncidenceMatrix([[1, 2, 3], [3, 4, 1]]))
Polyhedral fan in ambient dimension 2

julia> lineality_space(PF)
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0]
```
"""
lineality_space(PF::_FanLikeType{T}) where T<:scalar_types = SubObjectIterator{RayVector{T}}(pm_object(PF), _lineality_fan, lineality_dim(PF))

_lineality_fan(::Type{RayVector{T}}, PF::Polymake.BigObject, i::Base.Integer) where T<:scalar_types = RayVector{T}(view(PF.LINEALITY_SPACE, i, :))

_generator_matrix(::Val{_lineality_fan}, PF::Polymake.BigObject; homogenized=false) = homogenized ? homogenize(PF.LINEALITY_SPACE, 0) : PF.LINEALITY_SPACE

_matrix_for_polymake(::Val{_lineality_fan}) = _generator_matrix

###############################################################################
## Boolean properties
###############################################################################
@doc raw"""
    is_pointed(PF::PolyhedralFan)

Determine whether `PF` is pointed, i.e. all its cones are pointed.

# Examples
The normal fan of a non-fulldimensional polytope is not pointed.
```jldoctest
julia> C = convex_hull([0 0; 1 0])
Polyhedron in ambient dimension 2

julia> is_fulldimensional(C)
false

julia> nf = normal_fan(C)
Polyhedral fan in ambient dimension 2

julia> is_pointed(nf)
false

julia> lineality_dim(nf)
1
```
"""
is_pointed(PF::_FanLikeType) = pm_object(PF).POINTED::Bool


@doc raw"""
    is_smooth(PF::PolyhedralFan{QQFieldElem})

Determine whether `PF` is smooth.

# Examples
Even though the cones of this fan cover the positive orthant together, one of
these und thus the whole fan is not smooth.
```jldoctest
julia> PF = PolyhedralFan([0 1; 2 1; 1 0], IncidenceMatrix([[1, 2], [2, 3]]));

julia> is_smooth(PF)
false
```
"""
is_smooth(PF::_FanLikeType{QQFieldElem}) = pm_object(PF).SMOOTH_FAN::Bool

@doc raw"""
    is_regular(PF::PolyhedralFan)

Determine whether `PF` is regular, i.e. the normal fan of a polytope.

# Examples
This fan is not complete and thus not regular.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 -1], IncidenceMatrix([[1, 2], [3]]));

julia> is_regular(PF)
false
```
"""
is_regular(PF::_FanLikeType) = pm_object(PF).REGULAR::Bool

@doc raw"""
    is_complete(PF::PolyhedralFan)

Determine whether `PF` is complete, i.e. its support, the set-theoretic union
of its cones, covers the whole space.

# Examples
Normal fans of polytopes are complete.
```jldoctest
julia> is_complete(normal_fan(cube(3)))
true
```
"""
is_complete(PF::_FanLikeType) = pm_object(PF).COMPLETE::Bool


@doc raw"""
    is_simplicial(PF::PolyhedralFan)

Determine whether `PF` is simplicial, i.e. every cone should be generated by a
basis of the ambient space.

# Examples
The `normal_fan` of the cube is simplicial, while the `face_fan` is not.
```jldoctest
julia> is_simplicial(normal_fan(cube(3)))
true

julia> is_simplicial(face_fan(cube(3)))
false
```
"""
is_simplicial(PF::_FanLikeType) = pm_object(PF).SIMPLICIAL::Bool


###############################################################################
## Primitive collections
###############################################################################

@doc raw"""
    primitive_collections(PF::PolyhedralFan)

Return the primitive collections of a polyhedral fan.

# Examples
```jldoctest
julia> primitive_collections(normal_fan(simplex(3)))
1-element Vector{Set{Int64}}:
 Set([4, 2, 3, 1])
```
"""
function primitive_collections(PF::_FanLikeType)
    @req is_simplicial(PF) "PolyhedralFan must be simplicial."
    I = ray_indices(maximal_cones(PF))
    K = SimplicialComplex(I)
    return minimal_nonfaces(K)
end


###############################################################################
## Star subdivision
###############################################################################

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
 [0, 0, 1]
 [1, 1, 1]
 [0, 1, 0]
 [-1, -1, -1]

julia> ray_indices(maximal_cones(star))
6×5 IncidenceMatrix
[1, 2, 3]
[2, 3, 4]
[1, 3, 4]
[2, 4, 5]
[1, 2, 5]
[1, 4, 5]
```
"""
function star_subdivision(PF::_FanLikeType{T}, n::Int) where T<:scalar_types
  
  # check if n-th cone exist
  @req n <= n_cones(PF) "Cannot subdivide cone $n as it does not exist"
  
  # check if n-th cone can be subdivided
  cone_list = cones(PF)
  nthcone = Polymake.row(cone_list, n)
  @req length(nthcone) > 1 "Cannot subdivide cone $n as it is generated by a single ray"
  
  # check if nth cone is smooth
  @req _cone_is_smooth(PF, nthcone) "Cannot subdivide maximal cone $n as it is not smooth"
  
  # compute new ray
  R = Polymake.common.primitive(Oscar.pm_object(PF).RAYS)
  newindex = size(R,1) + 1
  newray = sum([R[i,:] for i in nthcone])
  
  # compute new rays
  newrays = [R; transpose(newray)]
  
  # compute the new maximal cones
  max_cone_indices = ray_indices(maximal_cones(PF))
  cone_indices = Polymake.row(cone_list, n)
  newmaxcones = (Vector{Int})[]
  for i in 1:n_maximal_cones(PF)
    indices_ith_max_cone = Polymake.row(max_cone_indices, i)
    if issubset(cone_indices, indices_ith_max_cone)
      @req _cone_is_smooth(PF, indices_ith_max_cone) "All cones containing sigma need to be smooth"
      for subset in subsets(Vector{Int64}(cone_indices), length(cone_indices)-1)
        tmp = Base.setdiff(indices_ith_max_cone, cone_indices)
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
  return PolyhedralFan{T}(newrays, newmaxcones)
  
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
