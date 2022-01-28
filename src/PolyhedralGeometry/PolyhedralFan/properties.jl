###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

@doc Markdown.doc"""
    rays(PF::PolyhedralFan)

Return the rays of `PF`.

# Examples
The rays of a normal fan of a cube point in every positive and negative unit
direction.
```jldoctest
julia> C = cube(3);

julia> NF = normal_fan(C);

julia> rays(NF)
6-element SubObjectIterator{RayVector{Polymake.Rational}}:
 [1, 0, 0]
 [-1, 0, 0]
 [0, 1, 0]
 [0, -1, 0]
 [0, 0, 1]
 [0, 0, -1]
```
"""
rays(PF::PolyhedralFan) = SubObjectIterator{RayVector{Polymake.Rational}}(pm_object(PF), _ray_fan, pm_object(PF).N_RAYS)

_ray_fan(::Type{RayVector{Polymake.Rational}}, PF::Polymake.BigObject, i::Base.Integer) = RayVector(PF.RAYS[i, :])

_vector_matrix(::Val{_ray_fan}, PF::Polymake.BigObject; homogenized=false) = PF.RAYS

_matrix_for_polymake(::Val{_ray_fan}) = _vector_matrix

_maximal_cone(::Type{Cone}, PF::Polymake.BigObject, i::Base.Integer) = Cone(Polymake.fan.cone(PF, i - 1))

@doc Markdown.doc"""
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
maximal_cones(PF::PolyhedralFan) = SubObjectIterator{Cone}(pm_object(PF), _maximal_cone, n_maximal_cones(PF))

_ray_indices(::Val{_maximal_cone}, obj::Polymake.BigObject) = obj.MAXIMAL_CONES

@doc Markdown.doc"""
    cones(PF::PolyhedralFan, cone_dim::Int)

Return an iterator over the cones of `PF` of dimension `cone_dim`.

# Examples
The 12 edges of the 3-cube correspond to the 2-dimensional cones of its face fan:
```jldoctest
julia> PF = face_fan(cube(3));

julia> cones(PF, 2)
12-element SubObjectIterator{Cone}:
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
 A polyhedral cone in ambient dimension 3
```
"""
function cones(PF::PolyhedralFan, cone_dim::Int)
    l = cone_dim  - length(lineality_space(PF))
    l < 1 && return nothing
    return SubObjectIterator{Cone}(pm_object(PF), _cone_of_dim, size(Polymake.fan.cones_of_dim(pm_object(PF), l), 1), (c_dim = l,))
end

function _cone_of_dim(::Type{Cone}, PF::Polymake.BigObject, i::Base.Integer; c_dim::Int = 0)
    return Cone(Polymake.polytope.Cone(RAYS = PF.RAYS[collect(Polymake.row(Polymake.fan.cones_of_dim(PF, c_dim), i)),:], LINEALITY_SPACE = PF.LINEALITY_SPACE))
end

_ray_indices(::Val{_cone_of_dim}, PF::Polymake.BigObject; c_dim::Int = 0) = Polymake.fan.cones_of_dim(PF, c_dim)

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################

@doc Markdown.doc"""
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
dim(PF::PolyhedralFan) = pm_object(PF).FAN_DIM

@doc Markdown.doc"""
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
n_maximal_cones(PF::PolyhedralFan) = pm_object(PF).N_MAXIMAL_CONES

@doc Markdown.doc"""
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
ambient_dim(PF::PolyhedralFan) = pm_object(PF).FAN_AMBIENT_DIM

@doc Markdown.doc"""
    nrays(PF::PolyhedralFan)

Return the number of rays of `PF`.

# Examples
The 3-cube has 8 vertices. Accordingly, its face fan has 8 rays.
```jldoctest
julia> nrays(face_fan(cube(3)))
8
```
"""
nrays(PF::PolyhedralFan) = pm_object(PF).N_RAYS


@doc Markdown.doc"""
    f_vector(PF::PolyhedralFan)

Compute the vector $(f₁,f₂,...,f_{dim(PF)-1})$` where $f_i$ is the number of
faces of $PF$ of dimension $i$.

# Examples
The f-vector of the normal fan of a polytope is the reverse of the f-vector of
the polytope.
```jldoctest
julia> c = cube(3)
A polyhedron in ambient dimension 3

julia> f_vector(c)
3-element Vector{Int64}:
  8
 12
  6

julia> nfc = normal_fan(c)
A polyhedral fan in ambient dimension 3

julia> f_vector(nfc)
3-element Vector{Polymake.Integer}:
  6
 12
  8
```
"""
function f_vector(PF::PolyhedralFan)
    pmf = pm_object(PF)
    ldim = pmf.LINEALITY_DIM
    return vcat(fill(0,ldim),pmf.F_VECTOR)
end


@doc Markdown.doc"""
    lineality_dim(PF::PolyhedralFan)

Return the dimension of the lineality space of the polyhedral fan `PF`, i.e.
the dimension of the largest linear subspace.

# Examples
The dimension of the lineality space is zero if and only if the fan is pointed.
```jldoctest
julia> C = convex_hull([0 0; 1 0])
A polyhedron in ambient dimension 2

julia> isfulldimensional(C)
false

julia> nf = normal_fan(C)
A polyhedral fan in ambient dimension 2

julia> ispointed(nf)
false

julia> lineality_dim(nf)
1
```
"""
lineality_dim(PF::PolyhedralFan) = pm_object(PF).LINEALITY_DIM

###############################################################################
## Points properties
###############################################################################

@doc Markdown.doc"""
    lineality_space(PF::PolyhedralFan)

Return a non-redundant matrix whose rows are generators of the lineality space
of `PF`.

# Examples
This fan consists of two cones, one containing all the points with $y ≤ 0$ and
one containing all the points with $y ≥ 0$. The fan's lineality is the common
lineality of these two cones, i.e. in $x$-direction.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 0; 0 -1], IncidenceMatrix([[1, 2, 3], [3, 4, 1]]))
A polyhedral fan in ambient dimension 2

julia> lineality_space(PF)
1-element SubObjectIterator{RayVector{Polymake.Rational}}:
 [1, 0]
```
"""
lineality_space(PF::PolyhedralFan) = SubObjectIterator{RayVector{Polymake.Rational}}(pm_object(PF), _lineality_fan, lineality_dim(PF))

_lineality_fan(::Type{RayVector{Polymake.Rational}}, PF::Polymake.BigObject, i::Base.Integer) = RayVector(PF.LINEALITY_SPACE[i, :])

_generator_matrix(::Val{_lineality_fan}, PF::Polymake.BigObject; homogenized=false) = PF.LINEALITY_SPACE

_matrix_for_polymake(::Val{_lineality_fan}) = _generator_matrix

###############################################################################
## Boolean properties
###############################################################################
@doc Markdown.doc"""
    ispointed(PF::PolyhedralFan)

Determine whether `PF` is pointed, i.e. all its cones are pointed.

# Examples
The normal fan of a non-fulldimensional polytope is not pointed.
```jldoctest
julia> C = convex_hull([0 0; 1 0])
A polyhedron in ambient dimension 2

julia> isfulldimensional(C)
false

julia> nf = normal_fan(C)
A polyhedral fan in ambient dimension 2

julia> ispointed(nf)
false

julia> lineality_dim(nf)
1
```
"""
ispointed(PF::PolyhedralFan) = pm_object(PF).POINTED


@doc Markdown.doc"""
    issmooth(PF::PolyhedralFan)

Determine whether `PF` is smooth.

# Examples
Even though the cones of this fan cover the positive orthant together, one of
these und thus the whole fan is not smooth.
```jldoctest
julia> PF = PolyhedralFan([0 1; 2 1; 1 0], IncidenceMatrix([[1, 2], [2, 3]]));

julia> issmooth(PF)
false
```
"""
issmooth(PF::PolyhedralFan) = pm_object(PF).SMOOTH_FAN

@doc Markdown.doc"""
    isregular(PF::PolyhedralFan)

Determine whether `PF` is regular, i.e. the normal fan of a polytope.

# Examples
This fan is not complete and thus not regular.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 -1], IncidenceMatrix([[1, 2], [3]]));

julia> isregular(PF)
false
```
"""
isregular(PF::PolyhedralFan) = pm_object(PF).REGULAR

@doc Markdown.doc"""
    iscomplete(PF::PolyhedralFan)

Determine whether `PF` is complete, i.e. its support, the set-theoretic union
of its cones, covers the whole space.

# Examples
Normal fans of polytopes are complete.
```jldoctest
julia> iscomplete(normal_fan(cube(3)))
true
```
"""
iscomplete(PF::PolyhedralFan) = pm_object(PF).COMPLETE


###############################################################################
## Primitive collections
###############################################################################

@doc Markdown.doc"""
    primitive_collections(PF::PolyhedralFan)

Return the primitive collections of a polyhedral fan.

# Examples
```jldoctest
julia> primitive_collections(normal_fan(Oscar.simplex(3)))
1-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
```
"""
function primitive_collections(PF::PolyhedralFan)
    # collect data
    cones = [findall(x->x!=0, l) for l in eachrow(pm_object(PF).MAXIMAL_CONES)]
    all_points = [i for i in 1 : pm_object(PF).N_RAYS]
    d_max = maximum([length(i) for i in cones]) + 1
    # identify and return the primitive collections
    collections = Vector{Int}[]
    for d in 1:d_max
        checked  = Vector{Int}[]
        for cone in cones
            d <= length(cone) || continue
            for I_minus_j in Oscar.Hecke.subsets(cone, d)
                scanner = setdiff(all_points, I_minus_j)
                for j in scanner
                    I = vcat(I_minus_j, [j])
                    I in checked && continue
                    push!(checked, I)
                    # (1) I is contained in the primitive collections iff it is not contained in any cone
                    if !any(test_cone -> issubset(I, test_cone), cones)
                        # (2) I is generator of the primitive collections iff primitive_collections does not contain a "smaller" generator
                        if !any(prim -> issubset(prim, I), collections)
                            push!(collections, I) # add new generator
                        end
                    end
                end
            end
        end
    end
    # return the computed primitive collections
    return collections
end


###############################################################################
## Star subdivision
###############################################################################

@doc Markdown.doc"""
    starsubdivision(PF::PolyhedralFan, n::Int)

Return the star subdivision of a polyhedral fan at its n-th maximal torus orbit.

# Examples
```jldoctest
julia> star = starsubdivision(normal_fan(simplex(3)), 1)
A polyhedral fan in ambient dimension 3

julia> rays(star)
5-element SubObjectIterator{RayVector{Polymake.Rational}}:
 [0, 1, 0]
 [0, 0, 1]
 [-1, -1, -1]
 [1, 0, 0]
 [1, 1, 1]

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
function starsubdivision(PF::PolyhedralFan, n::Int)
    # extract defining information on the fan
    maxcones = IncidenceMatrix(pm_object(PF).MAXIMAL_CONES)
    R = Polymake.common.primitive(pm_object(PF).RAYS)
    
    # check if n-th maximal cone does exist
    if length(maxcones) < n
        throw(ArgumentError("Cannot subdivide maximal cone $n as it does not exist!"))
    end
    nthmaxcone = Polymake.row(maxcones, n)
    
    # construct this cone and check if it is smooth
    cone = Polymake.fan.cone(pm_object(PF), n-1)
    if !cone.SMOOTH_CONE
        throw(ArgumentError("Cannot subdivide maximal cone $n as it is not smooth!"))
    end
    
    # compute new rays to be added from star subdivision
    newindex = size(R,1) + 1
    newray = sum([R[i,:] for i in nthmaxcone])
    
    # add this ray to form list of the new rays
    newrays = [R; transpose(newray)]
    
    # identify all maximal cones in the new fan
    d = Polymake.polytope.dim(cone)
    newmaxcones = [Vector{Int64}(Polymake.row(maxcones, i)) for i in 1:(Polymake.nrows(maxcones)) if i!= n]
    for subset in Oscar.Hecke.subsets(Vector{Int64}(nthmaxcone), d-1)
        tmp = Vector{Int64}(subset)
        append!(tmp, newindex)
        push!(newmaxcones, tmp)
    end
    newmaxcones = IncidenceMatrix(newmaxcones)
    
    # return the new fan
    return PolyhedralFan(newrays, newmaxcones)
    
end

###############################################################################
## Cartesian/Direct product
###############################################################################

@doc Markdown.doc"""
    *(PF1::PolyhedralFan, PF2::PolyhedralFan)

Return the Cartesian/direct product of two polyhedral fans.

# Examples
```jldoctest
julia> normal_fan(Oscar.simplex(2))*normal_fan(Oscar.simplex(3))
A polyhedral fan in ambient dimension 5
```
"""
function Base.:*(PF1::PolyhedralFan, PF2::PolyhedralFan)
    return PolyhedralFan(Polymake.fan.product(pm_object(PF1), pm_object(PF2)))
end
