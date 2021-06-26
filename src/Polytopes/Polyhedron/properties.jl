###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################
struct Polyhedra
end
struct Points
end
@doc Markdown.doc"""

    Halfspaces

Dummy type used for specifying the desired output format.
One halfspace `H(a,b)` is given by a vector `a` and a value `b` such that
$$H(a,b) = \{ x | ax ≤ b \}.$$
"""
struct Halfspaces
end

#TODO: take into account lineality space

struct PolyhedronFacePolyhedronIterator
    p::Polyhedron
    face_dim::Int
end


#Note: the following generates the entire Hasse diagram and may be costly
function Base.iterate(iter::PolyhedronFacePolyhedronIterator, index = 1)
    if iter.face_dim<0
        return nothing
    end
    faces = Polymake.polytope.faces_of_dim(pm_polytope(iter.p),iter.face_dim)
    nfaces = length(faces)
    while true
        if index > nfaces
            return nothing
        end
        isfar=true
        for v in faces[index]
            #Checking that the first coordinate is zero is equivalent
            #  to being a vertex of the far face
            if !iszero(pm_polytope(iter.p).VERTICES[1+v[1],1])
                isfar=false
            end
        end
        if isfar==true
            index +=1
        else
            p = Polyhedron(Polymake.polytope.Polytope(VERTICES=pm_polytope(iter.p).VERTICES[[f+1 for f in faces[index]],:],LINEALITY_SPACE = pm_polytope(iter.p).LINEALITY_SPACE))
            return (p,index+1)
        end
    end
end
#Note: it is impossible to know the number of faces prior to computation (and extraction of far faces)
Base.IteratorSize(::Type{PolyhedronFacePolyhedronIterator}) = Base.SizeUnknown()

"""
    faces(as, P, face_dim)

Return the faces of `P` of dimension `face_dim` as an iterator over the type of object given
by `as`.
Optional arguments for `as` include
* `Polyhedron`/`Polyhedra`: Return the representation of a vertex as a point.

# Arguments
- `as::Type{T}`: Element type of the returned iterator.
- `P::Polyhedron`: A polyhedron.
- `face_dim::Int`: Dimension of the desired faces.

# Example
An `Array` containing the six sides of the 3-dimensional cube can be obtained via the following input:
```julia-repl
julia> F = faces(Polyhedron, cube(3), 2)
Oscar.PolyhedronFacePolyhedronIterator(A polyhedron in ambient dimension 3, 2)

julia> collect(F)
6-element Vector{Any}:
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
```
"""
function faces(as::Type{T}, P::Polyhedron, face_dim::Int) where {T}
    if as == Polyhedron || as == Polyhedra
        PolyhedronFacePolyhedronIterator(P,face_dim-size(lineality_space(P),1))
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
end
"""
    faces(P, face_dim)

Return the faces of `P` of dimension `face_dim` as an iterator over `Polyhedron` objects.

# Arguments
- `P::Polyhedron`: A polyhedron.
- `face_dim::Int`: Dimension of the desired faces.

# Example
An `Array` containing the six sides of the 3-dimensional cube can be obtained via the following input:
```julia-repl
julia> F = faces(cube(3),2)
Oscar.PolyhedronFacePolyhedronIterator(A polyhedron in ambient dimension 3, 2)

julia> collect(F)
6-element Vector{Any}:
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
```
"""
faces(P::Polyhedron, face_dim::Int) = faces(Polyhedron, P, face_dim)

struct VertexPointIterator
    p::Polyhedron
end

function Base.iterate(iter::VertexPointIterator, index = 1)
    vertices = pm_polytope(iter.p).VERTICES
    while true
        if size(vertices, 1) < index
            return nothing
        end

        if iszero(vertices[index, 1])
            index += 1
        else
            return (vertices[index, 2:end], index + 1)
        end
    end
end
Base.eltype(::Type{VertexPointIterator}) = Polymake.Vector{Polymake.Rational}
Base.length(iter::VertexPointIterator) = nvertices(iter.p)

"""
    vertices(as, P)

Return an iterator over the vertices of a polyhedron `P` in the format defined by `as`.
Optional arguments for `as` include
* `Points`: Returns the representation of a vertex as a point.

See also `vertices_as_point_matrix`.

# Arguments
- `as::Type{T}`: Element type of the returned iterator.
- `P::Polyhedron`: A polyhedron.

# Examples
The following code computes the vertices of the Minkowski sum of a triangle and a square:
```julia-repl
julia> P = simplex(2) + cube(2);

julia> collect(vertices(Points, P))
5-element Vector{Polymake.Vector{Polymake.Rational}}:
pm::Vector<pm::Rational>
-1 -1
pm::Vector<pm::Rational>
2 -1
pm::Vector<pm::Rational>
2 1
pm::Vector<pm::Rational>
-1 2
pm::Vector<pm::Rational>
1 2
```
"""
function vertices(as::Type{T}, P::Polyhedron) where {T}
    if as == Points
        VertexPointIterator(P)
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
end

"""
    vertices(P)

Return an iterator over the vertices of a polyhedron `P` as points.

See also `vertices_as_point_matrix`.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
The following code computes the vertices of the Minkowski sum of a triangle and a square:
```julia-repl
julia> P = simplex(2) + cube(2);

julia> collect(vertices(P))
5-element Vector{Polymake.Vector{Polymake.Rational}}:
pm::Vector<pm::Rational>
-1 -1
pm::Vector<pm::Rational>
2 -1
pm::Vector<pm::Rational>
2 1
pm::Vector<pm::Rational>
-1 2
pm::Vector<pm::Rational>
1 2
```
"""
vertices(P::Polyhedron) = vertices(Points, P)

"""
    `vertices_as_point_matrix(P)`

Return a matrix whose rows are the vertices of `P`.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
The following code computes the vertices of the Minkowski sum of a triangle and a square:
```julia-repl
julia> P = simplex(2) + cube(2);

julia> vertices_as_point_matrix(P)
pm::Matrix<pm::Rational>
-1 -1
2 -1
2 1
-1 2
1 2
```
"""
function vertices_as_point_matrix(P::Polyhedron)
    decompose_vdata(pm_polytope(P).VERTICES).vertices
end


struct PolyhedronRayIterator
    p::Polyhedron
end

function Base.iterate(iter::PolyhedronRayIterator, index = 1)
    vertices = pm_polytope(iter.p).VERTICES
    while true
        if size(vertices, 1) < index
            return nothing
        end

        if !iszero(vertices[index, 1])
            index += 1
        else
            return (vertices[index, 2:end], index + 1)
        end
    end
end
Base.eltype(::Type{PolyhedronRayIterator}) = Polymake.Vector{Polymake.Rational}
Base.length(iter::PolyhedronRayIterator) = nrays(iter.p)

"""
    nrays(P)

Return the number of rays of a `Polyhedron`.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
Reflecting the input, the upper half-plane indeed has one ray.
```julia-repl
julia> UH = convex_hull([0 0],[0 1],[1 0]);

julia> nrays(UH)
1
```
"""
nrays(P::Polyhedron) = length(pm_polytope(P).FAR_FACE)

"""
    nvertices(P)

Return the number of vertices of the polyhedron `P`.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
The 3-cube's number of vertices can be obtained with this input:
```julia-repl
julia> C = cube(3);

julia> nvertices(C)
8
```
"""
nvertices(P::Polyhedron) = pm_polytope(P).N_VERTICES - nrays(P)


"""
    rays(as, P)


Return minimal set of generators of the cone of unbounded directions of a polyhedron `P` (i.e. its rays)
 in the format defined by `as`. Optional arguments for `as` include
* `Points`: Returns a vector representation of a ray.

See also `rays_as_point_matrix`.

# Arguments
- `as::Type{T}`: Element type of the returned iterator.
- `P::Polyhedron`: A polyhedron.

# Examples
We can verify that the positive orthant of the plane is generated by the two rays in positive unit direction:
```julia-repl
julia> PO = convex_hull([0 0], [1 0; 0 1]);
julia> collect(rays(Points, PO))
2-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 0
 pm::Vector<pm::Rational>
0 1
```
"""
function rays(as::Type{T}, P::Polyhedron) where {T}
    if as == Points
        PolyhedronRayIterator(P)
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
end

"""
    rays(P)


Return minimal set of generators of the cone of unbounded directions of a polyhedron `P` (i.e. its rays)
 as points.

See also `rays_as_point_matrix`.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
We can verify that the positive orthant of the plane is generated by the two rays in positive unit direction:
```julia-repl
julia> PO = convex_hull([0 0], [1 0; 0 1]);
julia> collect(rays(PO))
2-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 0
 pm::Vector<pm::Rational>
0 1
```
"""
rays(P::Polyhedron) = rays(Points,P)

function rays_as_point_matrix(P::Polyhedron)
    decompose_vdata(pm_polytope(P).VERTICES).rays
end



struct PolyhedronFacetHalfspaceIterator
    p::Polyhedron
end

function Base.iterate(iter::PolyhedronFacetHalfspaceIterator, index = 1)
    facets = pm_polytope(iter.p).FACETS
    if size(facets, 1) < index
        return nothing
    end

    return ((Polymake.Vector(-facets[index, 2:end]), facets[index, 1]), index + 1)
end
Base.length(iter::PolyhedronFacetHalfspaceIterator) = nfacets(iter.p)
Base.eltype(::Type{PolyhedronFacetHalfspaceIterator}) =
    Tuple{Polymake.Vector{Polymake.Rational},Polymake.Rational}


struct PolyhedronFacetPolyhedronIterator
    p::Polyhedron
end

function Base.iterate(iter::PolyhedronFacetPolyhedronIterator, index = 1)
    n_facets = nfacets(iter.p)
    if index > n_facets
        return nothing
    end

    p = Polyhedron(Polymake.polytope.facet(pm_polytope(iter.p), index - 1))
    return (p, index + 1)
end
Base.length(iter::PolyhedronFacetPolyhedronIterator) = nfacets(iter.p)
# Base.eltype(::Type{PolyhedronFacetPolyhedronIterator}) = Polyhedron

"""
    nfacets(P)

Return the number of facets of the polyhedron `P`.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
The number of facets of the 5-dimensional cross polytope can be retrieved via the following line:
```julia-repl
julia> nfacets(cross(5))
32
```
"""
nfacets(P::Polyhedron) = pm_polytope(P).N_FACETS

@doc Markdown.doc"""
    facets(as::Type{T}, P::Polyhedron)

Return the facets of the polyhedron `P` in the format defined by `as`.
The allowed values for `as` are
* `Halfspaces`: Returns for each facet the tuple `(A, b)` describing the halfspace `dot(A,x) ≤ b`.
* `Polyhedron` or `Polyhedra`: Returns for each facet its realization as a polyhedron

See also `facets_as_halfspace_matrix_pair`.

# Arguments
- `as::Type{T}`: Element type of the returned iterator.
- `P::Polyhedron`: A polyhedron.

# Examples
We can retrieve the six facets of the 3-dimensional cube this way:
```julia-repl
julia> C = cube(3);

julia> collect(facets(Polyhedron, C))
6-element Vector{Any}:
A polyhedron in ambient dimension 3
A polyhedron in ambient dimension 3
A polyhedron in ambient dimension 3
A polyhedron in ambient dimension 3
A polyhedron in ambient dimension 3
A polyhedron in ambient dimension 3

julia> collect(facets(Halfspaces, C))
6-element Vector{Tuple{Polymake.Vector{Polymake.Rational},Polymake.Rational}}:
(pm::Vector<pm::Rational>
-1 0 0, 1)
(pm::Vector<pm::Rational>
1 0 0, 1)
(pm::Vector<pm::Rational>
0 -1 0, 1)
(pm::Vector<pm::Rational>
0 1 0, 1)
(pm::Vector<pm::Rational>
0 0 -1, 1)
(pm::Vector<pm::Rational>
0 0 1, 1)
```
"""
function facets(as::Type{T}, P::Polyhedron) where {T}
    if as == Halfspaces
        PolyhedronFacetHalfspaceIterator(P)
    elseif as == Polyhedra || as == Polyhedron
        PolyhedronFacetPolyhedronIterator(P)
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
end
@doc Markdown.doc"""
    facets(P::Polyhedron)

Return the facets of the polyhedron `P` as halfspaces.

See also `facets_as_halfspace_matrix_pair`.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
We can retrieve the six facets of the 3-dimensional cube this way:
```julia-repl
julia> C = cube(3);

julia> collect(facets(C))
6-element Vector{Tuple{Polymake.Vector{Polymake.Rational},Polymake.Rational}}:
(pm::Vector<pm::Rational>
-1 0 0, 1)
(pm::Vector<pm::Rational>
1 0 0, 1)
(pm::Vector<pm::Rational>
0 -1 0, 1)
(pm::Vector<pm::Rational>
0 1 0, 1)
(pm::Vector<pm::Rational>
0 0 -1, 1)
(pm::Vector<pm::Rational>
0 0 1, 1)
```
"""
facets(P::Polyhedron) = facets(Halfspaces, P)

#TODO: how do underscores work in markdown?
@doc Markdown.doc"""

    facets_as_halfspace_matrix_pair(P::Polyhedron)

Return `(A,b)` such that $P=P(A,b)$ where

$$P(A,b) = \{ x |  Ax ≤ b \}.$$

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> C = cube(3);

julia> facets_as_halfspace_matrix_pair(C)
(A = pm::SparseMatrix<pm::Rational, pm::NonSymmetric>
(3) (0 -1)
(3) (0 1)
(3) (1 -1)
(3) (1 1)
(3) (2 -1)
(3) (2 1)
, b = pm::SparseVector<pm::Rational>
1 1 1 1 1 1)
```
"""
function facets_as_halfspace_matrix_pair(P::Polyhedron)
    return decompose_hdata(pm_polytope(P).FACETS)
end


###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################
"""
    volume(P)

Return the (Euclidean) volume of a polyhedron.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> C = cube(2);

julia> volume(C)
4
```
"""
volume(P::Polyhedron) = (pm_polytope(P)).VOLUME

"""
    normalized_volume(P)

Return the (normalized) volume of a polyhedron.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> C = cube(2);

julia> normalized_volume(C)
8
```
"""
normalized_volume(P::Polyhedron) = factorial(dim(P))*(pm_polytope(P)).VOLUME

"""
    dim(P)

Return the dimension of a polyhedron.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> P = convex_hull(V)
A polyhedron in ambient dimension 3

julia> dim(P)
2
```
"""
dim(P::Polyhedron) = Polymake.polytope.dim(pm_polytope(P))


"""
    lattice_points(P)

Return the integer points contained in a bounded polyhedron.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> S = 2 * simplex(2);

julia> lattice_points(S)
6-element Vector{Polymake.VectorAllocated{Polymake.Integer}}:
 pm::Vector<pm::Integer>
0 0
 pm::Vector<pm::Integer>
0 1
 pm::Vector<pm::Integer>
0 2
 pm::Vector<pm::Integer>
1 0
 pm::Vector<pm::Integer>
1 1
 pm::Vector<pm::Integer>
2 0
```
"""
function lattice_points(P::Polyhedron)
    if pm_polytope(P).BOUNDED
        lat_pts = pm_polytope(P).LATTICE_POINTS_GENERATORS[1]
        return ([e[2:end] for e in eachrow(lat_pts)])
    else
        throw(ArgumentError("Polyhedron not bounded"))
    end
end
#TODO: should this be an iterator too? If so, we should probably find a
#      scalable way to construct these iterators for so many functions

"""
    ambient_dim(P)

Return the ambient dimension of a polyhedron.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> P = convex_hull(V)
A polyhedron in ambient dimension 3

julia> ambient_dim(P)
3
```
"""
ambient_dim(P::Polyhedron) = Polymake.polytope.ambient_dim(pm_polytope(P))

"""
    codim(P)

Return the codimension of a polyhedron.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> P = convex_hull(V)
A polyhedron in ambient dimension 3

julia> codim(P)
1
```
"""
codim(P::Polyhedron) = ambient_dim(P)-dim(P)

###############################################################################
## Points properties
###############################################################################


# Previously: This implementation is not correct. Ask Taylor.
# Taylor: lineality space generators always look like [0, v] so
#  v is a natural output.
"""
    lineality_space(P)

Return a matrix whose row span is the lineality space of a polyhedron.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> UH = convex_hull([0 0],[0 1; 1 0; -1 0]);

julia> lineality_space(UH)
pm::Matrix<pm::Rational>
1 0
```
"""
lineality_space(P::Polyhedron) = dehomogenize(pm_polytope(P).LINEALITY_SPACE)


"""
    recession_cone(P)

Return the recession cone of `P`.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> P = Polyhedron([1 -2; -1 1; -1 0; 0 -1],[2,1,1,1]);

julia> collect(vertices(P))
3-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
0 -1
 pm::Vector<pm::Rational>
-1 0
 pm::Vector<pm::Rational>
-1 -1

julia> recession_cone(P)
A polyhedral cone in ambient dimension 2

julia> collect(rays(recession_cone(P)))
2-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 1/2
 pm::Vector<pm::Rational>
1 1
```
"""
recession_cone(P::Polyhedron) = Cone(Polymake.polytope.recession_cone(pm_polytope(P)))

###############################################################################
## Boolean properties
###############################################################################
"""
    isfeasible(P)

Check whether a polyhedron is feasible, i.e. non-empty.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> P = Polyhedron([1 -1; -1 1; -1 0; 0 -1],[-1,-1,1,1])
A polyhedron in ambient dimension 2

julia> isfeasible(P)
false
```
"""
isfeasible(P::Polyhedron) = pm_polytope(P).FEASIBLE


"""
    issmooth(P)

Check whether a polyhedron is smooth.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> C = cube(8);

julia> issmooth(C)
true
```
"""
issmooth(P::Polyhedron) = pm_polytope(P).SMOOTH


"""
    isnormal(P)

Check whether a polyhedron is normal.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
TODO
"""
isnormal(P::Polyhedron) = pm_polytope(P).NORMAL


"""
    isbounded(P)

Check whether a polyhedron is bounded.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> P = Polyhedron([1 -3; -1 1; -1 0; 0 -1],[1,1,1,1]);

julia> isbounded(P)
false
```
"""
isbounded(P::Polyhedron) = pm_polytope(P).BOUNDED


"""
    isfulldimensional(P)

Check whether a polyhedron is full dimensional.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> isfulldimensional(convex_hull(V))
false
```
"""
isfulldimensional(P::Polyhedron) = pm_polytope(P).FULL_DIM

"""
    f_vector(P)

Compute the vector`(f_1,f_2,...,f_(dim(P)-1))` where `f_i` is the number of faces of `P` of dimension `i`.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
Here we compute the f-vector of the 5-cube:
```julia-repl
julia> f_vector(cube(5))
5-element Vector{Int64}:
 32
 80
 80
 40
 10
```
"""
function f_vector(P::Polyhedron)
    f_vec=[length(collect(faces(P,i))) for i in 0:dim(P)-1]
    return f_vec
end


"""
    support_function(P; convention = :max)

Produce a function `h(ω) = max{dot(x,ω) | x ∈ P}`. max may be changed
to min by setting convention = :min.

# Arguments
- `P::Polyhedron`: A polyhedron.
- `convention::Symbol`: Convention to be applied.

# Examples
```julia-repl
julia> P = cube(3) + simplex(3);

julia> φ = support_function(P);

julia> φ([1,2,3])
9

julia> ψ = support_function(P, convention = :min);

julia> ψ([1,2,3])
-6
```
"""
function support_function(P::Polyhedron; convention = :max)
    function h(ω::AbstractVector)
        lp=LinearProgram(P,ω; convention = convention)
        return solve_lp(lp)[1]
    end
    return h
end
