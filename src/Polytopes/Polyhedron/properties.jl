###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################
struct Polyhedra
end
struct Points
end
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
    faces(P::Polyhedron, face_dim::Int [, as::Type{T} = Polyhedron])

Returns the faces of `P` of dimension `face_dim` as an iterator over the type of object given
by `as`.

# Arguments
- `P::Polyhedron`: A polyhedron.
- `face_dim::Int`: Dimension of the desired faces.
- `as::Type{T}`: Object type which is to be contained in the returned iterator.

# Example
An `Array` containing the six sides of the 3-dimensional cube can be obtained via the following input:
```julia-repl
julia> F = faces(cube(3),2)
Oscar.PolyhedronFacePolyhedronIterator(A polyhedron in ambient dimension 3, 2)

julia> collect(F)
6-element Array{Any,1}:
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
```
"""
function faces(P::Polyhedron, face_dim::Int, as::Type{T} = Polyhedron) where {T}
    if as == Polyhedron || as == Polyhedra
        PolyhedronFacePolyhedronIterator(P,face_dim-size(lineality_space(P),1))
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
end

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
   vertices(P::Polyhedron, [,as::Type{T} = Points])

Return an iterator over the vertices of a polyhedron `P` in the format defined by `as`.
Optional arguments for `as` include
* `Points`: Returns the representation of a vertex as a point.

See also `vertices_as_point_matrix`.

# Arguments
- `P::Polyhedron`: A polyhedron.
- `as::Type{T}`: Object type which is to be contained in the returned iterator.

# Examples
The following code computes the vertices of the Minkowski sum of a triangle and a square:
```julia-repl
julia> P = simplex(2) + cube(2);

julia> collect(vertices(P))
5-element Array{Polymake.Vector{Polymake.Rational},1}:
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
function vertices(P::Polyhedron, as::Type{T} = Points) where {T}
    if as == Points
        VertexPointIterator(P)
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
end


"""
   `vertices_as_point_matrix(P::Polyhedron)`

Returns a matrix whose rows are the vertices of `P`.

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
    nrays(P::Polyhedron)

Return the number of rays of a `Polyhedron`.

# Arguments
- `P::Polyhedron`: A polyhedron.

# Examples
```julia-repl
julia> UH = convex_hull([0 0],[0 1],[1 0]);

julia> nrays(UH)
1
```
"""
nrays(P::Polyhedron) = length(pm_polytope(P).FAR_FACE)

"""
    nvertices(P::Polyhedron)

Returns the number of vertices of `P`.
"""
nvertices(P::Polyhedron) = pm_polytope(P).N_VERTICES - nrays(P)


"""
   rays(P::Polyhedron, [,as::Type{T} = Points])


Returns minimal set of generators of the cone of unbounded directions of a polyhedron `P` (i.e. its rays)
 in the format defined by `as`. Optional arguments for `as` include
* `Points`: Returns a vector representation of a ray.

See also `rays_as_point_matrix`.
"""
function rays(P::Polyhedron, as::Type{T} = Points) where {T}
    if as == Points
        PolyhedronRayIterator(P)
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
end

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
   nfacets(P::Polyhedron)

Returns the number of facets of the polyhedron `P`.

# Examples
Compute number of facets of the 5-dimensional cross polytope
```julia-repl
julia> nfacets(cross(5))
32
```
"""
nfacets(P::Polyhedron) = pm_polytope(P).N_FACETS

@doc Markdown.doc"""
   facets(P::Polyhedron, as = :halfspaces)

Returns the facets of the polyhedron `P` in the format defined by `as`.
The allowed values for `as` are
* `Halfspaces`: Returns for each facet the tuple `(A, b)` describing the halfspace `dot(A,x) ≤ b`.
* `Polyhedron` or `Polyhedra`: Returns for each facet its realization as a polyhedron

See also `facets_as_halfspace_matrix_pair`.

# Examples
```julia-repl
julia> C=cube(3);

julia> collect(facets(C,Polyhedron))
6-element Array{Any,1}:
A polyhedron in ambient dimension 3
A polyhedron in ambient dimension 3
A polyhedron in ambient dimension 3
A polyhedron in ambient dimension 3
A polyhedron in ambient dimension 3
A polyhedron in ambient dimension 3

julia> collect(facets(C,Halfspaces))
6-element Array{Tuple{Polymake.Vector{Polymake.Rational},Polymake.Rational},1}:
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
function facets(P::Polyhedron,  as::Type{T} = Halfspaces) where {T}
    if as == Halfspaces
        PolyhedronFacetHalfspaceIterator(P)
    elseif as == Polyhedra || as == Polyhedron
        PolyhedronFacetPolyhedronIterator(P)
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
end

#TODO: how do underscores work in markdown?
@doc Markdown.doc"""

   `facets_as_halfspace_matrix_pair(P::Polyhedron)`

Returns `(A,b)` such that $P=P(A,b)$ where

$$P(A,b) = \{ x |  Ax ≤ b \}.$$
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
   volume(P::Polyhedron)

Returns the (Euclidean) volume of a polyhedron.
"""
volume(P::Polyhedron) = (pm_polytope(P)).VOLUME

"""
   normalized_volume(P::Polyhedron)

Returns the (normalized) volume of a polyhedron.
"""
normalized_volume(P::Polyhedron) = factorial(dim(P))*(pm_polytope(P)).VOLUME

"""
   dim(P::Polyhedron)

Returns the dimension of a polyhedron.
"""
dim(P::Polyhedron) = Polymake.polytope.dim(pm_polytope(P))


"""
   lattice_points(P::Polyhedron)

Returns the integer points contained in a bounded polyhedron.
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
   ambient_dim(P::Polyhedron)

Returns the ambient dimension of a polyhedron.
"""
ambient_dim(P::Polyhedron) = Polymake.polytope.ambient_dim(pm_polytope(P))

"""
   codim(P::Polyhedron)

Returns the codimension of a polyhedron.
"""
codim(P::Polyhedron) = ambient_dim(P)-dim(P)

###############################################################################
## Points properties
###############################################################################


# Previously: This implementation is not correct. Ask Taylor.
# Taylor: lineality space generators always look like [0, v] so
#  v is a natural output.
"""
   lineality_space(`P`::Polyhedron)

Returns a matrix whose row span is the lineality space of a `P`.
"""
lineality_space(P::Polyhedron) = dehomogenize(pm_polytope(P).LINEALITY_SPACE)


"""
    recession_cone(P::Polyhedron)

Returns the recession cone of `P`.
"""
recession_cone(P::Polyhedron) = Cone(Polymake.polytope.recession_cone(pm_polytope(P)))

###############################################################################
## Boolean properties
###############################################################################
"""
   isfeasible(P::Polyhedron)

   Check whether a polyhedron is feasible, i.e. non-empty.
"""
isfeasible(P::Polyhedron) = pm_polytope(P).FEASIBLE


"""
   issmooth(P::Polyhedron)

   Check whether a polyhedron is smooth.
"""
issmooth(P::Polyhedron) = pm_polytope(P).SMOOTH


"""
   isnormal(P::Polyhedron)

   Check whether a polyhedron is normal.
"""
isnormal(P::Polyhedron) = pm_polytope(P).NORMAL


"""
   isbounded(P::Polyhedron)

   Check whether a polyhedron is bounded.
"""
isbounded(P::Polyhedron) = pm_polytope(P).BOUNDED


"""
   isfulldimensional(P::Polyhedron)

   Check whether a polyhedron is full dimensional.
"""
isfulldimensional(P::Polyhedron) = pm_polytope(P).FULL_DIM

"""
    f_vector(P::Polyhedron)

Compute the vector`(f_1,f_2,...,f_(dim(P)-1))` where `f_i` is the number of faces of `P` of dimension `i`.

# Examples
Compute the f-vector of the 5-cube:
```julia-repl
julia> f_vector(cube(5))
5-element Array{Int64,1}:
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
    support_function(P::Polyhedron; convention = :max)

Produce a function `h(ω) = max{dot(x,ω) | x ∈ P}`. max may be changed
    to min by setting convention = :min.
"""
function support_function(P::Polyhedron; convention = :max)
    function h(ω::AbstractVector)
        lp=LinearProgram(P,ω; convention = convention)
        return solve_lp(lp)[1]
    end
    return h
end
