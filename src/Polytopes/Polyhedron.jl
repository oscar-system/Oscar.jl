###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################


@doc Markdown.doc"""

    Polyhedron(A, b)

The (metric) polyhedron defined by

$$P(A,b) = \{ x |  Ax ≤ b \}.$$

see Def. 3.35 and Section 4.1. of Joswig, M. and Theobald, T. "Polyhedral and Algebraic Methods in Computational Geometry", Springer 2013.

# Arguments
- `A::Matrix`: Matrix corresponding to the linear coefficients of the inequalilites that describe P.
- `b::Vector`: Vector corresponding to the constant term of the inequalilites that describe P.

# Examples
```julia-repl
julia> A=[1 0; 0 1; -1 0 ; 0 -1];
julia> b=[1,1,0,0];
julia> Polyhedron([1 0;0 1;-1 0;0 -1],[1 1 0 0])
A polyhedron of dimension 2
```
""" struct Polyhedron #a real polymake polyhedron
    pm_polytope::Polymake.BigObject
    boundedness::Symbol # Values: :unknown, :bounded, :unbounded
end
function Polyhedron(pm_polytope::Polymake.BigObject)
    Polyhedron(pm_polytope, :unknown)
end
function Polyhedron(A::Union{Oscar.MatElem,AbstractMatrix}, b)
    Polyhedron(Polymake.polytope.Polytope{Polymake.Rational}(
        INEQUALITIES = matrix_for_polymake(remove_zero_rows([b -A])),
    ))
end

"""
    pm_polytope(P::Polyhedron)

Get the underlying polymake `Polytope`.
"""
pm_polytope(P::Polyhedron) = P.pm_polytope

==(P0::Polyhedron, P1::Polyhedron) = Polymake.polytope.equal_polyhedra(pm_polytope(P0), pm_polytope(P1))

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, P::Polyhedron)
    print(io, "A polyhedron in ambient dimension $(ambient_dim(P))")
end

Polymake.visual(P::Polyhedron; opts...) = Polymake.visual(pm_polytope(P); opts...)


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
by `as`. Optional arguments for `as` include
* `Polyhedron` or `Polyhedra`: Returns for each face its realization as a polyhedron
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

Returns an iterator over the vertices of a polyhedron `P` in the format defined by `as`.
Optional arguments for `as` include
* `Points`: Returns the representation of a vertex as a point.

See also `vertices_as_point_matrix`.
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

Returns the number of rays of `P`.
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
"""
nfacets(P::Polyhedron) = pm_polytope(P).N_FACETS

@doc Markdown.doc"""
   facets(P::Polyhedron, as = :halfspaces)

Returns the facets of the polyhedron `P` in the format defined by `as`.
The allowed values for `as` are
* `Halfspaces`: Returns for each facet the tuple `(A, b)` describing the halfspace `dot(A,x) ≤ b`.
* `Polyhedron` or `Polyhedra`: Returns for each facet its realization as a polyhedron

See also `facets_as_halfspace_matrix_pair`.
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

###############################################################################
###############################################################################
### Standard constructions
###############################################################################
###############################################################################

@doc Markdown.doc"""
   orbit_polytope(V::AbstractVecOrMat, G::PermGroup)

Construct the convex hull of the orbit of the point(s) in $V$ under the action of $G$.
"""
function orbit_polytope(V::AbstractMatrix, G::PermGroup)
   if size(V)[2] != degree(G)
      throw(ArgumentError("Dimension of points and group degree need to be the same."))
   end
   generators = PermGroup_to_polymake_array(G)
   pmGroup = Polymake.group.PermutationAction(GENERATORS=generators)
   pmPolytope = Polymake.polytope.orbit_polytope(homogenize(V,1), pmGroup)
   return Polyhedron(pmPolytope)
end
function orbit_polytope(V::AbstractVector, G::PermGroup)
   return orbit_polytope(Matrix(reshape(V,(1,length(V)))), G)
end

@doc Markdown.doc"""
   cube(d [, u, l])

Construct the $[-1,1]$-cube in dimension $d$. If $u$ and $l$ are given, the $[l,u]$-cube in dimension $d$ is returned.
""" cube(d) = Polyhedron(Polymake.polytope.cube(d))
cube(d, u, l) = Polyhedron(Polymake.polytope.cube(d, u, l))


@doc Markdown.doc"""
    convex_hull(V [, R [, L]])

The polytope given as the convex hull of the columns of V. Optionally, rays (R)
and generators of the lineality space (L) can be given as well.

see Def. 2.11 and Def. 3.1  of Joswig, M. and Theobald, T. "Polyhedral and Algebraic Methods in Computational Geometry", Springer 2013.
"""
function convex_hull(V::AnyVecOrMat; non_redundant::Bool=false)
    if !non_redundant
        pm_polytope =
            Polymake.polytope.Polytope{Polymake.Rational}(POINTS = matrix_for_polymake(homogenize(V, 1)))
        return Polyhedron(pm_polytope)
    else
        pm_polytope =
            Polymake.polytope.Polytope{Polymake.Rational}(VERTICES = matrix_for_polymake(homogenize(V, 1)))
            return Polyhedron(pm_polytope)
    end
end
#We want to be able to input trivial arguments. So we allow for R and L to be nothing
#  and we check that the length (# entries) of these matrices is positive, else we
#  call the function again with that argument = nothing
function convex_hull(V::AnyVecOrMat, R::Union{AnyVecOrMat,Nothing}; non_redundant::Bool=false)
    if R!=nothing && length(R) == 0
        return convex_hull(V,nothing;non_redundant=non_redundant)
    end
    points = stack(homogenize(V, 1), homogenize(R, 0))
    if !non_redundant
        pm_polytope = Polymake.polytope.Polytope{Polymake.Rational}(POINTS = matrix_for_polymake(points))
        return Polyhedron(pm_polytope)
    else
        pm_polytope = Polymake.polytope.Polytope{Polymake.Rational}(VERTICES = matrix_for_polymake(points))
        return Polyhedron(pm_polytope)
    end
end
function convex_hull(V::AnyVecOrMat, R::Union{AnyVecOrMat,Nothing}, L::AnyVecOrMat; non_redundant::Bool=false)
    if R!=nothing && length(R) == 0
        return convex_hull(V,nothing,L;non_redundant=non_redundant)
    end
    if L!=nothing && length(L) == 0
        return convex_hull(V,R,nothing; non_redundant=non_redundant)
    end
    points = stack(homogenize(V, 1), homogenize(R, 0))
    lineality = homogenize(L, 0)
    if !non_redundant
        pm_polytope = Polymake.polytope.Polytope{Polymake.Rational}(
            POINTS = matrix_for_polymake(points),
            INPUT_LINEALITY = matrix_for_polymake(lineality),
        )
        return Polyhedron(pm_polytope)
    else
        pm_polytope = Polymake.polytope.Polytope{Polymake.Rational}(
            VERTICES = matrix_for_polymake(points),
            LINEALITY_SPACE = matrix_for_polymake(lineality),
        )
        return Polyhedron(pm_polytope)
    end
end

"""
    newton_polytope(poly)

Compute the Newton polytope of the given polynomial `poly`.
"""
function newton_polytope(f)
    exponents = reduce(hcat, Oscar.exponent_vectors(f))'
    convex_hull(exponents)
end

function pm_far_face(P::Polyhedron)
    #We use cone instead of Polytope so the convex hull alg. doesn't give
    #  warning about no leading 1's
    Polymake.polytope.Cone(RAYS=pm_polytope(P).VERTICES[[a+1 for a in Array{Int,1}(pm_polytope(P).FAR_FACE)],2:end])
end


"""
    f_vector(P::Polyhedron)

Computes the vector`(f_1,f_2,...,f_(dim(P)-1))` where `f_i` is the number of faces of `P` of dimension `i`.
"""
function f_vector(P::Polyhedron)
    f_vec=[length(collect(faces(P,i))) for i in 0:dim(P)-1]
    return f_vec
end


"""
    support_function(P::Polyhedron; convention = :max)

Produces a function `h(ω) = max{dot(x,ω) | x ∈ P}`. max may be changed
    to min by setting convention = :min.
"""
function support_function(P::Polyhedron; convention = :max)
    function h(ω::AbstractVector)
        lp=LinearProgram(P,ω; convention = convention)
        return solve_lp(lp)[1]
    end
    return h
end


"""
   intersect(P::Polyhedron, Q::Polyhedron)

   Intersect two polyhedra.
"""
function intersect(P::Polyhedron, Q::Polyhedron)
   return Polyhedron(Polymake.polytope.intersection(pm_polytope(P), pm_polytope(Q)))
end


"""
   minkowski_sum(P::Polyhedron, Q::Polyhedron)

   Minkowski sum of two polyhedra.
"""
function minkowski_sum(P::Polyhedron, Q::Polyhedron; algorithm::Symbol=:standard)
   if algorithm == :standard
      return Polyhedron(Polymake.polytope.minkowski_sum(pm_polytope(P), pm_polytope(Q)))
   elseif algorithm == :fukuda
      return Polyhedron(Polymake.polytope.minkowski_sum_fukuda(pm_polytope(P), pm_polytope(Q)))
   else
      throw(ArgumentError("Unknown minkowski sum `algorithm` argument :" * string(algorithm)))
   end
end




#TODO: documentation  + extend to different fields.

"""
   +(P::Polyhedron, Q::Polyhedron)

   Minkowski sum of two polyhedra.
"""
+(P::Polyhedron, Q::Polyhedron) = minkowski_sum(P,Q)


#TODO: extend to different fields

"""
   *(k::Int, Q::Polyhedron)

   Returns the scaled polyhedron `kQ`.
"""
*(k::Int, P::Polyhedron) = Polyhedron(Polymake.polytope.scale(pm_polytope(P),k))


"""
   *(P::Polyhedron, k::Int)

   Returns the scaled polyhedron `kP`.
"""
*(P::Polyhedron,k::Int) = k*P


"""
   +(P::Polyhedron, v::AbstractVector)

   Returns the translation `P+v` of `P` by the vector `v`.
"""
function +(P::Polyhedron,v::AbstractVector)
    if ambient_dim(P) != length(v)
        throw(ArgumentError("Translation vector not correct dimension"))
    else
        return Polyhedron(Polymake.polytope.translate(pm_polytope(P),Polymake.Vector{Polymake.Rational}(v)))
    end
end


"""
   +(v::AbstractVector,P::Polyhedron)

   Returns the translation `v+P` of `P` by the vector `v`.
"""
+(v::AbstractVector,P::Polyhedron) = P+v

@doc Markdown.doc"""

   simplex(d[,n])

Construct the simplex which is the convex hull of the standard basis vectors
along with the origin in R^$d$, optionally scaled by $n$.
"""
simplex(d::Int64,n) = Polyhedron(Polymake.polytope.simplex(d,n))
simplex(d::Int64) = Polyhedron(Polymake.polytope.simplex(d))


@doc Markdown.doc"""

   cross(d[,n])

Construct a $d$-dimensional cross polytope around origin with vertices located at $\pm e_i$ for each unit vector $e_i$ of $R^d$.
If $n$ is not given, construct the unit cross polytope around origin.
"""
cross(d::Int64,n) = Polyhedron(Polymake.polytope.cross(d,n))
cross(d::Int64) = Polyhedron(Polymake.polytope.cross(d))
