@doc Markdown.doc"""
    Polyhedron(A, b)

The (metric) polyhedron defined by

$$P(A,b) = \{ x |  Ax ≤ b \}.$$

see Def. 3.35 and Section 4.1.
""" struct Polyhedron #a real polymake polyhedron
    pm_polytope::Polymake.BigObjectAllocated
    boundedness::Symbol # Values: :unknown, :bounded, :unbounded
end
function Polyhedron(pm_polytope::Polymake.BigObjectAllocated)
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

property_is_computed(P::Polyhedron, S::Symbol) = property_is_computed(pm_polytope(P), S)
==(P0::Polyhedron, P1::Polyhedron) = Polymake.polytope.equal_polyhedra(pm_polytope(P0), pm_polytope(P1))

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, P::Polyhedron)
    print(io, "A Polyhedron of dimension $(dim(P))")
end



###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################


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
Base.eltype(::Type{VertexPointIterator}) = Polymake.VectorAllocated{Polymake.Rational}
Base.length(iter::VertexPointIterator) = n_vertices(iter.p)

"""
   vertices(H::Polyhedron, as = :points)

Returns the vertices of a polyhedron.
"""
function vertices(P::Polyhedron; as::Symbol = :points)
    if as == :points
        VertexPointIterator(P)
    elseif as == :point_matrix
        return decompose_vdata(pm_polytope(P).VERTICES).vertices
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
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
Base.eltype(::Type{PolyhedronRayIterator}) = Polymake.VectorAllocated{Polymake.Rational}
Base.length(iter::PolyhedronRayIterator) = n_rays(iter.p)

"""
    n_rays(P::Polyhedron)

Returns the number of rays of `P`.
"""
n_rays(P::Polyhedron) = length(pm_polytope(P).FAR_FACE)

"""
    n_vertices(P::Polyhedron)

Returns the number of vertices of `P`.
"""
n_vertices(P::Polyhedron) = pm_polytope(P).N_VERTICES - n_rays(P)

"""
   rays(P::Polyhedron)

Returns minimal set of generators of the cone of unbounded directions of a polyhedron.
"""
function rays(P::Polyhedron; as::Symbol = :points)
    if as == :points
        PolyhedronRayIterator(P)
    elseif as == :point_matrix
        return decompose_vdata(pm_polytope(P).VERTICES).rays
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
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
Base.length(iter::PolyhedronFacetHalfspaceIterator) = n_facets(iter.p)
Base.eltype(::Type{PolyhedronFacetHalfspaceIterator}) =
    Tuple{Polymake.VectorAllocated{Polymake.Rational},Polymake.RationalAllocated}


struct PolyhedronFacetPolyhedronIterator
    p::Polyhedron
end

function Base.iterate(iter::PolyhedronFacetPolyhedronIterator, index = 1)
    nfacets = n_facets(iter.p)
    if index > nfacets
        return nothing
    end

    p = Polyhedron(Polymake.polytope.facet(pm_polytope(iter.p), index - 1))
    return (p, index + 1)
end
Base.length(iter::PolyhedronFacetPolyhedronIterator) = n_facets(iter.p)
# Base.eltype(::Type{PolyhedronFacetPolyhedronIterator}) = Polyhedron

"""
   n_facets(P::Polyhedron)

Returns the number of facets of the polyhedron `P`.
"""
n_facets(P::Polyhedron) = pm_polytope(P).N_FACETS

@doc Markdown.doc"""
   facets(P::Polyhedron, as = :halfspaces)

Returns the facets of the polyhedron `P` in the format defined by `as`.
The allowed values for `as` are
* `halfspaces`: Returns for each facet the tuple `(A, b)` describing the halfspace `dot(A,x) ≤ b`.
* `polyhedra`: Returns for each facet its realization as a polyhedron
* `halfspace_matrix_pair`: Returns `(A,b)` such `P={x | Ax ≦ b }`
"""
function facets(P::Polyhedron; as::Symbol = :halfspaces)
    if as == :halfspaces
        PolyhedronFacetHalfspaceIterator(P)
    elseif as == :polyhedra
        PolyhedronFacetPolyhedronIterator(P)
    elseif as == :halfspace_matrix_pair
        return decompose_hdata(pm_polytope(P).FACETS)
    else
        throw(ArgumentError("Unsupported `as` argument :" * string(as)))
    end
end



###############################################################################
###############################################################################
### Access properties
###############################################################################
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
   ambient_dim(P::Polyhedron)

Returns the ambient dimension of a polyhedron.
"""
ambient_dim(P::Polyhedron) = Polymake.polytope.ambient_dim(pm_polytope(P))

"""
   codim(P::Polyhedron)

Returns the codimension of a polyhedron.
"""
codim(P::Polyhedron) = ambient_dim(P)-dim(P)



@doc Markdown.doc"""
   cube(d [, u, l])

Construct the $[-1,1]$-cube in dimension $d$. If $u$ and $l$ are given, the $[l,u]$-cube in dimension $d$ is returned.
""" cube(d) = Polyhedron(Polymake.polytope.cube(d))
cube(d, u, l) = Polyhedron(Polymake.polytope.cube(d, u, l))

# TODO: This implementation is not correct. Ask Taylor.
# Taylor: lineality space generators always look like [0, v] so
#  v is a natural output.
"""
   lineality_space(H::Polyhedron)

Returns a basis of the lineality space of a polyhedron.
"""
lineality_space(H::Polyhedron) = dehomogenize(pm_polytope(H).LINEALITY_SPACE)

Polymake.visual(P::Polyhedron; opts...) = Polymake.visual(pm_polytope(P); opts...)

"""
    recession_cone(P::Polyhedron)

Returns the recession cone of `P`.
"""
recession_cone(P::Polyhedron) = Cone(Polymake.polytope.recession_cone(pm_polytope(P)))

###############################################################################
###############################################################################
### Standard constructions
###############################################################################
###############################################################################

@doc Markdown.doc"""
    convex_hull(V [, R [, L]])

The polytope given as the convex hull of the columns of V. Optionally, rays (R)
and generators of the lineality space (L) can be given as well.

see Def. 2.11 and Def. 3.1.
""" function convex_hull(V::AbstractVecOrMat)
    pm_polytope =
        Polymake.polytope.Polytope{Polymake.Rational}(POINTS = matrix_for_polymake(homogenize(V, 1)))
    return Polyhedron(pm_polytope)
end
function convex_hull(V::AbstractVecOrMat, R::AbstractVecOrMat)
    points = stack(homogenize(V, 1), homogenize(R, 0))
    pm_polytope = Polymake.polytope.Polytope{Polymake.Rational}(POINTS = matrix_for_polymake(points))
    return Polyhedron(pm_polytope)
end
function convex_hull(V::AbstractVecOrMat, R::AbstractVecOrMat, L::AbstractVecOrMat)
    points = stack(homogenize(V, 1), homogenize(R, 0))
    lineality = homogenize(L, 0)
    pm_polytope = Polymake.polytope.Polytope{Polymake.Rational}(
        POINTS = matrix_for_polymake(points),
        INPUT_LINEALITY = matrix_for_polymake(lineality),
    )
    return Polyhedron(pm_polytope)
end

"""
    newton_polytope(poly)

Compute the Newton polytope of the given polynomial `poly`.
"""
function newton_polytope(f)
    exponents = reduce(hcat, Oscar.exponent_vectors(f))'
    convex_hull(exponents)
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
