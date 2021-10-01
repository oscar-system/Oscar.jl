###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

@doc Markdown.doc"""

    Polyhedron(P::Polymake.BigObject)

Construct a `Polyhedron` corresponding to a `Polymake.BigObject` of type `Polytope`.
"""
function Polyhedron(pm_polytope::Polymake.BigObject)
    Polyhedron(pm_polytope, :unknown)
end

@doc Markdown.doc"""

    Polyhedron(A::Union{Oscar.MatElem,AbstractMatrix}, b)

The (convex) polyhedron defined by

$$P(A,b) = \{ x |  Ax ≤ b \}.$$

see Def. 3.35 and Section 4.1. of [JT13](@cite)

# Examples
The following lines define the square $[0,1]^2 \subset \mathbb{R}^2$:
```jldoctest
julia> A = [1 0; 0 1; -1 0 ; 0 -1];

julia> b = [1, 1, 0, 0];

julia> Polyhedron(A,b)
A polyhedron in ambient dimension 2
```
"""
Polyhedron(A::Union{Oscar.MatElem,AbstractMatrix}, b) = Polyhedron((A, b))

function Polyhedron(I::Union{HalfspaceIterator, Tuple{<:Union{Oscar.MatElem, AbstractMatrix}, Any}}, E::Union{Nothing, HalfspaceIterator, Tuple{<:Union{Oscar.MatElem, AbstractMatrix}, Any}} = nothing; non_redundant::Bool = false)
    IM = -affine_matrix_for_polymake(I)
    EM = isnothing(E) || _isempty_halfspace(E) ? Polymake.Matrix{Polymake.Rational}(undef, 0, size(IM, 2)) : affine_matrix_for_polymake(E)

    if non_redundant
        return Polyhedron(Polymake.polytope.Polytope{Polymake.Rational}(FACETS = IM, AFFINE_HULL = EM))
    else
        return Polyhedron(Polymake.polytope.Polytope{Polymake.Rational}(INEQUALITIES = IM, EQUATIONS = EM))
    end
end

"""
    pm_polytope(P::Polyhedron)

Get the underlying polymake `Polytope`.
"""
pm_polytope(P::Polyhedron) = P.pm_polytope

==(P0::Polyhedron, P1::Polyhedron) = Polymake.polytope.equal_polyhedra(pm_polytope(P0), pm_polytope(P1))


### Construct polyhedron from V-data, as the convex hull of points, rays and lineality.
@doc Markdown.doc"""
    convex_hull(V::Matrix [, R::Matrix [, L::Matrix]]; non_redundant::Bool = false)

Construct the convex hull of the vertices `V`, rays `R`, and lineality `L`. If
`R` or `L` are omitted, then they are assumed to be zero.

# Arguments
- `V::Matrix`: Points whose convex hull is to be computed; encoded as row vectors.
- `R::Matrix`: Rays completing the set of points; encoded row-wise as representative vectors.
- `L::Matrix`: Generators of the Lineality space; encoded as row vectors.

`R` can be given as an empty matrix or as `nothing` if the user wants to compute
the convex hull only from `V` and `L`.

If it is known that `V` and `R` only contain extremal points and that the
description of the lineality space is complete, set `non_redundant =
true` to avoid unnecessary redundancy checks.

See Def. 2.11 and Def. 3.1  of [JT13](@cite).

# Examples
The following lines define the square $[0,1]^2 \subset \mathbb{R}^2$:
```jldoctest
julia> Square = convex_hull([0 0; 0 1; 1 0; 1 1])
A polyhedron in ambient dimension 2
```
To construct the positive orthant, rays have to be passed:
```jldoctest
julia> V = [0 0];

julia> R = [1 0; 0 1];

julia> PO = convex_hull(V, R)
A polyhedron in ambient dimension 2
```
The closed-upper half plane can be constructed by passing rays and a lineality space:
```jldoctest
julia> V = [0 0];

julia> R = [0 1];

julia> L = [1 0];

julia> UH = convex_hull(V, R, L)
A polyhedron in ambient dimension 2
```
To obtain the x-axis in $\mathbb{R}^2$:
```jldoctest
julia> V = [0 0];

julia> R = nothing;

julia> L = [1 0];

julia> XA = convex_hull(V, R, L)
A polyhedron in ambient dimension 2
```
"""
function convex_hull(V::Union{VectorIterator{PointVector}, AnyVecOrMat, Oscar.MatElem}, R::Union{VectorIterator{RayVector}, AnyVecOrMat, Oscar.MatElem, Nothing} = nothing, L::Union{VectorIterator{RayVector}, AnyVecOrMat, Oscar.MatElem, Nothing} = nothing; non_redundant::Bool = false)
    # we access the matrices which polymake can work with.
    VM = matrix_for_polymake(V)
    RM = isnothing(R) || isempty(R) ? Polymake.Matrix{Polymake.Rational}(undef, 0, size(VM, 2)) : matrix_for_polymake(R)
    LM = isnothing(L) || isempty(L) ? Polymake.Matrix{Polymake.Rational}(undef, 0, size(VM, 2)) : matrix_for_polymake(L)

    # Rays and Points are homogenized and combined and
    # Lineality is homogenized
    points = stack(homogenize(VM, 1), homogenize(RM, 0))
    lineality = homogenize(LM, 0)

    # These matrices are in the right format for polymake.
    # given non_redundant can avoid unneccessary redundancy checks
    if non_redundant
        return Polyhedron(Polymake.polytope.Polytope{Polymake.Rational}(VERTICES = points, LINEALITY_SPACE = lineality))
    else
        return Polyhedron(Polymake.polytope.Polytope{Polymake.Rational}(POINTS = points, INPUT_LINEALITY = lineality))
    end
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, P::Polyhedron)
    ad = ambient_dim(P)
    if ad == -1.0
        print(io, "A polyhedron without ambient dimension")
    else
        print(io, "A polyhedron in ambient dimension $(ad)")
    end
end

Polymake.visual(P::Polyhedron; opts...) = Polymake.visual(pm_polytope(P); opts...)
