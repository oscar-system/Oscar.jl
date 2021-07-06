###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct Polyhedron #a real polymake polyhedron
    pm_polytope::Polymake.BigObject
    boundedness::Symbol # Values: :unknown, :bounded, :unbounded
end

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
```julia-repl
julia> A = [1 0; 0 1; -1 0 ; 0 -1];
julia> b = [1, 1, 0, 0];
julia> Polyhedron(A,b)
A polyhedron in ambient dimension 2
```
"""
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


### Construct polyhedron from V-data, as the convex hull of points, rays and lineality.
@doc Markdown.doc"""
    convex_hull(V::Matrix [, R::Matrix [, L::Matrix]]; non_redundant::Bool = false)

The polytope given as the convex hull of the rows of a set of points.

The matrices rows are the points, the rays and the generators of the lineality space,
respectively.

`R` can be given as an empty matrix or as `nothing` if the user wants to compute
the convex hull only from `V` and `L`.

If the user is sure that `V` and `R` only contains extreme points and that
the description of the lineality space is complete, they can set
`non_redundant = true` to avoid unneccessary redundancy checks.

See Def. 2.11 and Def. 3.1  of [JT13](@cite).

# Examples
The following lines define the square $[0,1]^2 \subset \mathbb{R}^2$:
```julia-repl
julia> Square = convex_hull([0 0; 0 1; 1 0; 1 1])
A polyhedron in ambient dimension 2
```
To construct the positive orthant, rays have to be passed:
```julia-repl
julia> V = [0 0];
julia> R = [1 0; 0 1];
julia> PO = convex_hull(V, R)
A polyhedron in ambient dimension 2
```
The closed-upper half plane can be constructed by passing rays and a lineality space:
```julia-repl
julia> V = [0 0];
julia> R = [0 1];
julia> L = [1 0];
julia> UH = convex_hull(V, R, L)
A polyhedron in ambient dimension 2
```
To obtain the x-axis in $\mathbb{R}^2$:
```julia-repl
julia> V = [0 0];
julia> R = [];
julia> L = [1 0];
julia> XA = convex_hull(V, R, L)
A polyhedron in ambient dimension 2
```
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

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, P::Polyhedron)
    print(io, "A polyhedron")
end

Polymake.visual(P::Polyhedron; opts...) = Polymake.visual(pm_polytope(P); opts...)
