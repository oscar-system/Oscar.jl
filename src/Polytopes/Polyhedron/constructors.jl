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


### Construct polyhedron from V-data, as the convex hull of points, rays and lineality.
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

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, P::Polyhedron)
    print(io, "A polyhedron")
end

Polymake.visual(P::Polyhedron; opts...) = Polymake.visual(pm_polytope(P); opts...)

