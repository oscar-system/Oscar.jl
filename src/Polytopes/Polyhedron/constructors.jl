###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

@doc Markdown.doc"""

    Polyhedron(P::Polymake.BigObject)

Construct a `Polyhedron` corresponding to a polymake `Polytope`.
"""
function Polyhedron(pm_polytope::Polymake.BigObject)
    Polyhedron(pm_polytope, :unknown)
end

@doc Markdown.doc"""

    Polyhedron(A::Union{Oscar.MatElem,AbstractMatrix}, b)

The (convex) polyhedron defined by

$$P(A,b) = \{ x |  Ax ≤ b \}.$$

see Def. 3.35 and Section 4.1. of [JT13](@cite)

# Arguments
- `A::Matrix`: Matrix corresponding to the linear coefficients of the inequalilites that describe P.
- `b::Vector`: Vector corresponding to the constant term of the inequalilites that describe P.

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

Polyhedron(H::HalfSpaceIterator) = Polyhedron(H.A, H.b)

"""
    pm_polytope(P::Polyhedron)

Get the underlying polymake `Polytope`.
"""
pm_polytope(P::Polyhedron) = P.pm_polytope

==(P0::Polyhedron, P1::Polyhedron) = Polymake.polytope.equal_polyhedra(pm_polytope(P0), pm_polytope(P1))


### Construct polyhedron from V-data, as the convex hull of points, rays and lineality.
@doc Markdown.doc"""
    convex_hull(V [, R [, L]])

The polytope given as the convex hull of the rows of a set of points.

see Def. 2.11 and Def. 3.1  of [JT13](@cite)

# Arguments
- `V::Matrix`: Points whose convex hull is to be computed; encoded as row vectors.
- `R::Matrix`: Rays completing the set of points; encoded row-wise as representative vectors.
- `L::Matrix`: Generators of the Lineality space; encoded as row vectors.

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
function convex_hull(V::Union{PointIterator{Points}, AnyVecOrMat, Oscar.MatElem}, R::Union{PointIterator{Ray}, AnyVecOrMat, Oscar.MatElem, Nothing} = nothing, L::Union{PointIterator, AnyVecOrMat, Oscar.MatElem, Nothing} = nothing; non_redundant::Bool = false)
    # we access the matrices which polymake can work with.
    VM = V isa PointIterator{Points} ? V.m : matrix_for_polymake(V)
    RM = R isa PointIterator{Ray} ? R.m : R isa Union{AnyVecOrMat, Oscar.MatElem} ? matrix_for_polymake(R) : Polymake.Matrix{Polymake.Rational}(undef, 0, size(V.m, 2))
    LM = L isa PointIterator ? L.m : L isa Union{AnyVecOrMat, Oscar.MatElem} ? matrix_for_polymake(L) : Polymake.Matrix{Polymake.Rational}(undef, 0, size(V.m, 2))

    # Rays and Points are homogenized and combined and
    # Lineality is homogenized
    points = stack(homogenize(VM, 1), homogenize(RM, 0))
    lineality = homogenize(LM, 0)

    # These matrices are in the right format for polymake.
    # given non_redundant can avoid unneccessary redundancy checks
    if non_redundant
        return Polyhedron(Polymake.polytope.Polytope{Polymake.Rational}(VERTICESS = points, LINEALITY_SPACE = lineality))
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
    print(io, "A polyhedron")
end

Polymake.visual(P::Polyhedron; opts...) = Polymake.visual(pm_polytope(P); opts...)
