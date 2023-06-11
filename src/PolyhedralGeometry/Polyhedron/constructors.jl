###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct Polyhedron{T<:scalar_types} #a real polymake polyhedron
    pm_polytope::Polymake.BigObject

    # only allowing scalar_types;
    # can be improved by testing if the template type of the `BigObject` corresponds to `T`

    @doc raw"""
        Polyhedron{T}(P::Polymake.BigObject) where T<:scalar_types

    Construct a `Polyhedron` corresponding to a `Polymake.BigObject` of type `Polytope`.
    The type parameter `T` is optional but recommended for type stability.
    """
    Polyhedron{T}(p::Polymake.BigObject) where T<:scalar_types = new{T}(p)
end

# default scalar type: `QQFieldElem`
polyhedron(A, b) = polyhedron(QQFieldElem, A, b)
polyhedron(A) = polyhedron(QQFieldElem, A)

# Automatic detection of corresponding OSCAR scalar type;
# Avoid, if possible, to increase type stability
polyhedron(p::Polymake.BigObject) = Polyhedron{detect_scalar_type(Polyhedron, p)}(p)

@doc raw"""
    polyhedron(::Type{T}, A::AnyVecOrMat, b) where T<:scalar_types

The (convex) polyhedron defined by

$$P(A,b) = \{ x |  Ax ≤ b \}.$$

see Def. 3.35 and Section 4.1. of [JT13](@cite)

# Examples
The following lines define the square $[0,1]^2 \subset \mathbb{R}^2$:
```jldoctest
julia> A = [1 0; 0 1; -1 0 ; 0 -1];

julia> b = [1, 1, 0, 0];

julia> polyhedron(A,b)
Polyhedron in ambient dimension 2
```
"""
polyhedron(::Type{T}, A::AnyVecOrMat, b::AbstractVector) where T<:scalar_types = polyhedron(T, (A, b))

polyhedron(::Type{T}, A::AbstractVector, b::Any) where T<:scalar_types = polyhedron(T, ([A], [b]))

polyhedron(::Type{T}, A::AbstractVector, b::AbstractVector) where T<:scalar_types = polyhedron(T, ([A], b))

polyhedron(::Type{T}, A::AbstractVector{<:AbstractVector}, b::Any) where T<:scalar_types = polyhedron(T, (A, [b]))

polyhedron(::Type{T}, A::AbstractVector{<:AbstractVector}, b::AbstractVector) where T<:scalar_types = polyhedron(T, (A, b))

polyhedron(::Type{T}, A::AnyVecOrMat, b::Any) where T<:scalar_types = polyhedron(T, A, [b])

@doc raw"""
    polyhedron(::Type{T}, I::Union{Nothing, AbstractCollection[AffineHalfspace]}, E::Union{Nothing, AbstractCollection[AffineHyperplane]} = nothing) where T<:scalar_types

The (convex) polyhedron obtained intersecting the halfspaces `I` (inequalities)
and the hyperplanes `E` (equations).

# Examples
The following lines define the square $[0,1]^2 \subset \mathbb{R}^2$:
```jldoctest
julia> A = [1 0; 0 1; -1 0 ; 0 -1];

julia> b = [1, 1, 0, 0];

julia> polyhedron((A,b))
Polyhedron in ambient dimension 2
```

As an example for a polyhedron constructed from both inequalities and
equations, we construct the polytope $[0,1]\times\{0\}\subset\mathbb{R}^2$
```jldoctest
julia> P = polyhedron(([-1 0; 1 0], [0,1]), ([0 1], [0]))
Polyhedron in ambient dimension 2

julia> is_feasible(P)
true

julia> dim(P)
1

julia> vertices(P)
2-element SubObjectIterator{PointVector{QQFieldElem}}:
 [1, 0]
 [0, 0]
```
"""
function polyhedron(::Type{T}, I::Union{Nothing, AbstractCollection[AffineHalfspace]}, E::Union{Nothing, AbstractCollection[AffineHyperplane]} = nothing) where T<:scalar_types
    if isnothing(I) || _isempty_halfspace(I)
        EM = affine_matrix_for_polymake(E)
        IM = Polymake.Matrix{scalar_type_to_polymake[T]}(undef, 0, size(EM, 2))
    else
        IM = -affine_matrix_for_polymake(I)
        EM = isnothing(E) || _isempty_halfspace(E) ? Polymake.Matrix{scalar_type_to_polymake[T]}(undef, 0, size(IM, 2)) : affine_matrix_for_polymake(E)
    end

    return Polyhedron{T}(Polymake.polytope.Polytope{scalar_type_to_polymake[T]}(INEQUALITIES = remove_zero_rows(IM), EQUATIONS = remove_zero_rows(EM)))
end

"""
    pm_object(P::Polyhedron)

Get the underlying polymake `Polytope`.
"""
pm_object(P::Polyhedron) = P.pm_polytope

function ==(P0::Polyhedron, P1::Polyhedron)
    # TODO: Remove the following 3 lines, see #758
    for pair in Iterators.product([P0, P1], ["RAYS", "FACETS"])
        Polymake.give(pm_object(pair[1]),pair[2])
    end
    Polymake.polytope.equal_polyhedra(pm_object(P0), pm_object(P1))
end


### Construct polyhedron from V-data, as the convex hull of points, rays and lineality.
@doc raw"""
    convex_hull([::Type{T} = QQFieldElem,] V [, R [, L]]; non_redundant::Bool = false)

Construct the convex hull of the vertices `V`, rays `R`, and lineality `L`. If
`R` or `L` are omitted, then they are assumed to be zero.

# Arguments
- `V::AbstractCollection[PointVector]`: Points whose convex hull is to be computed.
- `R::AbstractCollection[RayVector]`: Rays completing the set of points.
- `L::AbstractCollection[RayVector]`: Generators of the Lineality space.

If an argument is given as a matrix, its content has to be encoded row-wise.

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
Polyhedron in ambient dimension 2
```
To construct the positive orthant, rays have to be passed:
```jldoctest
julia> V = [0 0];

julia> R = [1 0; 0 1];

julia> PO = convex_hull(V, R)
Polyhedron in ambient dimension 2
```
The closed-upper half plane can be constructed by passing rays and a lineality space:
```jldoctest
julia> V = [0 0];

julia> R = [0 1];

julia> L = [1 0];

julia> UH = convex_hull(V, R, L)
Polyhedron in ambient dimension 2
```
To obtain the x-axis in $\mathbb{R}^2$:
```jldoctest
julia> V = [0 0];

julia> R = nothing;

julia> L = [1 0];

julia> XA = convex_hull(V, R, L)
Polyhedron in ambient dimension 2
```
"""
function convex_hull(::Type{T}, V::AbstractCollection[PointVector], R::Union{AbstractCollection[RayVector], Nothing} = nothing, L::Union{AbstractCollection[RayVector], Nothing} = nothing; non_redundant::Bool = false) where T<:scalar_types
    # Rays and Points are homogenized and combined and
    # Lineality is homogenized
    points = stack(homogenized_matrix(V, 1), homogenized_matrix(R, 0))
    lineality = isnothing(L) || isempty(L) ? zero_matrix(QQ, 0, size(points,2)) : homogenized_matrix(L, 0)

    # These matrices are in the right format for polymake.
    # given non_redundant can avoid unnecessary redundancy checks
    if non_redundant
        return Polyhedron{T}(Polymake.polytope.Polytope{scalar_type_to_polymake[T]}(VERTICES = points, LINEALITY_SPACE = lineality))
    else
        return Polyhedron{T}(Polymake.polytope.Polytope{scalar_type_to_polymake[T]}(POINTS = remove_zero_rows(points), INPUT_LINEALITY = remove_zero_rows(lineality)))
    end
end

convex_hull(V::AbstractCollection[PointVector], R::Union{AbstractCollection[RayVector], Nothing} = nothing, L::Union{AbstractCollection[RayVector], Nothing} = nothing; non_redundant::Bool = false) = convex_hull(QQFieldElem, V, R, L; non_redundant=non_redundant)

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, P::Polyhedron{T}) where T<:scalar_types
    try
        ad = ambient_dim(P)
        print(io, "Polyhedron in ambient dimension $(ad)")
        T != QQFieldElem && print(io, " with $T type coefficients")
    catch e
        print(io, "Polyhedron without ambient dimension")
    end
end

Polymake.visual(P::Polyhedron; opts...) = Polymake.visual(pm_object(P); opts...)
