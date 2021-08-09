###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

@doc Markdown.doc"""
    Cone(R [, L])

# Arguments
- `R::Matrix`: Rays generating the cone; encoded row-wise as representative vectors.
- `L::Matrix`: Generators of the Lineality space; encoded as row vectors.

A polyhedral cone, not necessarily pointed, defined by the positive hull
of the rays, Rays.

# Examples
To construct the positive orthant as a `Cone`, you can write:
```julia-repl
julia> R = [1 0; 0 1];

julia> PO = Cone(R)
A polyhedral cone in ambient dimension 2
```

To obtain the upper half-space of the plane:
```julia-repl
julia> R = [0 1];

julia> L = [1 0];

julia> HS = Cone(R, L)
A polyhedral cone in ambient dimension 2
```
"""
function Cone(R::Union{PointIterator{Ray}, Oscar.MatElem, AbstractMatrix}, L::Union{PointIterator{Ray}, Oscar.MatElem, AbstractMatrix, Nothing} = nothing; non_redundant::Bool = false)
    RM = R isa PointIterator{Ray} ? R.m : R isa Union{AnyVecOrMat, Oscar.MatElem} ? matrix_for_polymake(R) : Polymake.Matrix{Polymake.Rational}(undef, 0, size(V.m, 2))
    LM = L isa PointIterator ? L.m : L isa Union{AnyVecOrMat, Oscar.MatElem} ? matrix_for_polymake(L) : Polymake.Matrix{Polymake.Rational}(undef, 0, size(V.m, 2))

    if non_redundant
        return Cone(Polymake.polytope.Cone{Polymake.Rational}(RAYS = RM, LINEALITY_SPACE = LM,))
    else
        return Cone(Polymake.polytope.Cone{Polymake.Rational}(INPUT_RAYS = RM, INPUT_LINEALITY = LM,))
    end
end

==(C0::Cone, C1::Cone) = Polymake.polytope.equal_polyhedra(pm_cone(C0), pm_cone(C1))


"""
    positive_hull(generators)

A polyhedral cone, not necessarily pointed, defined by the positive hull
of the `generators`. Redundant rays are allowed in the generators.

# Arguments
- `generators::Matrix`: Rays generating the cone; encoded row-wise as representative vectors.

# Examples
```julia-repl
julia> R = [1 0; 0 1];

julia> PO = positive_hull(R)
A polyhedral cone in ambient dimension 2
```
"""
function positive_hull(generators::Union{Oscar.MatElem,AbstractMatrix})
    # TODO: Filter out zero rows
    C=Polymake.polytope.Cone{Polymake.Rational}(INPUT_RAYS =
      matrix_for_polymake(remove_zero_rows(generators)))
    Cone(C)
end


"""
    pm_cone(C::Cone)

Get the underlying polymake `Cone`.
"""
pm_cone(C::Cone) = C.pm_cone


###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################

function Base.show(io::IO, C::Cone)
    print(io,"A polyhedral cone in ambient dimension $(ambient_dim(C))")
end

Polymake.visual(C::Cone; opts...) = Polymake.visual(pm_cone(C); opts...)
