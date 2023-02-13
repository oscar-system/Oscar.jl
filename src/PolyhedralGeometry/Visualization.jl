@doc Markdown.doc"""
    visualize(P::Union{Polyhedron, Cone, PolyhedralFan, PolyhedralComplex})

Visualize a polyhedral object of dimension at most four (in 3-space).
In dimensions up to 3 a usual embedding is shown.
Four-dimensional polytopes are visualized as a Schlegel diagram, which is a projection onto one of the facets; e.g., see Chapter 5 of [Zie95](@cite).

In higher dimensions there is no standard method; use projections to lower dimensions or try ideas from [GJRW10](@cite).
"""

function visualize(P::Union{Polyhedron, Cone, PolyhedralFan, PolyhedralComplex})
    d = ambient_dim(P)
    b = P isa Polyhedron
    if d < 4 || (d == 4 && b && dim(P) == 4)
        pmo = pm_object(P)
        Polymake.visual(pmo)
    else
        d == 4 && b && throw(ArgumentError(string("Can only visualize full-dimensional ", typeof(P), " of ambient dimension ", d, ".")))
        throw(ArgumentError(string("Can not visualize ", typeof(P), " of ambient dimension ", d, ". Supported range: 1 <= d <= ", 3 + b)))
    end
end

