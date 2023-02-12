@doc Markdown.doc"""
    visualize(P::Polyhedron)

Visualize a polyhedron of dimension at most four (in 3-space).
In dimensions up to 3 a usual embedding is shown.
Four-dimensional polytopes are visualized as a Schlegel diagram, which is a projection onto one of the facets; e.g., see Chapter 5 of [Zie95](@cite).

In higher dimensions there is no standard method; use projections to lower dimensions or try ideas from [GJRW10](@cite).

"""
function visualize(P::Polyhedron)
    pmP = pm_object(P)
    Polymake.visual(pmP)
end


@doc Markdown.doc"""
    visualize(C::Cone)

Visualize a cone.
"""
function visualize(C::Cone)
    pmC = pm_object(C)
    Polymake.visual(pmC)
end


@doc Markdown.doc"""
    visualize(PF::PolyhedralFan)

Visualize a polyhedral fan.
"""
function visualize(PF::PolyhedralFan)
    pmF = pm_object(PF)
    Polymake.visual(pmF)
end


@doc Markdown.doc"""
    visualize(PF::PolyhedralComplex)

Visualize a polyhedral complex.
"""
function visualize(PC::PolyhedralComplex)
    pmPC = pm_object(PC)
    Polymake.visual(pmPC)
end

