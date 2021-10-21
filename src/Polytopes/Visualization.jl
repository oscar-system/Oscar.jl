@doc Markdown.doc"""
    visualize(P::Polyhedron)

Visualize a polyhedron.
"""
function visualize(P::Polyhedron)
    pmP = pm_polytope(P)
    Polymake.visual(pmP)
end


@doc Markdown.doc"""
    visualize(C::Cone)

Visualize a cone.
"""
function visualize(C::Cone)
    pmC = pm_cone(C)
    Polymake.visual(pmC)
end


@doc Markdown.doc"""
    visualize(PF::PolyhedralFan)

Visualize a polyhedral fan.
"""
function visualize(PF::PolyhedralFan)
    pmF = pm_fan(PF)
    Polymake.visual(pmF)
end
