@doc Markdown.doc"""
    visualize(P::Polyhedron)

Visualize a polyhedron.
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
    visualize(G::Graphs.Graph{T}) where {T <: Union{Polymake.Directed, Polymake.Undirected}}

Visualize a graph.
"""
function visualize(G::Graphs.Graph{T}) where {T <: Union{Polymake.Directed, Polymake.Undirected}}
    BigGraph = Polymake.graph.Graph(ADJACENCY=pm_object(G))
    Polymake.visual(BigGraph)
end
