abstract type GAPDigraph <: AbstractGraph{Directed} end

@doc raw"""
    Digraph

A digraph with vertices `1, ..., n`, backed by the GAP package Digraphs.
Every `Digraph` stores the underlying GAP object together with the numbers
of vertices and edges and a flag indicating mutability; see
[`digraph`](@ref) for the constructors, and [`nv`](@ref) and [`ne`](@ref)
for the basic attributes.

```jldoctest
julia> D = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges
```
"""
@attributes mutable struct Digraph <: GAPDigraph
    X::GapObj
    nv::Int64
    ne::Int64
    mut::Bool

    function Digraph(d::GapObj)
        @req DigraphWrap.IsDigraph(d) "the GAP object must be a digraph"
        nv = DigraphWrap.DigraphNrVertices(d)
        ne = DigraphWrap.DigraphNrEdges(d)
        mut = DigraphWrap.IsMutableDigraph(d)
        z = new(d, nv, ne, mut)
        return z
    end
end

GAP.@install GapObj(D::Digraph) = D.X

function Base.show(io::IO, ::MIME"text/plain", D::Digraph)
    @show_name(io, D)
    @show_special(io, D)
    if is_terse(io)
        print(io, "Digraph")
    else
        print(io, "Digraph with ", nv(D), nv(D) == 1 ? " vertex, " : " vertices, ",
              ne(D), ne(D) == 1 ? " edge" : " edges")
    end
end
