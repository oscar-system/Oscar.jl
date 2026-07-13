mutable struct Digraph{T<:GraphTypes}
    X::GapObj
    nv::Int64
    ne::Int64
    mut::Bool

    function Digraph{T}(d::GapObj) where T<:GraphTypes
        @req DigraphWrap.IsDigraph(d) "the GAP object must be a digraph"
        nv = DigraphWrap.DigraphNrVertices(d)
        ne = DigraphWrap.DigraphNrEdges(d)
        mut = DigraphWrap.IsMutableDigraph(d)
        z = new{T}(d, nv, ne, mut)
        return z
    end
end

GAP.@install GapObj(D::Digraph) = D.X

const DirectedDigraph = Digraph{Directed}
const UndirectedDigraph = Digraph{Undirected}

digraph(d::GapObj) = Digraph{Directed}(d)
undirected_digraph(d::GapObj) = Digraph{Undirected}(d)

function Base.show(io::IO, ::MIME"text/plain", D::Digraph)
    @show_name(io, D)
    @show_special(io, D)
    if is_terse(io)
        print(io, "Digraph")
    else
        print(io, "Digraph with ", nv(D), " vertices, ", ne(D), " edges")
    end
end
