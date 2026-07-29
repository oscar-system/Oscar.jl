abstract type GAPDigraph <: AbstractGraph{Directed} end

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
        print(io, "Digraph with ", nv(D), " vertices, ", ne(D), " edges")
    end
end

# Used when constructing a digraph to generate whether the digraph is mutable
function _filt(mut::Bool)
    if mut
        return GAP.Globals.IsMutableDigraph
    else
        return GAP.Globals.IsImmutableDigraph
    end
end

