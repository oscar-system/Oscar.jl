#exports
export QuiverRepresentation
export QuiverGrassmannian
export quiver_grassmannian
export quiver_representation
#checks that the dimensions of the input matrices match the dimension labels on vertices
function check_matrix_dimensions(quiver::Graph{Directed}, ambient_dims::Vector{Int}, input_matrices::AbstractVector{<:MatElem})
    for (e, A) in zip(edges(quiver), input_matrices)
        u = src(e)
        v = dst(e)

        @req nrows(A) == ambient_dims[u] begin
            "Matrix for edge $u → $v has $(nrows(A)) rows, " *
            "but the source vertex has dimension $(ambient_dims[u]). " *
            "Check input matrices satisfy OSCAR's freemodule morphism convention."
        end

        @req ncols(A) == ambient_dims[v] begin
            "Matrix for edge $u → $v has $(ncols(A)) columns, " *
            "but the target vertex has dimension $(ambient_dims[v]). " *
            "Check input matrices satisfy OSCAR's freemodule morphism convention."
        end
    end
end
#Create quiver representation from directed graphs, ambient vector space dimensions, linear maps
struct QuiverRepresentation{C <: FieldElem}
    quiver::Graph{Directed}
    ambient_dims::Vector{Int}
    vertex_vector_spaces::Vector{Generic.FreeModule{C}}
    edge_morphisms::Vector{Generic.ModuleHomomorphism{C}}
    base_field::Field # elem_type(C)
    function QuiverRepresentation(quiver::Graph{Directed}, ambient_dims::Vector{Int}, input_matrices::AbstractVector{<:MatElem}, base_field::Field)
        @req n_vertices(quiver) == length(ambient_dims) "each vertex needs an ambient dimension"
        @req n_edges(quiver) == length(input_matrices) "each edge needs a linear map"
        try
            input_matrices = change_base_ring.(Ref(base_field), input_matrices)
        catch
            error("The matrix entries cannot be coerced into the specified ring")
        end
        check_matrix_dimensions(quiver, ambient_dims, input_matrices)
        vertex_vector_spaces = [vector_space(base_field,ambient_dims[i]) for i in vertices(quiver)]
        new{elem_type(base_field)}(quiver, ambient_dims, vertex_vector_spaces, edge_morphisms(quiver, input_matrices, vertex_vector_spaces), base_field)
    end
end

#convert matrices to free module morphisms
function edge_morphisms(G::Graph{Directed}, As::AbstractVector{<:MatElem}, vertex_vector_spaces::Vector{Generic.FreeModule{C}}) where {C <: FieldElem}
    Vs = vertex_vector_spaces
    return [
        hom(Vs[src(e)],
            Vs[dst(e)],
            A)
        for (A, e) in zip(As, edges(G))
    ]
end
function Base.show(io::IO, qR::QuiverRepresentation)
    print(io, "Quiver representation over ",qR.base_field," with ambient dimensions ",qR.ambient_dims, " and ",length(qR.edge_morphisms)," arrows")
end

@doc raw"""
    quiver_representation(quiver::Graph{Directed}, ambient_dims::Vector{Int}, maps::AbstractVector{<:MatElem}, base_field::Field)

Return a `QuiverRepresentation` object corresponding to a directed graph `quiver`, ambient dimension vector `ambient_dims` corresponding to vertices, and list of linear maps `maps` corresponding to the edges of the graph.

# Examples
```jldoctest
julia> G = graph_from_edges(Directed, [[1,2]])
Directed graph with 2 nodes and the following edges:
(1, 2)

julia> A = matrix(QQ, [1 0 0 0; 0 1 0 0])
[1   0   0   0]
[0   1   0   0]

julia> Q = quiver_representation(G, [2,4], [A], QQ)
Quiver representation over Rational field with ambient dimensions [2, 4] and 1 arrows
```
"""
function quiver_representation(quiver::Graph{Directed}, ambient_dims::Vector{Int}, maps::AbstractVector{<:MatElem}, base_field::Field)
    return QuiverRepresentation(quiver, ambient_dims, maps, base_field)
end

####Internal functions for Quiver Grassmannian
function sign_j(j::Int, I::Vector{Int}, J::Vector{Int})
    g = count(>(j), J) + count(>(j), I)  
    return (-1)^(g)
end

#returns generator for (I,J) pair associated with an edge
function P_gen(A::MatElem, I::Vector{Int}, J::Vector{Int}, e::Edge, n1::Int, xdict::AbstractDict)
    N = 1:n1
    return sum(sign_j(j, I, J)*A[i,j]*xdict[(src(e),sort(union(I,j)))]*xdict[(dst(e),setdiff(J,i))] for 
                j in setdiff(N,I), i in J);# init=0)
end

#generators associated with an edge
function edge_gens(e::Edge, nsi::Vector{Int}, dsi::Vector{Int}, A::MatElem, xdict::AbstractDict)
    #variable labels
    L1 = subsets(nsi[1], dsi[1]-1)
    L2 = subsets(nsi[2], dsi[2]+1)
    #pairs
    X = [(a, b) for a in L1, b in L2 if length(b) - length(a) >= 2]
    #generators
    T = [P_gen(A, I, J, e, nsi[1], xdict) for (I,J) in X]
    return unique(filter!(!iszero, T))
end
#creates weights for graded ring
function grading_weights(Ls::Vector{Tuple{Int, Vector{Int}}}, e::Vector{Int})
    return [begin
        z = zeros(Int, length(e))
        z[l[1]] = 1
        z
    end for l in Ls]
end

struct QuiverGrassmannian{C <: FieldElem}
    quiver_representation::QuiverRepresentation{C}
    ambient_ring::MPolyRing
    defining_ideal::MPolyIdeal
    dimension_vector::Vector{Int}
    base_field::Field
end

function Base.show(io::IO, Q::QuiverGrassmannian)
    print(io, "Quiver Grassmannian over ",  Q.base_field," with subspace dimensions ", Q.dimension_vector, " defined by ", Q.defining_ideal)
end

@doc raw"""
     quiver_grassmannian(Q::QuiverRepresentation, dims::Vector{Int})

Return the data of the quiver Grassmannian parameterizing subrepresentations of a given quiver representation `Q` with subspace dimension vector `dims` on each node. Explicitly, the function returns a `QuiverGrassmannian` object, containing the underlying quiver representation, ambient ring, defining ideal equations, and the dimension vector.
# Examples
```jldoctest
julia> G = graph_from_edges(Directed, [[1,2]])
Directed graph with 2 nodes and the following edges:
(1, 2)

julia> A = matrix(QQ, [1 0 0 0; 0 1 0 0])
[1   0   0   0]
[0   1   0   0]

julia> Q = quiver_representation(G, [2,4], [A], QQ)
Quiver representation over Rational field with ambient dimensions [2, 4] and 1 arrows

julia> Qsr = quiver_grassmannian(Q, [1,2])
Quiver Grassmannian over Rational field with subspace dimensions [1, 2] defined by Ideal with 5 generators
```
"""
function quiver_grassmannian(Q::QuiverRepresentation, dims::Vector{Int})
    #quiver rep data
    G = Q.quiver
    ns = Q.ambient_dims
    @req length(dims) == length(ns) "each vertex needs a subspace dimension"
    @req all([dims[i] <= ns[i] for i in 1:length(ns)]) "subspace dimension of vertex must be less than or
                                                        equal to the ambient dimension"
    F = Q.base_field
    #create labels for ambient ring variables
    Ls = [(i, s) for i in 1:length(ns) for s in subsets(ns[i], dims[i])]
    sort!(Ls)
    #create ring
    R, x = graded_polynomial_ring(F, :x=>Ls; weights = grading_weights(Ls, ns))
    #index dictionary
    xdict = Dict(Ls[i] => x[i] for i in 1:length(Ls))
    #create ideal generators for each edge
    Gs = elem_type(R)[]
    for (e, A) in zip(edges(G), Q.edge_morphisms)
        #quiver generators for node
        nsi = [ns[src(e)], ns[dst(e)]]
        dsi = [dims[src(e)], dims[dst(e)]]
        Ge = edge_gens(e, nsi, dsi, transpose(matrix(A)), xdict)
        append!(Gs, Ge)
    end
    #create grassmann generators for each node
    for j in 1:length(ns)
        if ns[j]-1 > dims[j] > 1
            Jj = filter(t -> t[1] == j, Ls)
            xx = [xdict[t] for t in Jj]
            Gr_di_ni = gens(grassmann_pluecker_ideal(dims[j], ns[j]))
            phi_i = hom(parent(Gr_di_ni[1]), R, xx)  
            phiG = phi_i.(Gr_di_ni)
            append!(Gs, phiG)
        end
    end
    return QuiverGrassmannian(Q, R, ideal(Gs), dims, F)
end
