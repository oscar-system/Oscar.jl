#exports
export QuiverRepresentation
export QuiverGrassmannian
export quiver_grassmannian
export quiver_representation

#Create quiver representation from directed graphs, ambient vector space dimensions, linear maps
struct QuiverRepresentation
    quiver::Graph{Directed}
    ambient_dims::Vector{Int}
    edge_morphisms::Vector
    function QuiverRepresentation(quiver, ambient_dims::Vector{Int}, input_matrices::Vector)
        @req n_vertices(quiver) == length(ambient_dims) "each vertex needs an ambient dimension"
        @req n_edges(quiver) == length(input_matrices) "each edge needs a linear map"
        new(quiver, ambient_dims, input_matrices, edge_morphisms(quiver,input_matrices,ambient_dims))
    end
end

#convert matrices to free module morphisms
function edge_morphisms(G, As, ns)
    R = base_ring(As[1])
    return [
        hom(free_module(R, ns[src(e)]),
            free_module(R, ns[dst(e)]),
            transpose(A))
        for (A, e) in zip(As, edges(G))
    ]
end

@doc raw"""
    quiver_representation(quiver::Graph{Directed}, ambient_dims::Vector{Int}, maps::Vector)

Returns `QuiverRepresentation` object corresponding to a directed graph `quiver`, ambient dimension vector `ambient_dims` corresponding to vertices, and list of linear maps `maps` corresponding to the edges of the graph.

# Examples
```jldoctest
julia> G = graph_from_edges(Directed, [[1,2]])
Directed graph with 2 nodes and the following edges:
(1, 2)

julia> A = transpose(matrix(QQ,[1 0 0 0;0 1 0 0]))
[1   0]
[0   1]
[0   0]
[0   0]

julia> Q = quiver_representation(G,[2,4],[A])
QuiverRepresentation(Directed graph with 2 nodes and 1 edges, [2, 4], QQMatrix[[1 0; 0 1; 0 0; 0 0]])
```
"""
function quiver_representation(quiver, ambient_dims, maps)
    return QuiverRepresentation(quiver, ambient_dims, maps)
end

####Internal functions for Quiver Grassmannian
function sign_j(j,I,J)
    g = count(>(j), J) + count(>(j), I)  
    return (-1)^(g)
end

#returns generator for (I,J) pair associated with an edge
function P_gen(A,I,J,e,n1,xdict)
    N = 1:n1
    return sum(sign_j(j, I, J)*A[i,j]*xdict[(src(e),sort(union(I,j)))]*xdict[(dst(e),setdiff(J,i))] for 
                j in setdiff(N,I), i in J);# init=0)
end

#generators associated with an edge
function edge_gens(e, nsi, dsi, A, xdict)
    #variable labels
    L1 = subsets(nsi[1], dsi[1]-1)
    L2 = subsets(nsi[2], dsi[2]+1)
    #pairs
    X = [(a, b) for a in L1, b in L2 if length(b) - length(a) >= 2]
    #generators
    T = [P_gen(A, I, J, e, nsi[1], xdict) for (I,J) in X]
    return unique(filter!(!iszero, T))
end


#Creat coordinate ring of Quiver Grassmannian
struct QuiverGrassmannian#add types
    quiver_representation::QuiverRepresentation
    ambient_ring::MPolyRing
    defining_ideal::MPolyIdeal
    dimension_vector::Vector{Int64}
end

#generate free module homomorphisms from matrices
function edge_morphisms(A,G,)
    R = base_ring(A[1])
    return [hom(A,free_module(R,ns))]
end


@doc raw"""
     quiver_grassmannian(Q::QuiverRepresentation, dims::Vector{Int})

Return the data of the quiver Grassmannian parameterizing subrepresentations of a given quiver representation `Q` with subspace dimension vector `dims` on each node. Explicitly, the function returns a `QuiverGrassmannian` object, containing the underlying quiver representation, ambient ring, defining ideal equations, and the dimension vector.
# Examples
```jldoctest
julia> G = graph_from_edges(Directed, [[1,2]])
Directed graph with 2 nodes and the following edges:
(1, 2)

julia> A = transpose(matrix(QQ,[1 0 0 0;0 1 0 0]))
[1   0]
[0   1]
[0   0]
[0   0]

julia> Q = quiver_representation(G,[2,4],[A])
QuiverRepresentation(Directed graph with 2 nodes and 1 edges, [2, 4], QQMatrix[[1 0; 0 1; 0 0; 0 0]])

julia> Qsr = quiver_grassmannian(Q,[1,2])
QuiverGrassmannian(QuiverRepresentation(Directed graph with 2 nodes and 1 edges, [2, 4], QQMatrix[[1 0; 0 1; 0 0; 0 0]]), Multivariate polynomial ring in 8 variables over QQ, Ideal with 5 generators, [1, 2])
```
"""
function quiver_grassmannian(Q::QuiverRepresentation,dims::Vector{Int})
    #quiver rep data
    G = Q.quiver
    ns = Q.ambient_dims
    As = Q.input_matrices
    @req length(dims) == length(ns) "each vertex needs a subspace dimension"
    @req all([dims[i]<= ns[i] for i in 1:length(ns)]) "dimension on vertex must be less or
                                                        equal to than ambient dimension"
    #ambient ring
    F = base_ring(As[1])
    Ls = [(i,s) for i in 1:length(ns) for s in subsets(ns[i],dims[i])]
    Is = sort!(Ls)
    R,x = polynomial_ring(F,:x=>Is)
    #index dictionary
    xdict = Dict(Is[i] => x[i] for i in 1:length(Is))
    #create ideal generators for each edge
    Gs = elem_type(R)[]
    for (e, A) in zip(edges(G), As)
        #quiver generators for node
        nsi = [ns[src(e)],ns[dst(e)]]
        dsi = [dims[src(e)],dims[dst(e)]]
        Ge = edge_gens(e,nsi,dsi,A,xdict)
        append!(Gs,Ge)
    end
    #create grassmann generators for each node
    for j in 1:length(ns)
        if ns[j]-1 >dims[j]>1
            Jj = filter(t -> t[1] == j, Is)
            xx = [xdict[t] for t in Jj]
            Gr_di_ni = gens(grassmann_pluecker_ideal(dims[j],ns[j]))
            phi_i = hom(parent(Gr_di_ni[1]),R,xx)  
            phiG = phi_i.(Gr_di_ni)
            append!(Gs,phiG)
        end
    end
    return QuiverGrassmannian(Q,R,ideal(Gs),dims)
end
