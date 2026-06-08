#exports
export QuiverRepresentation
export sign_j
export P_gen
export edge_gens
export QuiverGrassmannian
export quiver_grassmannian
#Create quiver representation from directed graphs, ambient vector space dimensions, linear maps

struct QuiverRepresentation
    quiver::Graph{Directed}
    ambient_dims::Vector{Int}
    maps::Vector
    function QuiverRepresentation(quiver, ambient_dims::Vector{Int}, maps::Vector)
        @req n_vertices(quiver) == length(ambient_dims) "each vertex needs an ambient dimension"
        @req n_edges(quiver) == length(maps) "each edge needs a linear map"
        new(quiver, ambient_dims, maps)
    end
end

@doc raw"""
    quiver_representaion(quiver::Graph{Directed}, ambient_dims::Vector{Int}, maps::Vector)

Returns quiver representaion corresponding to a directed graph, ambient dimension vector, and list of linear maps.
"""

function quiver_representation(quiver, ambient_dims,maps)
    return QuiverRepresentation(quiver, ambient_dims,maps)
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
    quiver_representation
    ambient_ring
    defining_ideal
    dimension_vector
end


@doc raw"""
    quiver_grassmannian(quiver::Graph{Directed}, ambient_dims::Vector{Int}, maps::Vector)

Returns coordinate ring of moduli space of representations of the given quiver.
"""
function quiver_grassmannian(Q::QuiverRepresentation,dims::Vector{Int})
    #quiver rep data
    G = Q.quiver
    ns = Q.ambient_dims
    As = Q.maps
    @req length(dims) == length(ns) "each vertex needs a subspace dimension"
    @req all([dims[i]<= ns[i] for i in 1:length(ns)]) "dimension on vertex must be less or
                                                        equal to than ambient dimension"
    #ambient ring
    F = base_ring(As[1])
    Ls = [(i,s) for i in 1:length(ns) for s in subsets(ns[i],dims[i])]
    Is = sort!(Ls)
    R,x = polynomial_ring(F,:x=>Is)
    #index dictionary
    xdict = Dict([Is[i] => x[i] for i in 1:length(Is)])
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
    return QuiverGrassmannian(Q,R,Gs,dims)
end
