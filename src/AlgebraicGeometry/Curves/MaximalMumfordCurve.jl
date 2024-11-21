function edge_pairing_candidates(G::Graph, e::Edge)
    v1, v2 = src(e), dst(e)
    adj_v1 = filter(v -> v != v2, neighbors(G,v1))
    adj_v2 = filter(v -> v != v1, neighbors(G,v2))
    pairing1 = ( (Edge(v1,adj_v1[1]),Edge(v2,adj_v2[1])), (Edge(v1,adj_v1[2]),Edge(v2,adj_v2[2])) )
    pairing2 = ( (Edge(v1,adj_v1[1]),Edge(v2,adj_v2[2])), (Edge(v1,adj_v1[1]),Edge(v2,adj_v2[2])) )
    return [pairing1,pairing2]
end

function maximal_edge_pairing(G::Graph)
    pGfaces = faces(IncidenceMatrix, tutte_lifting(G), 2)
    EP = Dict{Edge, Tuple{Tuple{Edge, Edge},Tuple{Edge, Edge}}}()
    for edge in edges(G)
        pairing1, pairing2 =  edge_pairing_candidates(G,edge)
        # take the first pair of edges of each pairing candidate
        # the second pair of edges is not required
        # as the first pair is correct if and only if the second pair is
        edge11 = pairing1[1][1]
        edge12 = pairing1[1][2]
        edge21 = pairing2[1][1]
        edge22 = pairing2[1][2]
        # for each pair of edges, construct the set of four vertices
        pairingVertices1 = [src(edge11),dst(edge11),src(edge12),dst(edge12)]
        pairingVertices2 = [src(edge21),dst(edge21),src(edge22),dst(edge22)]
        # and check whether there is a tutte lifting facet containing all of them
        for i in 1:nrows(pGfaces)
            if all(pGfaces[i,pairingVertices1])
                EP[edge] = pairing1
                break
            elseif all(pGfaces[i,pairingVertices2])
                EP[edge] = pairing2
                break
            end
        end
    end
    return EP
end
