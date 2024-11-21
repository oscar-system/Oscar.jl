################################################################################
#
# OSCAR late night complaints list
# - edges should know if they are directed or undirected
# - if undirected, Edge(1,3) should be the same as Edge(3,1)
# - groebner basis of the zero ideal should not be empty (?)
#
################################################################################


################################################################################
#
# Graph curves
#
################################################################################
@doc raw"""
    graph_curve(G::Graph)

Return the graph curve of `G`, i.e., a union of lines whose dual graph is `G`, see [BE91](@cite).

Assumes that `G` is trivalent and 3-connected.

# Examples
```jldoctest
julia> G1 = edgegraph(simplex(3))
Undirected graph with 4 nodes and the following edges:
(2, 1)(3, 1)(3, 2)(4, 1)(4, 2)(4, 3)

julia> C1 = graph_curve(G1)
Projective curve
  in projective 2-space over QQ with coordinates [x1, x2, x3]
defined by ideal (x1^2*x2*x3 - x1*x2^2*x3 + x1*x2*x3^2)

julia> G2 = edgegraph(cube(3))
Undirected graph with 8 nodes and the following edges:
(2, 1)(3, 1)(4, 2)(4, 3)(5, 1)(6, 2)(6, 5)(7, 3)(7, 5)(8, 4)(8, 6)(8, 7)

julia> C2 = graph_curve(G2)
Projective curve
  in projective 4-space over QQ with coordinates [x1, x2, x3, x4, x5]
defined by ideal with 5 generators

```
"""
"""
TODO: Move check functions to src/Combinatorics/Graphs/functions.jl file
"""
function check_i_valency(G::Graph, i::int)
  return all(v -> degree(G, v) == i, vertices(G))
end

function check_i_connectivity(G::Graph, i::int)
  for combination in AbstractAlgebra.combinations(collect(edges(G)),i-1)
    G1 =deepcopy(G)
    rem_edge!.(Ref(G1),combination)
    if !is_connected(G1)
      return false 
    end
  end
  return true 
end

function graph_curve(G::Graph; check::Bool=true)
    C,_ = graph_curve_with_vertex_ideals(G; check=check)
    return C
end

function graph_curve_with_vertex_ideals(G::Graph; check::Bool=true)
    if check
        @req check_i_valency(G,3) "G is not trivalent"
        @req check_i_connectivity(G,3) "G is not three-connected"
    end

    R,_ = graded_polynomial_ring(QQ,Int(nv(G)//2+1))
    rowOfVariables = matrix(R,[gens(R)])

    cycleMatrix = kernel(matrix(QQ,signed_incidence_matrix(G)); side=:right)
    cycleMatrix = matrix(R,cycleMatrix) # converting to matrix over R for vcat below

    edgesG = collect(edges(G)) # indexes all edges of G

    vertexIdeals = MPolyIdeal[]
    for v in vertices(G)
        edgesContainingV = findall(edge->(v in edge),edgesG)
        cycleMatrix_v = cycleMatrix[edgesContainingV,:]
        push!(vertexIdeals,ideal(minors(vcat(cycleMatrix_v,rowOfVariables),3)))
    end
    graphCurveIdeal = reduce(intersect,vertexIdeals)
    graphCurve = projective_curve(graphCurveIdeal)

    return graphCurve, vertexIdeals
end

function graph_curve_with_vertex_and_edge_ideals(G::Graph)
    C, vertexIdeals = graph_curve_with_vertex_ideals(G)
    edgeIdeals = Dict{Edge,MPolyIdeal}()
    for e in edges(G)
        edgeIdeals[e] = radical(vertexIdeals[src(e)] + vertexIdeals[dst(e)]) # is radical really necessary?
    end
    return C, vertexIdeals, edgeIdeals
end



################################################################################
#
# Edge pairings
#
################################################################################
function edge_pairing_candidates(G::Graph, e::Edge)
    v1, v2 = src(e), dst(e)
    adj_v1 = filter(v -> v != v2, neighbors(G,v1))
    adj_v2 = filter(v -> v != v1, neighbors(G,v2))
    pairing1 = ( (Edge(max(v1,adj_v1[1]),min(v1,adj_v1[1])),
                  Edge(max(v2,adj_v2[1]),min(v2,adj_v2[1]))),
                 (Edge(max(v1,adj_v1[2]),min(v1,adj_v1[2])),
                  Edge(max(v2,adj_v2[2]),min(v2,adj_v2[2]))) )
    pairing2 = ( (Edge(max(v1,adj_v1[1]),min(v1,adj_v1[1])),
                  Edge(max(v2,adj_v2[2]),min(v2,adj_v2[2]))),
                 (Edge(max(v1,adj_v1[1]),min(v1,adj_v1[1])),
                  Edge(max(v2,adj_v2[2]),min(v2,adj_v2[2]))) )
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


################################################################################
#
# Edge deformations
#
################################################################################
function edge_deformation(e::Edge,edgeIdeals,edgePairing,vertexIdeals::Vector{MPolyIdeal})
    v0 = src(e)
    v1 = dst(e)
    lineIdeal = intersect(vertexIdeals[v0], vertexIdeals[v1])
    println("lineIdeal: ",gens(lineIdeal))
    quadraticGenerator = first(filter(g -> (total_degree(g) == 2), gens(lineIdeal)))
    edgePlane = edge_plane(e,vertexIdeals)
    println("quadraticGenerator: ",quadraticGenerator)
    println("edgePlane: ",groebner_basis(edgePlane))
    if length(groebner_basis(edgePlane))>0
        quadraticGenerator = reduce(quadraticGenerator,groebner_basis(edgePlane)) # why is reduce necessary?
    end
    l1,l2 = edge_pairing_lines(e,edgeIdeals,edgePairing,vertexIdeals)
    l12 = l1*l2
    if length(groebner_basis(edgePlane))>0
        l12 = reduce(l12,groebner_basis(edgePlane)) # why is reduce necessary?
    end
    l12 = deformation_sign(quadraticGenerator,l12)*l12
    return quadraticGenerator,l12
end

function edge_plane(e::Edge,vertexIdeals::Vector{MPolyIdeal})
    v0 = src(e)
    v1 = dst(e)
    planeIdeal = radical(vertexIdeals[v0]*vertexIdeals[v1])
    return ideal_of_linear_generators(planeIdeal)
end

function ideal_of_linear_generators(I::MPolyIdeal)
    return ideal(vcat([zero(base_ring(I))],[g for g in gens(I) if total_degree(g) == 1]))
end

function edge_pairing_lines(e::Edge,edgeIdeals,edgePairing,vertexIdeals::Vector{MPolyIdeal})
    if src(e) < dst(e)
        e = Edge(dst(e),src(e))
    end
    edgePair1, edgePair2 = edgePairing[e]
    edge11 = edgePair1[1]
    edge12 = edgePair1[2]
    line1 = radical(edgeIdeals[edge11]*edgeIdeals[edge12])
    line1 = ideal_of_linear_generators(line1)
    edge21 = edgePair2[1]
    edge22 = edgePair2[2]
    line2 = radical(edgeIdeals[edge21]*edgeIdeals[edge22])
    line2 = ideal_of_linear_generators(line2)
    ePlane = edge_plane(e,vertexIdeals)
    linearForm1 = first(filter(g -> !(g in ePlane), gens(line1)))
    linearForm2 = first(filter(g -> !(g in ePlane), gens(line2)))
    return linearForm1,linearForm2
end

function deformation_sign(q1::MPolyRingElem,q2::MPolyRingElem)
    R = parent(q1)
    K = coefficient_ring(R)
    Rt,tx = polynomial_ring(K,vcat([:t],symbols(R)))
    iota = hom(R,Rt,tx[2:end])
    H1 = iota.(hessian_matrix(q1))
    H2 = iota.(hessian_matrix(q2))
    t = first(tx)
    resultantIdeal = ideal(minors(H1+t*H2,3))
    resultantIdeal = saturation(resultantIdeal,ideal([t]))
    @req ngens(resultantIdeal) == 1 "resultantIdeal should have exactly one generator"
    resultant = first(gens(resultantIdeal))
    return -coeff(resultant,t)
end
