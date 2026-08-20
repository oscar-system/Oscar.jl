# Functions for generating simple higher-level phylogenetic networks 

using Oscar


#function which generates level-2 sunlets obtained by gluing level-1 sunlets with j and k leaves along an edge 
#TODO: mod out by action of S_n on the leaves 
# or more generally, filter by graphs which give the same model (permutation of coordinates)

function glued_sunlets(j,k) 
    n = j+k
    leaves = 1:n
    interior_nodes = (n+1):2*n
    #construct graph of only leaves 
    G = graph_from_edges(Directed, [[interior_nodes[i], leaves[i]] for i in 1:n])
    hybs1 = (n+1):(n+j-1)
    hybs2 = (n+j+2):2*n
    partial_sunlets = [] 
    for h1 in hybs1
        rts = [i for i in (n+1):(n+j+1) if i != h1]
        for rt in rts 
            H = copy(G)
            orient_cycle!(H, (n+1),(n+j+1), rt, h1)
            push!(partial_sunlets, H)
        end 
           
    end 
    sunlets = [] 
    for h2 in hybs2 
        rts = [i for i in (n+j):2*n if i!= h2]
        for H in partial_sunlets
            for rt in rts 
                F = copy(H)
                orient_cycle!(F,(n+j), 2*n, rt, h2)
                #ensure gluing doesn't create any loops and/or hybrid vertices 
                if is_acyclic(F) && length(filter(v -> indegree(F,v) > 1,  vertices(F))) <3 
                    push!(sunlets, F)
                end 
            end 

        end 
    end 


    return unique(sunlets) 
end 

function orient_cycle!(G, a, b, root, hybrid)
    cycle = collect(a:b)
    m = length(cycle)

    ir = findfirst(==(root), cycle)
    ih = findfirst(==(hybrid), cycle)
    #ix = middle == 0 ? 0 : findfirst(==(middle),cycle)

    # Walk from root clockwise until hybrid
    i = ir
    while i != ih
        j = i == m ? 1 : i + 1
        add_edge!(G, cycle[i], cycle[j])
        i = j
    end

    # Walk from root counter-clockwise until hybrid
    i = ir
    while i != ih
        j = i == 1 ? m : i - 1
        add_edge!(G, cycle[i], cycle[j])
        i = j
    end
end


function max_outdeg(G)
        V = vertices(G)
        L = [outdegree(G,v) for v in V ]
    return maximum(L)
end

function max_indeg(G)
        V = vertices(G)
        L = [indegree(G,v) for v in V ]
    return maximum(L)
end 

#some individual graphs

#= 

### Level 1 

# 4 sunlet

G_1_4_1 = graph_from_edges(Directed,[[5,1],[6,2], [7,3],[8,4],[5,6],[6,7],[7,8],[5,8]])


### Level 2 

# 2 leaves  

#diamond 

G_2_2 = graph_from_edges(Directed, [[3,1],[4,2],[5,4],[6,4],[5,3],[6,3],[7,5],[7,6]])

# 3 leaves 

#diamond with "interior leaf" 

#FIXME: doesn't have the right leaf set, apply a permutation to fix this! 
G_2_3_1 =  graph_from_edges(Directed, [[3,1],[4,2],[5,4],[6,4],[5,3],[6,3],[7,5],[7,6],[7,8]])

A = [[3,1],[4,2],[5,4],[6,4],[5,3],[6,3],[7,5],[7,6],[7,8]] 

perm = [1,2,3,7,5,6,4]
[permute!(a,perm) for a in A]
# pentagon 

G_2_3_2 = graph_from_edges(Directed, [[8,1],[4,8],[5,8],[5,4],[4,6],[5,7],[7,6],[6,2],[7,3]])

# 4 leaves 

#pentagon with "interior leaf" 

#FIXME: doesn't have the right leaf set, apply a permutation to fix this! 

G_2_4_1 = graph_from_edges(Directed, [[8,1],[4,8],[5,8],[5,9],[9,10],[9,4],[4,6],[5,7],[7,6],[6,2],[7,3]] ) 

# unbalanced hexagon 

G_2_4_2 = graph_from_edges(Directed,[[]])
=# 