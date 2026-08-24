# here are functions relating to the enumeration of level-2 (hopefully level-k in general) phylogenetic networks.

# given number of leaves n generates the possible simple level-2 network topologies, which are uniquely given by partitions of n.
# note that n must be at least 2
function lvl2_leaf_partition(n::Int64)

    partitions = []

    #get partitions of the form [a,b,0]

    for i in 1:(div(n,2))

        push!(partitions, [n-i, i, 0])
    end
    #get partitions of the form [a,b,c]

    for i in 0:(n - 2 - Int64(ceil(n/3)))
        for j in 0:(div(i,2))

            if  (3 + 2*i - n) > j
                continue
            end

            push!(partitions, [n-(2+i), 1+i-j, 1+j])
        end
    end

    return partitions
end

# function for creating the (undirected) graph corresponding to network topology given by partition (a,b,c)

function graph_from_partition(v::Vector{Int64})

    path_1 = [[i,i+1] for i in 1 : v[1]+1]
    path_2 = unique!(vcat([[1,(v[1]+3)],[v[2]+v[1]+2, (v[1]+2)]], [[i,i+1] for i in (v[1]+3) : v[1]+v[2]+1]))

    if v[3] != 0

        path_3 = unique!(vcat([[1,(v[2]+v[1]+3)],[v[1]+v[2]+v[3]+2, (v[1]+2)]], [[i,i+1] for i in (v[1]+v[2]+1)+3 : v[1]+v[2]+v[3]+1]))

    else

        path_3 = [1,(v[1]+2)]

    end

    g = graph_from_edges(vcat(path_1, path_2, path_3))

    return g

end

#given partition produces all symmetries of the network topology
function list_of_perm(v::Vector{Int64})
    
    #useful for vertical reflection symmetry

    no_inner_vertex = sum(v)+2
    G = symmetric_group(no_inner_vertex)
    perm_list = [G([i for i in 1:no_inner_vertex])]

    check1 = false
    check2 = false

    if v[1] == v[2]

        check1 = true
        #this function requires the Permutations.jl package
        # Permutation(perm_vec1)
        perm_vec1 = [ vcat([1], [i for i in (v[1]+3):(v[1]+v[2]+2)], [v[1]+2], [i for i in 2:(v[1]+1)]) ]
        push!(perm_list, G( perm_vec1[1] ) )

    end

    if v[2]==v[3]

        check2 = true
        perm_vec2 = [ vcat([i for i in 1:(v[1]+2)], [i for i in v[1]+v[2]+3:v[1]+v[2]+v[3]+2], [i for i in v[1]+3:v[1]+v[2]+2]) ]
        push!(perm_list, G( perm_vec2[1] ) )

    end

    if check1==true && check2==true
        
        push!(perm_list, perm_list[1]*perm_list[2])
        
    end

    perm_vec3 = [ vcat( 
                [v[1]+2], 
                reverse([i for i in 2:(v[1]+1)]), 
                [1], 
                reverse([i for i in (v[1]+3):(v[1]+v[2]+2)]), 
                reverse([i for i in (v[1]+v[2]+3):(v[1]+v[2]+v[3]+2)])
                )]
    push!(perm_list, G(perm_vec3[1])) 

    return perm_list

end

#creates all possible pairs of inner vertices for reticulation and mods out by symmetries of network topology
function hybrid_vertex_mod_sym(v::Vector{Int64})

    perms = list_of_perm(v)
    my_list = []

    p1 = [i for i in 2:v[1]+1]
    p2 = [i for i in v[1]+3:v[1]+v[2]+2]
    p3 = [i for i in v[1]+v[2]+3:v[1]+v[2]+v[3]+2]

    paths = [p1,p2,p3]

    n = sum(v)+2
    my_list = []

    for i in 1:n-1 
        for j in i+1:n

            check = false

            if i == v[1] + 2 || j == v[1] + 2
                check = true
            end

            for p in paths

                if i in p && j in p

                    check = true
                    continue
                end
            end


            if !check
                push!(my_list, [i,j])
            end
        end
    end

    #h_pair = [1,3]

    # for tau in perms
    #     push!(deletion_list, [tau(h_pair[1]) tau(h_pair[2])])
    # end

    rep_list = []
    
    for h_pair in my_list

        orbit_list = []
        min_candidate_list = []

        for tau in perms
            if tau(h_pair[1]) > tau(h_pair[2])
                push!(orbit_list, [tau(h_pair[2]), tau(h_pair[1])])
            else
                push!(orbit_list, [tau(h_pair[1]), tau(h_pair[2])])
            end
        end

        min1 = minimum([orbit_list[i][1] for i in 1:length(orbit_list)])

        for rep in orbit_list
                    
            if rep[1] == min1

                push!(min_candidate_list, rep)

            end        
        end

        min2 = minimum([min_candidate_list[i][2] for i in 1:length(min_candidate_list)])
        min_candidate_list_val = [min_candidate_list[i][2] for i in 1:length(min_candidate_list) ]

        ind = findfirst(==(min2), min_candidate_list_val)

        push!(rep_list, min_candidate_list[ind])
    end

    return unique!(rep_list)

end

