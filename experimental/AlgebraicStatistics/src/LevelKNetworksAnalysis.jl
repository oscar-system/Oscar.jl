"""This part serves for drawing all 2-level networks with n leaves"""

function draw_path(p)
    return [[p[i], p[i+1]] for i in 1:length(p)-1]
end

function add_leaves(sum, path)
    path = path[2: length(path) - 1]
    edges = []
    for i in 1:length(path)
        push!(edges, [path[i], sum + i])
    end
    return edges
end

function draw_network(ret, par)
    edges = Vector{Vector{Int}}()
    v_1 = ret[1]
    v_2 = ret[2]

    if v_1 > v_2
        temp = v_2
        v_2 = v_1
        v_1 = temp
    end

    p_1 = [i for i in 1:par[1] + 2]
    p_2 = vcat([1], [i for i in par[1] + 3: par[1] + par[2] + 2], [par[1] + 2]) 
    p_3 = vcat([1], [i for i in par[1] + par[2] + 3: par[1] + par[2] + par[3]+ 2], [par[1] + 2])
    paths = [p_1, p_2, p_3]

    if v_1 == 1
        q_1, q_2, q_3 = [], [], []
        if v_2 in p_1
            q_1 = p_3
            q_2 = p_2
            q_3 = p_1
        end
        if v_2 in p_2
            q_1 = p_1
            q_2 = p_3
            q_3 = p_2
        end
        if v_2 in p_3
            q_1 = p_1
            q_2 = p_2
            q_3 = p_3
        end
        append!(edges, draw_path(q_1))
        append!(edges, draw_path(q_2))
        i = findfirst(==(v_2), q_3)
        A1 = vcat(q_3[1:i-1], [v_2])
        A2 = q_3[i:end]
        append!(edges, draw_path(A1))
        append!(edges, draw_path(reverse(A2)))
    else
        if (v_1 in p_1) && (v_2 in p_2)
            q_1 = p_3
            q_2 = p_1
            q_3 = p_2
        end

        if (v_1 in p_1) && (v_2 in p_3)
            q_1 = p_2
            q_2 = p_1
            q_3 = p_3
        end

        if (v_1 in p_2) && (v_2 in p_3)
            q_1 = p_1
            q_2 = p_2
            q_3 = p_3
        end

        append!(edges, draw_path(q_1))

        i = findfirst(==(v_1), q_2)
        A1 = vcat(q_2[1:i-1], [v_1])
        A2 = q_2[i:end]
        append!(edges, draw_path(A1))
        append!(edges, draw_path(reverse(A2)))

        j = findfirst(==(v_2), q_3)
        B1 = vcat(q_3[1:j-1], [v_2])
        B2 = q_3[j:end]
        append!(edges, draw_path(B1))
        append!(edges, draw_path(reverse(B2)))
    end

    sum = par[1] + par[2] + par[3] + 2

    for i in 1:3
        append!(edges, add_leaves(sum, paths[i]))
        sum = sum + par[i]
    end
    g = graph_from_edges(Directed, edges)
    n = phylogenetic_network(g)
    M = cavender_farris_neyman_model(n)

    return M
end

draw_network([5,6], [2,1,1])

"""here are functions relating to the enumeration of level-2 (hopefully level-k in general) phylogenetic networks.

given number of leaves n generates the possible simple level-2 network topologies, which are uniquely given by partitions of n.
note that n must be at least 2"""

function lvl2_leaf_partition(n::Int64)

    partitions = Vector{Vector{Int}}()

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

"""This part serves for comparing two networks"""
function pair_list(numbers)
    numbers = vec(numbers')
    if length(numbers) % 2 != 0
        error("The length of the input array must be even.")
        end

    return [[numbers[i], numbers[i + 1]] for i in 1:2:length(numbers)]
end

function cfn_model_from_the_internet_graph(numbers)
    my_edges = pair_list(numbers)
    my_graph = graph_from_edges(Directed, my_edges)
    N = phylogenetic_network(my_graph)

    return cavender_farris_neyman_model(N)
end

function degree_two_component_stats(M, name)
    stats = deserialize("dir_project_data/stats")
    net = M.phylo_model.graph.graph
    if !haskey(stats, name)
        phi = parametrization(M)
        H = components_of_kernel(2, phi, show_progress = true)
        if isempty(H)
            println("0")
            return nothing
        end
        I = Oscar.ideal(reduce(vcat, collect(values(H))))

        dimension = dim(I)
        I_degree = degree(I)
        is_I_prime = is_prime(I)

        Oscar.save("dir_project_data/ideals/"*name, I)

        stats[name] = [name, dimension, I_degree, is_I_prime] 
        serialize("dir_project_data/stats", stats)
    end
    data = stats[name]
    my_ideal = Oscar.load("dir_project_data/ideals/"*data[1])
    data[1] = my_ideal
    return data
end

function print_stats(M, name)
    stats = degree_two_component_stats(M, name)
    if stats == nothing
        println("This network has zero_ideal")
    else
        println("dim(I) = ", stats[2])
        println("deg(I) = ", stats[3])
        println("Is I prime? ", stats[4])
        println("\n\n")
    end
end

function evaluate_at_points(f, M_2)
    sample = 1000
    F = GF(32003)

    phi = parametrization(M_2)
    R, y = parameter_ring(M_2)
    S, x = model_ring(M_2)

    RF, yF = polynomial_ring(F, ["y$i" for i in 1:ngens(R)])
    SF, xF = polynomial_ring(F, ["x$i" for i in 1:ngens(S)])

    M_points = [[rand(F) for _ in 1:ngens(RF)] for _ in 1:sample]

    phi_images = [phi(xi) for xi in gens(domain(phi))]

    phi_F = [RF(im) for im in phi_images]

    points = [[evaluate(im, p) for im in phi_F] for p in M_points]

    f_F = SF(f)

    return any(evaluate(f_F, p) != 0 for p in points)
end

function check_polynomials(I, M_2)
    any(evaluate_at_points(f, M_2) for f in generators(I))
end

function compare_two_networks(net_1, net_2)
    # first one means left inclusion, the other means right inclusion
    result = [false, false]
    M_1 = net_1[1]
    name_1 = net_1[2]
    M_2 = net_2[1]
    name_2 = net_2[2]
    
    stats_1 = degree_two_component_stats(M_1, name_1)
    stats_2 = degree_two_component_stats(M_2, name_2)

    if stats_1 == nothing && stats_2 == nothing
        return false
    end

    if stats_1 == nothing || stats_2 == nothing
        return true
    end

    I_1, dim_1, deg_1, prime_1 = stats_1
    I_2, dim_2, deg_2, prime_2 = stats_2

    if dim_1 > dim_2
        result[1] = true
    end 

    if dim_2 > dim_1
        result[2] = true
    end

    result[1] = check_polynomials(I_1, M_2)
    if result[1] == false
        return false
    end

    result[2] = check_polynomials(I_2, M_1)
    if result[2] == true
       return true
    else
        return false
    end
end

function compare_networks(M)
    result = fill(true, length(M), length(M))
    for (i,m) in enumerate(M)
        print_stats(m[1], m[2])
    end

    for i in 1:length(M)-1
        for j in i+1:length(M)
            println(i, j)
            if !compare_two_networks(M[i], M[j])
                result[i,j] = false
            end
        end
        Oscar.save("dir_project_data/temp_data/calc_leaves_data", result)
    end
    return result
end
