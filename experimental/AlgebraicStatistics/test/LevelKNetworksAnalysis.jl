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
    println(my_edges)
    N = phylogenetic_network(my_graph)

    return cavender_farris_neyman_model(N)
end

function degree_two_component_stats(M)
    stats = deserialize("stats")
    net = M.phylo_model.graph
    if !haskey(stats, net)
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

        stats[net] = [I, dimension, I_degree, is_I_prime, phi] 
        serialize("stats", stats)
    end
    return stats[net]   
end

function print_stats(M)
    stats = degree_two_component_stats(M)
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

    phi = parametrization(M_2)
    R, y = parameter_ring(M_2)
    S, x = model_ring(M_2)

    M_points = [[QQ(rand(-100:100)) for _ in 1:ngens(R)] for _ in 1:sample]

    points = [[evaluate(im, p) for im in [phi(x) for x in gens(domain(phi))]]  for p in M_points]

    return any(evaluate(f, p) != 0 for p in points)
end

function check_polynomials(I, M_2)
    any(evaluate_at_points(f, M_2) for f in generators(I))
end

function compare_two_networks(M_1, M_2)
    # first one means left inclusion, the other means right inclusion
    result = [false, false]

    stats_1 = degree_two_component_stats(M_1)
    stats_2 = degree_two_component_stats(M_2)

    if stats_1 == nothing && stats_2 == nothing
        println("These networks are NOT distinguishable")
        return false
    end

    if stats_1 == nothing || stats_2 == nothing
        println("These networks are distinguishable!")
        return true
    end

    I_1, dim_1, deg_1, prime_1, phi_1 = stats_1
    I_2, dim_2, deg_2, prime_2, phi_2 = stats_2
    println("Stats calculated!")

    if dim_1 > dim_2
        println("First dimension bigger")
        result[1] = true
    end 

    if dim_2 > dim_1
        println("Second dimension bigger")
        result[2] = true
    end

    result[1] = check_polynomials(I_1, M_2)
    if result[1] == false
        println("These networks are NOT distinguishable!\n")
        return false
    end

    result[2] = check_polynomials(I_2, M_1)
    if result[2] == true
       println("These networks are distinguishable!")
       return true
    else
        println("These networks are NOT distinguishable!\n")
        return false
    end
end

function compare_networks(M)
    result = true
    for (i,m) in enumerate(M)
        println("Network no. ", i)
        print_stats(m)
    end

    for i in 1:length(M)-1
        for j in i+1:length(M)
            print("I compare net ", i, " with net ", j, "\n")
            if !compare_two_networks(M[i], M[j])
                result = false
            end
        end
    end
    if result
        println("All networks were distinguishable!")
    end
end

M_1  = cfn_model_from_the_internet_graph([3 1 
4 6 
11 4 
11 3 
6 5 
3 5 
4 2 
1 2 
6 10 
5 8 
2 9 
1 7])

M_2 = cfn_model_from_the_internet_graph([7 8 
6 8 
3 7 
6 2 
10 6 
7 10 
10 11 
8 9 
3 5 
2 4 
1 3 
1 2])

M_3 = cfn_model_from_the_internet_graph([5 6 
6 4 
6 7 
3 5 
5 4 
4 2 
2 1 
2 8 
3 10 
1 9 
11 3 
11 1])

compare_two_networks(M_1, M_2)

compare_networks([M_1, M_2, M_3])

