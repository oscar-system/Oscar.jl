@testset "Examples.GITFans, the example" begin

    using Oscar
    using Oscar.GITFans

    Q = [
     1  1   0   0   0
     1  0   1   1   0
     1  0   1   0   1
     1  0   0   1   1
     0  1   0   0  -1
     0  1   0  -1   0
     0  1  -1   0   0
     0  0   1   0   0
     0  0   0   1   0
     0  0   0   0   1
     ]

    n = nrows(Q)
    D = free_abelian_group(ncols(Q))
    w = [D(Q[i, :]) for i = 1:n]
    R, T = graded_polynomial_ring(QQ, :T => 1:n, w)
    a = ideal([
        T[5]*T[10] - T[6]*T[9] + T[7]*T[8],
        T[1]*T[9]  - T[2]*T[7] + T[4]*T[5],
        T[1]*T[8]  - T[2]*T[6] + T[3]*T[5],
        T[1]*T[10] - T[3]*T[7] + T[4]*T[6],
        T[2]*T[10] - T[3]*T[9] + T[4]*T[8],
    ])

    perms_list = [ [1,3,2,4,6,5,7,8,10,9], [5,7,1,6,9,2,8,4,10,3] ];
    sym = symmetric_group(n);
    G, emb = sub([sym(x) for x in perms_list]...);

    fanobj = GITFans.git_fan(a, Q, G)
    @test f_vector(fanobj) == [20, 110, 240, 225, 76]

    collector_cones = GITFans.orbit_cones(a, Q, G)
    matrix_action = GITFans.action_on_target(Q, G)
    orbit_list = GITFans.orbit_cone_orbits(collector_cones, matrix_action)
    @test map(length, orbit_list) == [10, 15, 10, 1]

    perm_actions = GITFans.action_on_orbit_cone_orbits(orbit_list, matrix_action)
    q_cone = positive_hull(Q)

    (hash_list, edges) = GITFans.fan_traversal(orbit_list, q_cone, perm_actions)
    @test length(hash_list) == 6
    @test hash_list[1] == [
           BitSet([]),
           BitSet([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]),
           BitSet([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
           BitSet([1])]

    orbit_graph = Graph{Undirected}(length(hash_list));
    [add_edge!(orbit_graph, e...) for e in edges];
    # The following works in an interactive session.
    # visualize(orbit_graph)

    expanded = GITFans.orbits_of_maximal_GIT_cones(orbit_list, hash_list, matrix_action);
    orbit_lengths = map(length, expanded)
    @test sum(orbit_lengths) == 76

    maxcones = vcat( expanded... );
    full_edges = GITFans.edges_intersection_graph(maxcones, size(Q, 2) - 1);
    @test length(full_edges) == 180

    full_graph = Graph{Undirected}(length(maxcones));
    [add_edge!(full_graph, e...) for e in full_edges];
    # The following works in an interactive session.
    # visualize(full_graph)

    fanobj = GITFans.hashes_to_polyhedral_fan(orbit_list, hash_list, matrix_action)

    # Inspect the fan object.
    @test f_vector(fanobj) == [20, 110, 240, 225, 76]
    c = cones(fanobj, 5)[1]
    @test n_rays(fanobj) == 20
    @test dim(fanobj) == 5
    @test n_maximal_cones(fanobj) == 76
    @test n_cones(fanobj) == 671
    @test !is_complete(fanobj)
    @test is_pointed(fanobj)
    @test !is_regular(fanobj)
    @test !is_simplicial(fanobj)
    @test !is_smooth(fanobj)

#   # Now try the construction with trivial symmetry group (takes longer).
#   G2 = trivial_subgroup(G)[1]
#   fanobj2 = GITFans.git_fan(a, Q, G2)
#   @test f_vector(fanobj2) == [20, 110, 240, 225, 76]

end

@testset "Examples.GITFans, is_monomial_free" begin

    # The following input gave a wrong `true` result
    # with an earlier version of `is_monomial_free`.
    using Oscar
    using Oscar.GITFans

    R, T = polynomial_ring(QQ, 4)
    ideal_gens = [
        T[1]^2*T[2]*T[3]^2*T[4]^2+T[1]^2*T[2]^2*T[3],
        2*T[1]*T[2]^2*T[3]^2*T[4]+T[1]*T[3]^2*T[4]^2-T[1]^2,
        T[1]*T[2]*T[3]*T[4]-T[1]^2*T[4]^2+T[1]*T[2]-T[4],
        T[1]*T[2]*T[3]-T[2]*T[3]*T[4],
        T[1]^2*T[2]^2*T[4]^2-T[1]^2*T[3]^2*T[4]^2+T[1]*T[2]*T[3]*T[4]
    ]

    @test ! GITFans.is_monomial_free(ideal(ideal_gens))
    
end
