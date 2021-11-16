@testset "Examples.GITFans, the example" begin

    using Oscar
    using Oscar.GITFans

    Q = [
     1  1   0   0   0 ;
     1  0   1   1   0 ;
     1  0   1   0   1 ;
     1  0   0   1   1 ;
     0  1   0   0  -1 ;
     0  1   0  -1   0 ;
     0  1  -1   0   0 ;
     0  0   1   0   0 ;
     0  0   0   1   0 ;
     0  0   0   0   1 ];

    n = size(Q, 1)
    Qt, T = Oscar.PolynomialRing(Oscar.QQ, :T => 1:n)
    D = free_abelian_group(size(Q,2))
    w = [D(Q[i, :]) for i = 1:n]
    R = grade(Qt, w)
    a = ideal([
        T[5]*T[10] - T[6]*T[9] + T[7]*T[8],
        T[1]*T[9]  - T[2]*T[7] + T[4]*T[5],
        T[1]*T[8]  - T[2]*T[6] + T[3]*T[5],
        T[1]*T[10] - T[3]*T[7] + T[4]*T[6],
        T[2]*T[10] - T[3]*T[9] + T[4]*T[8],
    ])

    perms_list = [ [1,3,2,4,6,5,7,8,10,9], [5,7,1,6,9,2,8,4,10,3] ];
    sym10 = symmetric_group(n);
    G, emb = sub([sym10(x) for x in perms_list]...);

    fanobj = GITFans.git_fan(a, Q, G)
    @test fanobj.F_VECTOR == [20, 110, 240, 225, 76]

    collector_cones = GITFans.orbit_cones(a, Q, G)
    matrix_action = GITFans.action_on_target(Q, G)
    orbit_list = GITFans.orbit_cone_orbits(collector_cones, matrix_action)
    @test map(length, orbit_list) == [10, 15, 10, 1]

    perm_actions = GITFans.action_on_orbit_cone_orbits(orbit_list, matrix_action)
    q_cone = Polymake.polytope.Cone(INPUT_RAYS = Q)

    (hash_list, edges) = GITFans.fan_traversal(orbit_list, q_cone, perm_actions)
    @test length(hash_list) == 6
    @test hash_list[1] == [
           BitSet([]),
           BitSet([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]),
           BitSet([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
           BitSet([1])]

    intergraph = Polymake.graph.graph_from_edges(collect(edges));

    expanded = GITFans.orbits_of_maximal_GIT_cones(orbit_list, hash_list, matrix_action);
    orbit_lengths = map(length, expanded)
    @test sum(orbit_lengths) == 76

    maxcones = vcat( expanded... );
    full_edges = GITFans.edges_intersection_graph(maxcones, size(Q, 2) - 1);
    @test length(full_edges) == 180

    full_intergraph = Polymake.graph.graph_from_edges(collect(full_edges));

    fanobj = GITFans.hashes_to_polyhedral_fan(orbit_list, hash_list, matrix_action)
    @test fanobj.F_VECTOR == [20, 110, 240, 225, 76]

#   # Now try the construction with trivial symmetry group.
#   G2 = trivial_subgroup(G)[1]
#   fanobj2 = GITFans.git_fan(a, Q, G2)
#   @test fanobj2.F_VECTOR == [20, 110, 240, 225, 76]
#T This call would fail with `Segmentation fault`.

end

@testset "Examples.GITFans, is_monomial_free" begin

    # The following input gave a wrong `true` result
    # with an earlier version of `is_monomial_free`.
    using Oscar
    using Oscar.GITFans

    nr_variables = 4
    vars_strings = map( i -> "x"*string(i), 1:nr_variables )
    R, T = PolynomialRing(QQ,vars_strings)
    ideal_gens = [
        T[1]^2*T[2]*T[3]^2*T[4]^2+T[1]^2*T[2]^2*T[3],
        2*T[1]*T[2]^2*T[3]^2*T[4]+T[1]*T[3]^2*T[4]^2-T[1]^2,
        T[1]*T[2]*T[3]*T[4]-T[1]^2*T[4]^2+T[1]*T[2]-T[4],
        T[1]*T[2]*T[3]-T[2]*T[3]*T[4],
        T[1]^2*T[2]^2*T[4]^2-T[1]^2*T[3]^2*T[4]^2+T[1]*T[2]*T[3]*T[4]
    ]

    @test ! GITFans.is_monomial_free(ideal(ideal_gens))
    
end
