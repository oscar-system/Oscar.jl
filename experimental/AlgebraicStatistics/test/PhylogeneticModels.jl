@testset "PhylogeneticModels" begin # Changed from "Graphical Models tests"
  mktempdir() do path
    tree = graph_from_edges(Directed,[[4,1],[4,2],[4,3]])

    @testset "CFN model" begin
      model = cavender_farris_neyman_model(tree)
      @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
      @test model isa GroupBasedPhylogeneticModel
      @test phylogenetic_model(model) isa PhylogeneticModel

      # Model parameters PHYLO MODEL
      @test n_states(model) == 2
      @test [entry_root_distribution(model, i) for i in 1:2] == [1//2,1//2]
      for i in 1:2, j in 1:2
        @test entry_transition_matrix(
          model, 1, i == j ? 1 : 2, Edge(4, 2),) == entry_transition_matrix(model, i, j, 4, 2)
      end

      # Model parameters G-B PHYLO MODEL
      for e in edges(graph(model))
        @test allunique([entry_fourier_parameter(model, i, e) for i in 1:2])
      end
      G = group(model)
      @test is_isomorphic(unique(parent.(G))[1], abelian_group(2))

      # Polynomial rings
      @test ngens(parameter_ring(model)[1]) ==  
        ngens(parameter_ring(phylogenetic_model(model))[1]) == 2(n_edges(tree))
      @test ngens(full_model_ring(model)[1]) ==
        ngens(full_model_ring(phylogenetic_model(model))[1]) ==  n_states(model)^(n_leaves(tree))

      test_save_load_roundtrip(path, model) do loaded
        @test ngens(full_model_ring(loaded)[1]) ==
          ngens(full_model_ring(phylogenetic_model(loaded))[1]) ==  n_states(loaded)^(n_leaves(tree))
      end
    end

    @testset "Jukes Cantor" begin
      model = jukes_cantor_model(tree)
      @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
      @test model isa GroupBasedPhylogeneticModel
      @test phylogenetic_model(model) isa PhylogeneticModel

      # Model parameters PHYLO MODEL
      @test n_states(model) == 4
      @test [entry_root_distribution(model, i) for i in 1:4] == [1//4,1//4,1//4,1//4]
      for i in 1:4, j in 1:4
        @test entry_transition_matrix(
          model, 1, i == j ? 1 : 2, Edge(4, 2),) == entry_transition_matrix(model, i, j, 4, 2)
      end

      # Model parameters G-B PHYLO MODEL
      for e in edges(graph(model))
        @test entry_fourier_parameter(model, 2, e) ==  entry_fourier_parameter(model, 3, e) ==
          entry_fourier_parameter(model, 4, src(e), dst(e))
      end
      G = group(model)
      @test is_isomorphic(unique(parent.(G))[1], abelian_group(2,2))

      # Polynomial rings
      @test ngens(parameter_ring(model)[1]) ==  
        ngens(parameter_ring(phylogenetic_model(model))[1]) == 2(n_edges(tree))
      @test ngens(full_model_ring(model)[1]) ==
        ngens(full_model_ring(phylogenetic_model(model))[1]) ==  n_states(model)^(n_leaves(tree))
    end

    @testset "Kimura 2" begin
      model = kimura2_model(tree)
      @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
      @test model isa GroupBasedPhylogeneticModel
      @test phylogenetic_model(model) isa PhylogeneticModel

      # Model parameters PHYLO MODEL
      @test n_states(model) == 4
      @test [entry_root_distribution(model, i) for i in 1:4] == [1//4,1//4,1//4,1//4]
      for i in 1:4, j in 1:4
        @test entry_transition_matrix(
          model, 1, (i == j ? 1 : (isodd(i+j) ? 2 : 3)), Edge(4, 2)) == 
            entry_transition_matrix(model, i, j, 4, 2)
      end

      # Model parameters G-B PHYLO MODEL
      for e in edges(graph(model))
        @test entry_fourier_parameter(model, 3, e) ==
          entry_fourier_parameter(model, 4, src(e), dst(e))
      end
      G = group(model)
      @test is_isomorphic(unique(parent.(G))[1], abelian_group(2,2))

      # Polynomial rings
      @test ngens(parameter_ring(model)[1]) ==  
        ngens(parameter_ring(phylogenetic_model(model))[1]) == 3(n_edges(tree))
      @test ngens(full_model_ring(model)[1]) ==
        ngens(full_model_ring(phylogenetic_model(model))[1]) ==  n_states(model)^(n_leaves(tree))
    end

    @testset "Kimura 3" begin
      model = kimura3_model(tree)
      @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
      @test model isa GroupBasedPhylogeneticModel
      @test phylogenetic_model(model) isa PhylogeneticModel

      # Model parameters PHYLO MODEL
      @test n_states(model) == 4
      @test [entry_root_distribution(model, i) for i in 1:4] == [1//4,1//4,1//4,1//4]
      e = Edge(4,2)
      @test entry_transition_matrix(model, 1, 1, 4, 2) == entry_transition_matrix(model, 2, 2, e) == 
        entry_transition_matrix(model, 3, 3, e) == entry_transition_matrix(model, 4, 4, e)
      @test entry_transition_matrix(model, 2, 1, 4, 2) == entry_transition_matrix(model, 1, 2, e) == 
        entry_transition_matrix(model, 4, 3, e) == entry_transition_matrix(model, 3, 4, e)
      @test entry_transition_matrix(model, 1, 3, 4, 2) == entry_transition_matrix(model, 3, 1, e) == 
        entry_transition_matrix(model, 2, 4, e) == entry_transition_matrix(model, 4, 2, e)
      @test entry_transition_matrix(model, 1, 4, 4, 2) == entry_transition_matrix(model, 4, 1, e) == 
        entry_transition_matrix(model, 2, 3, e) == entry_transition_matrix(model, 3, 2, e)
      

      # Model parameters G-B PHYLO MODEL
      for e in edges(graph(model))
        @test allunique([entry_fourier_parameter(model, i, e) for i in 1:4])
      end

      G = group(model)
      @test is_isomorphic(unique(parent.(G))[1], abelian_group(2,2))

      # Polynomial rings
      @test ngens(parameter_ring(model)[1]) ==  
        ngens(parameter_ring(phylogenetic_model(model))[1]) == 4(n_edges(tree))
      @test ngens(full_model_ring(model)[1]) ==
        ngens(full_model_ring(phylogenetic_model(model))[1]) ==  n_states(model)^(n_leaves(tree))
    end

    @testset "general_markov_model" begin
      model = general_markov_model(tree)         
      @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
      @test model isa PhylogeneticModel

      # Model parameters PHYLO MODEL
      @test n_states(model) == 4
      # @test_skip @test Oscar.n_states(general_markov_model(tree; number_states = 20)) == 20
      @test model.root_distribution[1] isa Symbol
      @test typeof(entry_root_distribution(model, 2)) <: MPolyRingElem
      @test allunique([entry_root_distribution(model, i) for i in 1:4])
      @test allunique([entry_transition_matrix(model, i, j, Edge(4,3)) for i in 1:4 for j in 1:4])

      # Polynomial rings
      @test ngens(parameter_ring(model)[1]) ==  4 + 16(n_edges(tree))
      @test ngens(full_model_ring(model)[1]) == n_states(model)^(n_leaves(tree))
    end

    # Test parametrizations for a specific tree and model
    model = jukes_cantor_model(tree)

    @testset "parametrization PhylogeneticModel" begin
      pm = phylogenetic_model(model)

      R, p = full_model_ring(pm)
      f = full_parametrization(pm)

      @test codomain(f) == parameter_ring(pm)[1]
      @test domain(f) == Oscar._ring(R)

      @test length(f.img_gens) == n_states(model)^n_leaves(tree)


      # Test some entries of the joint distribution vector
      a = gens(parameter_ring(pm)[1])[1:3]
      b = gens(parameter_ring(pm)[1])[4:6]
      @test f(p[1,1,1]) == 1//4*a[1]*a[2]*a[3] + 3//4*b[1]*b[2]*b[3]
      @test f(p[1,2,3]) == 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//4*b[1]*b[2]*b[3]
      @test f(p[4,4,3]) == 1//4*a[1]*a[2]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//2*b[1]*b[2]*b[3]

      @testset "Equivalent classes probabilities" begin
        eq_classes = equivalent_classes(pm)
        @test length(eq_classes) == 5
        @test issetequal(collect(keys(eq_classes)),
                         [(1, 2, 1), (1, 1, 1), (1, 2, 2), (1, 2, 3), (1, 1, 2)])
        @test length(eq_classes[1, 1, 1]) == 4
        @test length(eq_classes[1, 1, 2]) == 12
        @test length(eq_classes[1, 2, 3]) == 24
        @test length(eq_classes[1, 2, 2]) == 12
        @test length(eq_classes[1, 2, 1]) == 12

        @test issetequal(eq_classes[1, 1, 1],
                         [p[1, 1, 1], p[2, 2, 2], p[3, 3, 3], p[4, 4, 4]])
      end
    end

    @testset "parametrization GroupBasedPhylogeneticModel" begin

      R, q = full_model_ring(model)
      f = full_parametrization(model)

      @test codomain(f) == parameter_ring(model)[1]
      @test domain(f) == Oscar._ring(R)

      @test length(f.img_gens) == n_states(model)^n_leaves(tree)

      
      # Test some entries of the joint distribution vector
      x = gens(parameter_ring(model)[1])[1:3]
      y = gens(parameter_ring(model)[1])[4:6]
      @test f(q[3,2,4]) == y[1]*y[2]*y[3]
      @test f(q[2,1,2]) == x[2]*y[1]*y[3]
      @test f(q[1,1,1]) == x[1]*x[2]*x[3]


      @testset "Equivalent classes fourier" begin
        eq_classes = equivalent_classes(model)
        @test length(eq_classes) == 5
        @test issetequal(collect(keys(eq_classes)),
                         [(1, 1, 1), (2, 2, 1), (2, 1, 2), (2, 3, 4), (1, 2, 2)])
        @test length(eq_classes[1, 1, 1]) == 1
        @test length(eq_classes[2, 2, 1]) == 3
        @test length(eq_classes[2, 1, 2]) == 3
        @test length(eq_classes[2, 3, 4]) == 6
        @test length(eq_classes[1, 2, 2]) == 3

        @test issetequal(eq_classes[2, 2, 1],
                         [q[2, 2, 1], q[3, 3, 1], q[4, 4, 1]])
      end
    end

    R, q = model_ring(model)
    @testset "Coordinate change Prob - Fourier" begin
      f = coordinate_change(model)
      f1 = inverse_coordinate_change(model)
      @test f1(f(q[1, 2, 2])) == q[1, 2, 2]
      @test f1(f(q[2, 3, 4])) == q[2, 3, 4]
    end

    @testset "vanishing ideal" begin
      vanishing_ideal(model) == ideal( -q[2, 3, 4]^2*q[1, 1, 1] + q[2, 2, 1]*q[2, 1, 2]*q[1, 2 , 2])
    end

    @testset "Affine parametrization" begin
      p = full_affine_parametrization(phylogenetic_model(model))
      @test sum(p.img_gens) == 1

      f = affine_parametrization(model)
      @test f(q[1,1,1]) == 1
    end

    @testset "Transition Matrix Types" begin
      tree = graph_from_edges(Directed,[[4,1],[4,2],[4,3]])
      @testset "M = VarName - T = FieldElem" begin
        m = jukes_cantor_model(tree)
        parameter_ring(phylogenetic_model(m))
        entry_transition_matrix(m, 1, 3, 4, 1)
      end
      @testset  "M = VarName - T = VarName" begin
        m = general_markov_model(tree)
        parameter_ring(m)
        entry_transition_matrix(m, 1, 3, 4, 1)
      end

      @testset "M = MPolyRingElem - T = MPolyRingElem" begin
        Rp, p = polynomial_ring(QQ, :p => 1:4)
        RM, (a, b, c, d) = polynomial_ring(Rp, [:a, :b, :c, :d])
        p = [p[1], p[2], p[3], p[4]]
        M = [a * p[1] c * p[2] b * p[3] b * p[4];
             c * p[1] a * p[2] b * p[3] b * p[4];
             b * p[1] b * p[2] a * p[3] d * p[4];
             b * p[1] b * p[2] d * p[3] a * p[4]]
        PM = PhylogeneticModel(tree, M, p)
        parameter_ring(PM)
        entry_transition_matrix(PM, 1, 3, 4, 1)

      end

      @testset "M = MPolyRingElem - T = FieldElem" begin 
        RM, (a, b, c, d) = polynomial_ring(QQ, [:a, :b, :c, :d])
        p = QQ.([5//8, 1//8, 1//8, 1//8])
        M = [a * p[1] c * p[2] b * p[3] b * p[4];
             c * p[1] a * p[2] b * p[3] b * p[4];
             b * p[1] b * p[2] a * p[3] d * p[4];
             b * p[1] b * p[2] d * p[3] a * p[4]]
        PM = PhylogeneticModel(tree, M, p)
        parameter_ring(PM)
        entry_transition_matrix(PM, 1, 3, 4, 1)
      end

      @testset "M = MPolyRingElem - T = RationalFunctionFieldElem" begin
        Rp, p = rational_function_field(QQ, :p => 1:4)
        RM, (lambda, mu) = polynomial_ring(Rp, [:lambda, :mu])
        M = [lambda * p[1]                              lambda * p[2]                              lambda * p[1] + mu * p[1] // (p[1] + p[3]) lambda * p[4];
             lambda * p[1]                              lambda * p[2]                              lambda * p[3] lambda * p[4] + mu * p[4] // (p[2] + p[4]);
             lambda * p[1] + mu * p[1] // (p[1] + p[3]) lambda * p[2]                              lambda * p[3]                              lambda * p[4];
             lambda * p[1]                              lambda * p[2] + mu * p[2] // (p[2] + p[4]) lambda * p[3]                              lambda * p[4]]
        
        PM = PhylogeneticModel(tree, M, p)
        parameter_ring(PM)
        entry_transition_matrix(PM, 1, 3, 4, 1)
      end
    end

    @testset "Networks" begin
      G1 = graph_from_edges(Directed,[[5,1], [6,5], [7,6], [7,8], [8,5], [6,2], [7,3], [8,4]]);
      N1 = phylogenetic_network(G1)

      @test level(N1) == 1
      @test n_hybrid(N1) == 1
      @test Set(hybrid_edges(N1)[1]) == Set([Edge(6, 5), Edge(8, 5)])

      G2 = graph_from_edges(Directed,[[7,2], [8,7], [7,6], [8,6], [9,8], [10,9], [11,10], [11,12], [12,10], [7,2], [6,1], [9,3], [11,4], [12,5]]);
      N2 = phylogenetic_network(G2)

      @test level(N2) == 1
      @test n_hybrid(N2) == 2
      @test Set(hybrids(N2)[6]) == Set([Edge(7, 6), Edge(8, 6)])
      @test Set(hybrids(N2)[10]) == Set([Edge(11, 10), Edge(12, 10)])
      
      G3 = graph_from_edges(Directed,[[6,2], [6,7], [6,5], [7,5], [5,4], [7,8], [8,4], [4,1], [8,3]]);
      N3 = phylogenetic_network(G3)
      @test level(N3) == 2
      @test n_hybrid(N3) == 2
      @test issetequal(hybrid_vertices(N3), [4,5])

      PM1 = cavender_farris_neyman_model(N1)
      PM2 = cavender_farris_neyman_model(N2)
      PM3 = cavender_farris_neyman_model(N3)
      
      @test n_generators(parameter_ring(PM1)[1]) == n_edges(N1) * 2 + n_hybrid(N1)*2
      @test n_generators(parameter_ring(PM2)[1]) == n_edges(N2) * 2 + n_hybrid(N2)*2
      @test n_generators(parameter_ring(PM3)[1]) == n_edges(N3) * 2 + n_hybrid(N3)*2

      @test ngens(full_model_ring(PM1)[1]) == 2^length(leaves(N1))
      @test ngens(full_model_ring(PM2)[1]) == 2^length(leaves(N2))
      @test ngens(full_model_ring(PM3)[1]) == 2^length(leaves(N3))

      @testset "parametrizations Networks" begin
        # N1
        SS = parameter_ring(PM1)
        RR = model_ring(PM1)
        q = RR[2]
        fN1 = parametrization(PM1)

        T1 = graph_from_edges(Directed,[[5,1], [6,5], [7,6], [7,8], [6,2], [7,3], [8,4]]);
        T2 = graph_from_edges(Directed,[[5,1], [7,6], [7,8], [8,5], [6,2], [7,3], [8,4]]);

        PMT1 = cavender_farris_neyman_model(T1)
        PMT2 = cavender_farris_neyman_model(T2)

        set_attribute!(PMT1, :parameter_ring, SS[1:2])
        set_attribute!(PMT2, :parameter_ring, SS[1:2])

        set_attribute!(PMT1, :model_ring, RR)
        set_attribute!(PMT2, :model_ring, RR)

        f1 = parametrization(PMT1)
        f2 = parametrization(PMT2)

        @test fN1(q[1, 1, 1, 1]) == SS[3][Edge(6,5)]*f1(q[1, 1, 1, 1]) + SS[3][Edge(8,5)]*f2(q[1, 1, 1, 1]) 
        @test fN1(q[1, 2, 1, 2]) == SS[3][Edge(6,5)]*f1(q[1, 2, 1, 2]) + SS[3][Edge(8,5)]*f2(q[1, 2, 1, 2]) 

        # N2
        SS = parameter_ring(PM2)
        RR = model_ring(PM2)
        q = RR[2]
        fN2 = parametrization(PM2)

        T1 = graph_from_edges(Directed,[[7,2], [8,7], [7,6], [9,8], [10,9], [11,10], [11,12], [7,2], [6,1], [9,3], [11,4], [12,5]]);
        T2 = graph_from_edges(Directed,[[7,2], [8,7], [7,6], [9,8], [10,9], [11,12], [12,10], [7,2], [6,1], [9,3], [11,4], [12,5]]);
        T3 = graph_from_edges(Directed,[[7,2], [8,7], [8,6], [9,8], [10,9], [11,10], [11,12], [7,2], [6,1], [9,3], [11,4], [12,5]]);
        T4 = graph_from_edges(Directed,[[7,2], [8,7], [8,6], [9,8], [10,9], [11,12], [12,10], [7,2], [6,1], [9,3], [11,4], [12,5]]);

        PMT1 = cavender_farris_neyman_model(T1)
        PMT2 = cavender_farris_neyman_model(T2)
        PMT3 = cavender_farris_neyman_model(T3)
        PMT4 = cavender_farris_neyman_model(T4)

        set_attribute!(PMT1, :parameter_ring, SS[1:2])
        set_attribute!(PMT2, :parameter_ring, SS[1:2])
        set_attribute!(PMT3, :parameter_ring, SS[1:2])
        set_attribute!(PMT4, :parameter_ring, SS[1:2])

        set_attribute!(PMT1, :model_ring, RR)
        set_attribute!(PMT2, :model_ring, RR)
        set_attribute!(PMT3, :model_ring, RR)
        set_attribute!(PMT4, :model_ring, RR)

        f1 = parametrization(PMT1)
        f2 = parametrization(PMT2)
        f3 = parametrization(PMT3)
        f4 = parametrization(PMT4)

        l11 = SS[3][Edge(7,6)]
        l12 = SS[3][Edge(8,6)]
        l21 = SS[3][Edge(11,10)]
        l22 = SS[3][Edge(12,10)]

        @test fN2(q[1, 1, 1, 1, 1]) == l11*l21*f1(q[1, 1, 1, 1, 1]) + l11*l22*f2(q[1, 1, 1, 1, 1]) + l12*l21*f3(q[1, 1, 1, 1, 1]) + l12*l22*f4(q[1, 1, 1, 1, 1]) 
        @test fN2(q[1, 1, 2, 2, 1]) == l11*l21*f1(q[1, 1, 2, 2, 1]) + l11*l22*f2(q[1, 1, 2, 2, 1]) + l12*l21*f3(q[1, 1, 2, 2, 1]) + l12*l22*f4(q[1, 1, 2, 2, 1]) 

      end

    end
  end
end
