@testset "PhylogeneticModels" begin # Changed from "Graphical Models tests"
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
    @test Oscar.n_states(model) == 4
    @test_skip @test Oscar.n_states(general_markov_model(tree; number_states = 20)) == 20
    @test model.root_distribution[1] isa Symbol
    @test typeof(entry_root_distribution(model, 2)) <: MPolyRingElem
    @test allunique([entry_root_distribution(model, i) for i in 1:4])
    @test allunique([Oscar.entry_transition_matrix(model, i, j, Edge(4,3)) for i in 1:4 for j in 1:4])

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

  R, q = Oscar.model_ring(model)
  @testset "Coordinate change Prob - Fourier" begin
    f = coordinate_change(model)
    f1 = inverse_coordinate_change(model)
    @test f1(f(q[1, 2, 2])) == q[1, 2, 2]
    @test f1(f(q[2, 3, 4])) == q[2, 3, 4]
  end

  @testset "vanishing ideal" begin
    vanishing_ideal(model) == ideal( -q[2, 3, 4]^2*q[1, 1, 1] + q[2, 2, 1]*q[2, 1, 2]*q[1, 2 , 2])
  end

  @test_skip @testset "Affine parametrization" begin
    #Probability map
    p = probability_map(affine_phylogenetic_model!(model))
    @test sum(collect(values(p))) == 1

    #Fourier map
    q = fourier_map(affine_phylogenetic_model!(model))
    @test q[1,1,1] == 1

    # GMM model
    p = probability_map(affine_phylogenetic_model!(general_markov_model(tree)))
    @test sum(collect(values(p))) == 1
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
end
