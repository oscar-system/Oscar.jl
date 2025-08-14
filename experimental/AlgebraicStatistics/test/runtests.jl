const ColoredGGM{Directed} = GaussianGraphicalModel{
  Directed,
  @NamedTuple{color::S}
} where S <: Oscar.GraphMap;
# 

@testset "GaussianGraphicalModels" begin
  DG = graph_from_edges(Directed, [[1,2],[2,3]])
  @testset "Directed" begin
    M1 = gaussian_graphical_model(DG)
    cov_mat = covariance_matrix(M1)
    V1 = vanishing_ideal(M1)
    @test V1 == ideal(
      [-cov_mat[1, 2] * cov_mat[2, 3] + cov_mat[1, 3] * cov_mat[2, 2]]
    )

    label!(DG,
           Dict((1, 2) => "pink", (2, 3) => "pink"),
           Dict(i => "green" for i in 1:3);
           name=:color)
    M2 = gaussian_graphical_model(DG)
    cov_mat = covariance_matrix(M2)
    V2 = vanishing_ideal(M2)
    @test V2 == ideal(
      [-cov_mat[1, 2] * cov_mat[2, 3] + cov_mat[1, 3] * cov_mat[2, 2]]
    )

  end

  UG = complete_bipartite_graph(1, 2)
  @testset "Undirected" begin
    M3 = gaussian_graphical_model(UG)
    V3 = vanishing_ideal(M3)
    cov_mat = covariance_matrix(M3)
    @test V3 == ideal([
      -cov_mat[1, 1] * cov_mat[2, 3] + cov_mat[1, 2] * cov_mat[1, 3]])
  end

  @testset "Mixed" begin
    MG = mixed_graph_from_edges(collect(edges(DG)), collect(edges(UG)))
    M4 = gaussian_graphical_model(MG)
  end
end

# exported items: experimental/GraphicalModels/src/GraphicalModels.jl

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
    @test ngens(model_ring(model)[1]) ==
          ngens(model_ring(phylogenetic_model(model))[1]) ==  n_states(model)^(n_leaves(tree))
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
    @test ngens(model_ring(model)[1]) ==
          ngens(model_ring(phylogenetic_model(model))[1]) ==  n_states(model)^(n_leaves(tree))
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
    @test ngens(model_ring(model)[1]) ==
          ngens(model_ring(phylogenetic_model(model))[1]) ==  n_states(model)^(n_leaves(tree))
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
    @test ngens(model_ring(model)[1]) ==
          ngens(model_ring(phylogenetic_model(model))[1]) ==  n_states(model)^(n_leaves(tree))
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
    @test ngens(model_ring(model)[1]) == n_states(model)^(n_leaves(tree))
  end

  @test_skip @testset "affine models" begin 
    model = jukes_cantor_model(tree)
    @test Oscar.affine_phylogenetic_model!(model) isa GroupBasedPhylogeneticModel
    #sum of each row of transition matrices should equal one 
    tr_mat = Oscar.transition_matrices(Oscar.affine_phylogenetic_model!(model))
    for j in 1:4
      rowsum = 0
      for i in 1:4
          rowsum = rowsum+tr_mat[Edge(4, 2)][j,i]
      end
      @test rowsum == 1
    end
    #sum of root distribution should equal one
    @test sum(Oscar.affine_phylogenetic_model!(model).phylo_model.root_distr) == 1
  end

  # Test parametrizations for a specific tree and model
  model = jukes_cantor_model(tree)

  @testset "parametrization PhylogeneticModel" begin
    pm = phylogenetic_model(model)
    f = parametrization(pm)

    @test f.codomain == parameter_ring(pm)[1]
    @test f.domain == model_ring(pm)[1]

    @test length(f.img_gens) == n_states(model)^n_leaves(tree)

    _, p = model_ring(pm)

    # Test some entries of the joint distribution vector
    a = gens(parameter_ring(pm)[1])[1:3]
    b = gens(parameter_ring(pm)[1])[4:6]
    @test f(p[1,1,1]) == 1//4*a[1]*a[2]*a[3] + 3//4*b[1]*b[2]*b[3]
    @test f(p[1,2,3]) == 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//4*b[1]*b[2]*b[3]
    @test f(p[4,4,3]) == 1//4*a[1]*a[2]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//2*b[1]*b[2]*b[3]

    @testset "Equivalent classes probabilities" begin
      @test length(equivalent_classes(pm)) == 5
      @test issetequal(collect(keys(equivalent_classes(pm))),
                        [p[1, 1, 1], p[2, 2, 1], p[3, 2, 1], p[2, 1, 1], p[1, 2, 1]])
      @test length(equivalent_classes(pm)[p[1, 1, 1]]) == 4
      @test length(equivalent_classes(pm)[p[2, 2, 1]]) == 12
      @test length(equivalent_classes(pm)[p[3, 2, 1]]) == 24
      @test length(equivalent_classes(pm)[p[2, 1, 1]]) == 12
      @test length(equivalent_classes(pm)[p[1, 2, 1]]) == 12

      @test issetequal(equivalent_classes(pm)[p[1, 1, 1]],
                        [p[1, 1, 1], p[2, 2, 2], p[3, 3, 3], p[4, 4, 4]])
    end
  end

  @testset "parametrization GroupBasedPhylogeneticModel" begin

    f = parametrization(model)

    @test f.codomain == parameter_ring(model)[1]
    @test f.domain == model_ring(model)[1]

    @test length(f.img_gens) == n_states(model)^n_leaves(tree)

    _, q = model_ring(model)

    # Test some entries of the joint distribution vector
    x = gens(parameter_ring(model)[1])[1:3]
    y = gens(parameter_ring(model)[1])[4:6]
    @test f(p[3,2,4]) == y[1]*y[2]*y[3]
    @test f(p[2,1,2]) == x[2]*y[1]*y[3]
    @test f(p[1,1,1]) == x[1]*x[2]*x[3]


    @testset "Equivalent classes probabilities" begin
      @test length(equivalent_classes(model)) == 5
      @test issetequal(collect(keys(equivalent_classes(model))),
                       [q[1, 1, 1], q[2, 2, 1], q[2, 1, 2], p[4, 3, 2], p[1, 2, 2]])
      @test length(equivalent_classes(model)[q[1, 1, 1]]) == 1
      @test length(equivalent_classes(model)[q[2, 2, 1]]) == 3
      @test length(equivalent_classes(model)[q[2, 1, 2]]) == 3
      @test length(equivalent_classes(model)[q[4, 3, 2]]) == 6
      @test length(equivalent_classes(model)[q[1, 2, 2]]) == 3

      @test issetequal(equivalent_classes(model)[p[2, 2, 1]],
                       [p[2, 2, 1], p[3, 3, 1], p[4, 4, 1]])
    end
  end

  
  @test_skip @testset "specialized (inverse) fourier transform" begin
    # model = jukes_cantor_model(tree)

    FT = probability_ring(model).([1 1 1 1 1
    1 -1//3 -1//3 1 -1//3
    1 -1//3 1 -1//3 -1//3
    1 1 -1//3 -1//3 -1//3
    1 -1//3 -1//3 -1//3 1//3])

    IFT = probability_ring(model).([1//16 3//16 3//16 3//16 3//8
    3//16 -3//16 -3//16 9//16 -3//8
    3//16 -3//16 9//16 -3//16 -3//8
    3//16 9//16 -3//16 -3//16 -3//8
    3//8 -3//8 -3//8 -3//8 3//4])

    @test specialized_fourier_transform(model) == FT
    @test inverse_specialized_fourier_transform(model) == IFT

    p_equivclasses = compute_equivalent_classes(probability_map(model))
    f_equivclasses = compute_equivalent_classes(fourier_map(model))
    @test specialized_fourier_transform(model, p_equivclasses.classes, f_equivclasses.classes) == FT
    @test inverse_specialized_fourier_transform(model, p_equivclasses.classes, f_equivclasses.classes) == IFT
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

end
