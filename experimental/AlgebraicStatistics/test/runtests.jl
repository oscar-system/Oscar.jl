@testset "GaussianGraphicalModels" begin
  G = graph_from_edges(Directed, [[1,2],[2,3]])
  M = gaussian_graphical_model(G)
  cov_mat = covariance_matrix(M)
  
  @test vanishing_ideal(M) == ideal(
    [-cov_mat[1, 2] * cov_mat[2, 3] + cov_mat[1, 3] * cov_mat[2, 2]]
  )
end

# exported items: experimental/GraphicalModels/src/GraphicalModels.jl
# TODOs:
# * how to check e.g. collect(values(compute_equivalent_classes(probability_map(model)))) against something hardcoded? (Marina)
# * specialized (inverse) fourier transform, two functions still as comments (Christiane)

@testset "Graphical Models tests" begin
  tree = graph_from_edges(Directed,[[4,1],[4,2],[4,3]])

  @test_skip @testset "cavender_farris_neyman_model" begin
    model = cavender_farris_neyman_model(tree)         
    @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
    @test model isa GroupBasedPhylogeneticModel
    #number of states
    @test Oscar.number_states(model) == 2
    #root distribution
    @test model.phylo_model.root_distr == [1//2,1//2]
    #transition matrices
    tr_mat = Oscar.transition_matrices(model)
    @test tr_mat[Edge(4, 2)][1,1] == tr_mat[Edge(4, 2)][2,2]
    @test tr_mat[Edge(4, 2)][1,2] == tr_mat[Edge(4, 2)][2,1]
    # generators of the polynomial ring
    @test length(gens(model.phylo_model.prob_ring)) == 2(length(collect(edges(graph(model)))))
    @test length(gens(model.fourier_ring)) == 2(length(collect(edges(graph(model))))) 
    # fourier parameters
    fp = model.fourier_params
    for i in 1:3
      @test fp[Edge(4, i)][1] != fp[Edge(4, i)][2]
    end
    # group of the model
    G = group_of_model(model)
    @test is_isomorphic(unique(parent.(G))[1], abelian_group(2))

    @test vanishing_ideal(QQ, model) == vanishing_ideal(QQ, model; algorithm=:markov)
  end

  @testset "Jukes Cantor" begin
    model = jukes_cantor_model(tree)
    @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
    @test model isa PhylogeneticModel
    #number of states
    @test Oscar.number_states(model) == 4
    #root distribution
    model.root_distribution == [1//4,1//4,1//4,1//4]
    #transition matrices
    # tr_mat = Oscar.transition_matrices(model)
    # for i in 1:4, j in 1:4
    #   @test tr_mat[Edge(4, 2)][1, i == j ? 1 : 2] == tr_mat[Edge(4, 2)][i, j] #
    # end
    # generators of the polynomial ring
    @test length(gens(parameter_ring(model))) == 2(length(collect(edges(graph(model))))) 
    @test_skip @test length(gens(model_ring(model))) == 2(length(collect(edges(graph(model))))) 
     # fourier parameters
    fp = model.trans_mat_signature
    for i in 1:3
      @test fp[4, i] == fp[4, i] == fp[4, i]
    end

  end

  @test_skip @testset "kimura2_model" begin
    model = kimura2_model(tree)         
    @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
    @test model isa GroupBasedPhylogeneticModel
    # number of states
    @test Oscar.number_states(model) == 4
    # root distribution
    @test model.phylo_model.root_distr == [1//4,1//4,1//4,1//4]
    # transition matrices
    tr_mat = Oscar.transition_matrices(model)
    for i in 1:4, j in 1:4
      @test tr_mat[Edge(4, 2)][1, (i == j ? 1 : (isodd(i+j) ? 2 : 3))] == tr_mat[Edge(4, 2)][i, j]
    end
    # generators of the polynomial ring
    @test length(gens(model.phylo_model.prob_ring)) == 3(length(collect(edges(graph(model)))))
    @test length(gens(model.fourier_ring)) == 3(length(collect(edges(graph(model))))) 
    # fourier parameters
    fp = model.fourier_params
    for i in 1:3
      @test fp[Edge(4, i)][3] == fp[Edge(4, i)][4]
    end
    # group of the model
    G = group_of_model(model)
    @test is_isomorphic(unique(parent.(G))[1], abelian_group(2,2))
  end
    
  @test_skip @testset "kimura3_model" begin
    model = kimura3_model(tree)         
    @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
    # number of states
    @test Oscar.number_states(model) == 4
    # root distribution
    @test model.phylo_model.root_distr == [1//4,1//4,1//4,1//4]
    # transition matrices
    tr_mat = Oscar.transition_matrices(model)
    @test tr_mat[Edge(4, 2)][1,1] == tr_mat[Edge(4, 2)][2,2] == tr_mat[Edge(4, 2)][3,3] == tr_mat[Edge(4, 2)][4,4]
    @test tr_mat[Edge(4, 2)][2,1] == tr_mat[Edge(4, 2)][1,2] == tr_mat[Edge(4, 2)][4,3] == tr_mat[Edge(4, 2)][3,4]
    @test tr_mat[Edge(4, 2)][1,3] == tr_mat[Edge(4, 2)][3,1] == tr_mat[Edge(4, 2)][2,4] == tr_mat[Edge(4, 2)][4,2]
    @test tr_mat[Edge(4, 2)][1,4] == tr_mat[Edge(4, 2)][4,1] == tr_mat[Edge(4, 2)][2,3] == tr_mat[Edge(4, 2)][3,2]
    # generators of the polynomial ring
    @test length(gens(model.phylo_model.prob_ring)) == 4(length(collect(edges(graph(model)))))
    @test length(gens(model.fourier_ring)) == 4(length(collect(edges(graph(model))))) 
    # fourier parameters
    fp = model.fourier_params
    for i in 1:3
      @test allunique(fp[Edge(4, i)])
    end
    # group of the model
    G = group_of_model(model)
    @test is_isomorphic(unique(parent.(G))[1], abelian_group(2,2))
  end

  @test_skip @testset "general_markov_model" begin
    model = general_markov_model(tree)         
    @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
    @test model isa PhylogeneticModel
    # number of states
    @test Oscar.number_states(model) == 4
    # testing with specified states (n=20)
    model = general_markov_model(tree; number_states = 20)
    @test Oscar.number_states(model) == 20
    model = general_markov_model(tree)
    # root distribution
    @test model.root_distr[1] isa QQMPolyRingElem
    # transition matrices
    tr_mat = Oscar.transition_matrices(model)
    @test tr_mat[Edge(4,1)] isa AbstractAlgebra.Generic.MatSpaceElem{QQMPolyRingElem}
    # generators of the polynomial ring
    @test length(gens(model.prob_ring)) == 148 
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

  @test_skip @testset "Probability parametrization" begin
    p = probability_map(model)
    @test length(p) == number_states(model)^length(Oscar.leaves(tree))

    R, a, b = polynomial_ring(QQ, :a => 1:n_edges(tree), :b => 1:n_edges(tree); cached=false)
    H = Oscar.hom(probability_ring(model), R, gens(R))

    # Test some entries of the joint distribution vector
    @test H(p[1,1,1]) == 1//4*a[1]*a[2]*a[3] + 3//4*b[1]*b[2]*b[3]
    @test H(p[1,2,3]) == 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//4*b[1]*b[2]*b[3]
    @test H(p[4,4,3]) == 1//4*a[1]*a[2]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//2*b[1]*b[2]*b[3]

    @testset "Equivalent classes probabilities" begin
      p_eqclasses = compute_equivalent_classes(p)
      
      classes = [1//4*a[1]*a[2]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//2*b[1]*b[2]*b[3],
                1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//4*b[1]*b[2]*b[3],
                1//4*a[1]*a[3]*b[2] + 1//4*a[2]*b[1]*b[3] + 1//2*b[1]*b[2]*b[3],
                1//4*a[1]*b[2]*b[3] + 1//4*a[2]*a[3]*b[1] + 1//2*b[1]*b[2]*b[3],
                1//4*a[1]*a[2]*a[3] + 3//4*b[1]*b[2]*b[3]]

      # Test that we get the expected number and polynomials in equivalent classes
      @test length(unique(vcat(H.(collect(values(p_eqclasses.parametrization))), classes))) == length(p_eqclasses.parametrization)

       # Test the keys for a specific class
      @test setdiff(p_eqclasses.classes[1,1,1], [(1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4)]) == []    
    end

    @testset "Sum equivalent classes Fourier" begin
      p_eqclasses = compute_equivalent_classes(p)
      sum_p = sum_equivalent_classes(p_eqclasses)
      @test H(sum_p[1,2,1]) == 3*a[1]*a[3]*b[2] + 3*a[2]*b[1]*b[3] + 6*b[1]*b[2]*b[3]
      @test H(sum_p[1,1,1]) == a[1]*a[2]*a[3] + 3*b[1]*b[2]*b[3]
      @test H(sum_p[1,2,2]) == 3*a[1]*b[2]*b[3] + 3*a[2]*a[3]*b[1] + 6*b[1]*b[2]*b[3]
      @test H(sum_p[1,2,3]) == 6*a[1]*b[2]*b[3] + 6*a[2]*b[1]*b[3] + 6*a[3]*b[1]*b[2] + 6*b[1]*b[2]*b[3]
      @test H(sum_p[1,1,2]) == 3*a[1]*a[2]*b[3] + 3*a[3]*b[1]*b[2] + 6*b[1]*b[2]*b[3]

    end
  end

  @test_skip @testset "Fourier parametrization" begin
    q = fourier_map(model)
    @test length(q) == number_states(model)^length(Oscar.leaves(tree))
    
    S, x = polynomial_ring(QQ, :x => (1:n_edges(tree), 1:2); cached=false)
    H = Oscar.hom(fourier_ring(model), S, gens(S))

    # Test zero entries of the fourier vector
    zerokeys = [(1,2,1), (3,1,1), (4,4,2), (1,2,3), (3,2,1), (2,1,4), (3,2,3), (2,1,1), 
            (1,3,2), (1,4,2), (2,1,3), (2,2,4), (4,3,4), (4,4,4), (4,3,1), (3,3,2),
            (4,1,2), (2,2,3), (4,3,3), (4,4,3), (4,2,2), (1,3,4), (2,3,2), (1,3,1), 
            (2,4,2), (1,1,2), (1,4,1), (3,3,4), (1,4,3), (3,4,4), (4,1,1), (3,1,2), 
            (3,4,1), (3,3,3), (4,1,3), (4,2,4), (3,4,3), (4,2,1), (3,2,2), (2,3,1), 
            (2,4,4), (1,1,4), (2,4,1), (2,3,3), (1,1,3), (1,2,4), (2,2,2), (3,1,4)]
    sum([q[i] for i in zerokeys]) == 0

    # Test some entries of the joint distribution vector
    @test H(q[3,2,4]) == x[1, 2]*x[2, 2]*x[3, 2]
    @test H(q[2,1,2]) == x[2, 1]*x[1, 2]*x[3, 2]
    @test H(q[1,1,1]) == x[1, 1]*x[2, 1]*x[3, 1]

    @testset "Equivalent classes Fourier" begin
      q_eqclasses = compute_equivalent_classes(q)

      classes = [x[2, 1]*x[1, 2]*x[3, 2], x[3, 1]*x[1, 2]*x[2, 2], 
                x[1, 1]*x[2, 2]*x[3, 2], x[1, 2]*x[2, 2]*x[3, 2], 
                x[1, 1]*x[2, 1]*x[3, 1]]
  
      # Test that we get the expected number and polynomials in equivalent classes
      @test length(unique(vcat(H.(collect(values(q_eqclasses.parametrization))), classes))) == length(q_eqclasses.parametrization)

      # Test all classes
      @test H(q_eqclasses.parametrization[1,1,1]) == x[1, 1]*x[2, 1]*x[3, 1]
      @test H(q_eqclasses.parametrization[1,2,2]) == x[1, 1]*x[2, 2]*x[3, 2]
      @test H(q_eqclasses.parametrization[2,1,2]) == x[2, 1]*x[1, 2]*x[3, 2]
      @test H(q_eqclasses.parametrization[2,2,1]) == x[3, 1]*x[1, 2]*x[2, 2]
      @test H(q_eqclasses.parametrization[2,3,4]) == x[1, 2]*x[2, 2]*x[3, 2]

      # Test that the keys for one class are the expected ones
      @test setdiff(q_eqclasses.classes[2,3,4],  [(2,3,4),(2,4,3),(3,2,4),(3,4,2),(4,2,3),(4,3,2)]) == []
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
