# exported items: experimental/GraphicalModels/src/GraphicalModels.jl
# TODOs:
# * how to check e.g. collect(values(compute_equivalent_classes(probability_map(model)))) against something hardcoded? (Marina)
# * Markov model: root distribution, transition matricies, question: why 148 gens of ring?
# * @testset affine 
# * specialized (inverse) fourier transform,two fuctions still as comments (Christiane)

@testset "Graphical Models tests" begin

  tree = graph_from_edges(Directed,[[4,1],[4,2],[4,3]])

  @testset "cavender_farris_neyman_model" begin
    model = cavender_farris_neyman_model(tree)         
    @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
    @test model isa GroupBasedPhylogeneticModel
    #number of states
    @test Oscar.number_states(model) == 2
    #root distribution
    @test model.phylo_model.root_distr == [1//2,1//2]
    #transition matricies
    tr_mat = Oscar.transition_matrices(model)
    @test tr_mat[Edge(4, 2)][1,1] == tr_mat[Edge(4, 2)][2,2]
    @test tr_mat[Edge(4, 2)][1,2] == tr_mat[Edge(4, 2)][2,1]
    # generators of the polynomial ring
    @test length(gens(model.phylo_model.prob_ring)) == 2(length(collect(edges(graph(model)))))
    @test length(gens(model.fourier_ring)) == 2(length(collect(edges(graph(model))))) 
    # fourier parametrisation
    fp = model.fourier_params
    for i in 1:3
      @test fp[Edge(4, i)][1] != fp[Edge(4, i)][2]
    end
    # group of the model
    @test model.group == [[0],[1]]
  end

  @testset "Jukes Cantor" begin
    model = jukes_cantor_model(tree)
    @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
    @test model isa GroupBasedPhylogeneticModel
    #number of states
    @test Oscar.number_states(model) == 4
    #root distribution
    model.phylo_model.root_distr == [1//4,1//4,1//4,1//4]
    #transition matricies
    tr_mat = Oscar.transition_matrices(model)
    for i in 1:4, j in 1:4
      @test tr_mat[Edge(4, 2)][1, i == j ? 1 : 2] == tr_mat[Edge(4, 2)][i, j] #
    end
    # generators of the polynomial ring
    @test length(gens(model.phylo_model.prob_ring)) == 2(length(collect(edges(graph(model))))) 
    @test length(gens(model.fourier_ring)) == 2(length(collect(edges(graph(model))))) 
    # fourier parametrisation
    fp = model.fourier_params
    for i in 1:3
      @test fp[Edge(4, i)][2] == fp[Edge(4, i)][3] == fp[Edge(4, i)][4]
    end
    #group of the model
    @test model.group == [[0,0],[0,1],[1,0],[1,1]]
  end

  @testset "kimura2_model" begin
    model = kimura2_model(tree)         
    @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
    @test model isa GroupBasedPhylogeneticModel
    # number of states
    @test Oscar.number_states(model) == 4
    # root distribution
    @test model.phylo_model.root_distr == [1//4,1//4,1//4,1//4]
    # transition matricies
    tr_mat = Oscar.transition_matrices(model)
    for i in 1:4, j in 1:4
      @test tr_mat[Edge(4, 2)][1, (i == j ? 1 : (isodd(i+j) ? 2 : 3))] == tr_mat[Edge(4, 2)][i, j]
    end
    # generators of the polynomial ring
    @test length(gens(model.phylo_model.prob_ring)) == 3(length(collect(edges(graph(model)))))
    @test length(gens(model.fourier_ring)) == 3(length(collect(edges(graph(model))))) 
    # fourier parametrisation
    fp = model.fourier_params
    for i in 1:3
      @test fp[Edge(4, i)][3] == fp[Edge(4, i)][4]
    end
    # group of the model
    @test model.group == [[0,0],[0,1],[1,0],[1,1]]
  end
    
  @testset "kimura3_model" begin
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
    # fourier parametrisation
    fp = model.fourier_params
    #HOW TO COMPARE SECOND INDEX OF VARIABLE?
    # group of the model
    @test model.group == [[0,0],[0,1],[1,0],[1,1]]
  end

  @testset "general_markov_model" begin
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
    # transition matrices
    # generators of the polynomial ring
    @test length(gens(model.prob_ring)) == 148 
  end

  @testset "affine models" begin 
  model = jukes_cantor_model(tree)
  @test Oscar.affine_phylogenetic_model(model) isa GroupBasedPhylogeneticModel
  #sum of each row of transition matrices should equal one 
  tr_mat = Oscar.transition_matrices(Oscar.affine_phylogenetic_model(model))
  for j in 1:4
    rowsum = 0
    for i in 1:4
        rowsum = rowsum+tr_mat[Edge(4, 2)][j,i]
    end
    @test rowsum == 1
  end
  #sum of root distribution should equal one
  @test sum(Oscar.affine_phylogenetic_model(model).phylo_model.root_distr) ==1
  end

  #Now we test functionalities for one tree and one model
  model = jukes_cantor_model(tree)

  @testset "probability map" begin 
    @test length(probability_map(model))== 64
    # p[1,2,3,]==... check about 2 random examples
  end

  @testset "fourrier map" begin
    @test length(fourier_map(model)) == 64
    # q[1,2,3,]==... check about 2 random examples
  end

  @testset "compute_equivalence_classes fuctions" begin
    #check that number of equivalence classes is correct
    @test length(compute_equivalent_classes(probability_map(model))) == 5
    #collect(values(compute_equivalent_classes(probability_map(model))))
    #take probabliyring(pm), define a polynomial in there: How? Marina writes something.  
    #Until than like for website, new ring and ringhomomorphism.
    #compare like with website (see code) unique! for testing if the same
    #check that number of equivalence classes in fourrier coordinates is correct
    @test length(compute_equivalent_classes(fourier_map(model))) == 6
    #same here
  end

  @testset "specialized (inverse) fourier transform" begin
    model = jukes_cantor_model(tree)
    p_equivclasses = sum_equivalent_classes(probability_map(model))
    f_equivclasses = sum_equivalent_classes(fourier_map(model))
    @test specialized_fourier_transform(model) == model.phylo_model.prob_ring.([1 1 1 1 1
    1 -1//3 -1//3 1 -1//3
    1 -1//3 1 -1//3 -1//3
    1 1 -1//3 -1//3 -1//3
    1 -1//3 -1//3 -1//3 1//3])
    #=@test specialized_fourier_transform(model,p_equivclasses,f_equivclasses) == model.phylo_model.prob_ring.([1 1 1 1 1
    1 -1//3 -1//3 1 -1//3
    1 -1//3 1 -1//3 -1//3
    1 1 -1//3 -1//3 -1//3
    1 -1//3 -1//3 -1//3 1//3])=#
    @test inverse_specialized_fourier_transform(model) == model.phylo_model.prob_ring.([1//16 3//16 3//16 3//16 3//8
    3//16 -3//16 -3//16 9//16 -3//8
    3//16 -3//16 9//16 -3//16 -3//8
    3//16 9//16 -3//16 -3//16 -3//8
    3//8 -3//8 -3//8 -3//8 3//4])
    #=@test inverse_specialized_fourier_transform(model,p_equivclasses,f_equivclasses) == model.phylo_model.prob_ring.([1//16 3//16 3//16 3//16 3//8
    3//16 -3//16 -3//16 9//16 -3//8
    3//16 -3//16 9//16 -3//16 -3//8
    3//16 9//16 -3//16 -3//16 -3//8
    3//8 -3//8 -3//8 -3//8 3//4])=#
    end
end
