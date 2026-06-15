@testset "root_overlattices" begin
 @test length(Oscar.root_overlattices(4))==6
end

@testset "invariant functions for integer quadratic lattices" begin
  A = [[2 -1 0 0 0 0; -1 2 -1 0 0 0; 0 -1 2 -1 0 0; 0 0 -1 2 -1 0; 0 0 0 -1 2 0; 0 0 0 0 0 20], [2 0 0 0 -1 -1; 0 2 0 -1 0 -1; 0 0 2 -1 1 0; 0 -1 -1 4 1 2; -1 0 1 1 4 1; -1 -1 0 2 1 4], [2 -1 1 0 0 0; -1 2 -1 0 0 0; 1 -1 2 0 0 0; 0 0 0 2 0 0; 0 0 0 0 2 1; 0 0 0 0 1 8], [2 1 -1 -1 0 0; 1 2 -1 -1 0 0; -1 -1 2 1 0 0; -1 -1 1 2 0 0; 0 0 0 0 2 0; 0 0 0 0 0 12], [2 -1 0 0 0 -1; -1 2 0 0 0 0; 0 0 2 0 1 0; 0 0 0 2 1 0; 0 0 1 1 4 0; -1 0 0 0 0 4], [2 -1 1 0 -1 -1; -1 2 -1 0 1 1; 1 -1 2 0 0 0; 0 0 0 2 0 0; -1 1 0 0 4 1; -1 1 0 0 1 6], [2 -1 1 1 -1 0; -1 2 -1 -1 0 0; 1 -1 2 0 0 0; 1 -1 0 2 -1 0; -1 0 0 -1 2 0; 0 0 0 0 0 30]]
  L = [integer_lattice(gram=matrix(QQ,6,6,i),cached=false) for i in A]
  @test 7 == length(unique!(Oscar.invariant_function_graph_hash.(L)))
  @test 7 == length(unique!(Oscar.oscar_invariant_function.(L)))
end 
@testset "Lattice Canonical Form" begin
  function create_gram(data)  # helper function to create gram matricies of latticies
    U = upper_triangular_matrix(data)
    return U + transpose(U) - diagonal_matrix(diagonal(U))
  end

  n = 10
  matrix1 = [2,1,-1,1,1,0,0,0,0,-1,2,-1,0,0,0,0,0,0,-1,2,0,0,0,0,0,0,0,2,1,0,0,0,0,-1,2,0,0,0,0,-1,4,-2,2,2,0,4,-2,-2,0,4,0,0,4,0,4]
  matrix2 = [2,1,1,-1,-1,1,-1,1,-1,-1,2,0,-1,-1,1,-1,1,-1,0,4,-2,1,2,-2,-1,-1,0,4,1,-2,2,1,2,1,4,1,-1,0,0,1,4,-2,1,-2,-1,4,-1,2,1,4,0,-1,4,1,4]
  G1 = create_gram(matrix1)
  G2 = create_gram(matrix2)
  U = hnf_with_transform(matrix(ZZ,n,n,rand(0:1,n^2)))[2];
  L1 = integer_lattice(gram = G1);
  L2 = integer_lattice(gram = G2);
  L3 = lattice_in_same_ambient_space(L1,U*basis_matrix(L1));
  @test canonical_form(L1) != canonical_form(L2)
  @test canonical_form(L1) == canonical_form(L3)

  n = 17
  matrix1 = [2,1,-1,0,0,0,1,1,0,1,1,0,-1,0,0,-1,-1,2,0,0,0,0,1,1,0,0,1,0,0,0,0,-1,-1,2,0,0,0,-1,-1,0,0,-1,0,1,0,0,1,0,2,0,0,1,1,0,-1,1,0,0,0,-1,0,0,2,0,1,-1,-1,1,-1,0,1,0,0,-1,1,2,1,1,0,-1,-1,1,1,-1,0,-1,1,4,2,-1,0,0,1,1,0,-1,-2,0,4,1,-1,1,1,0,0,-1,-1,-1,2,-1,1,0,0,0,0,1,0,4,-1,0,-1,1,0,0,-1,4,-1,-1,0,0,0,-1,2,0,0,0,0,0,4,-1,0,-1,1,2,0,1,-1,2,0,0,4,0,4]
  matrix2 = [2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,2,1,1,1,1,1,-1,1,0,0,1,0,1,1,0,2,1,1,1,1,-1,1,-1,1,0,1,0,0,-1,2,1,1,1,-1,1,-1,1,1,1,1,1,-1,2,1,1,-1,1,0,0,0,0,0,0,-1,2,1,-1,1,-1,1,1,0,1,1,-1,2,-1,1,0,0,0,1,0,0,-1,2,-1,1,-1,0,-1,0,0,1,2,-1,1,0,1,0,0,-1,4,-3,-1,-1,-1,-1,2,4,1,1,1,1,-2,4,0,3,3,0,4,-1,0,-1,4,3,0,4,0,4]
  G1 = create_gram(matrix1)
  G2 = create_gram(matrix2)
  U = hnf_with_transform(matrix(ZZ,n,n,rand(0:1,n^2)))[2];
  L1 = integer_lattice(gram = G1);
  L2 = integer_lattice(gram = G2);
  L3 = lattice_in_same_ambient_space(L1,U*basis_matrix(L1));
  #@test canonical_form(L1) != canonical_form(L2)
  #test canonical_form(L1) == canonical_form(L3)

  n = 18
  matrix1 = [4,0,0,2,1,2,2,1,0,0,0,0,0,0,0,0,0,0,4,0,2,3,1,0,2,0,0,0,0,0,0,0,0,0,0,4,2,3,3,1,1,0,0,0,0,0,0,0,0,0,0,4,4,3,2,2,0,0,0,0,0,0,0,0,0,0,6,4,2,3,0,0,0,0,0,0,0,0,0,0,4,2,2,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,4,0,0,2,1,2,2,1,0,0,4,0,2,3,1,0,2,0,0,4,2,3,3,1,1,0,0,4,4,3,2,2,0,0,6,4,2,3,0,0,4,2,2,0,0,2,1,0,0,2,0,0,4,0,8]
  matrix2 = [2,1,1,1,0,0,0,0,1,1,-1,1,-1,0,-1,1,1,1,2,1,0,0,0,0,0,0,0,-1,1,-1,0,-1,1,0,1,2,0,0,0,0,0,1,1,-1,0,0,0,0,0,1,0,2,0,0,0,0,1,1,0,0,0,0,-1,0,0,1,2,0,0,0,-1,-1,-1,-1,1,0,0,0,0,0,2,1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,1,2,-1,-1,-1,-1,1,-1,0,-1,-1,0,1,2,1,1,0,-1,1,0,0,0,1,-1,4,3,0,-1,1,0,1,1,2,0,4,1,-1,1,0,0,0,1,-1,4,0,1,-1,1,0,-1,-2,4,-3,0,-1,1,0,1,4,0,1,-1,0,-2,2,-1,-1,-1,0,4,1,1,-1,4,1,1,4,0,4]
  G1 = create_gram(matrix1)
  G2 = create_gram(matrix2)
  U = hnf_with_transform(matrix(ZZ,n,n,rand(0:1,n^2)))[2];
  L1 = integer_lattice(gram = G1);
  L2 = integer_lattice(gram = G2);
  L3 = lattice_in_same_ambient_space(L1,U*basis_matrix(L1));
  @info "first try of long calcs"
  @test canonical_form(L1) != canonical_form(L2)
  @info "second try of long calcs"
  @test canonical_form(L1) != canonical_form(L2)
  @test canonical_form(L1) == canonical_form(L3)
  @info "is_isometric_with_isometry 1"
  @time is_isometric_with_isometry(L1, L2)
  @info "is_isometric_with_isometry 2"
  @time is_isometric_with_isometry(L1, L2)
end