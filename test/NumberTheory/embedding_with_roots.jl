@testset "root lattice embeddings" begin 
  #the Nikulin lattice is an overlattice of A_1^8
  N=integer_lattice(gram=matrix([-4 1 1 1 1 1 1 1;1 -2 0 0 0 0 0 0; 1 0 -2 0 0 0 0 0; 1 0 0 -2 0 0 0 0; 1 0 0 0 -2 0 0 0; 1 0 0 0 0 -2 0 0; 1 0 0 0 0 0 -2 0; 1 0 0 0 0 0 0 -2 ]))
  A=matrix([-2 0 0 0 0 0 0 0; 0 -2 0 0 0 0 0 0; 0 0 -2 0 0 0 0 0; 0 0 0 -2 0 0 0 0; 0 0 0 0 -2 0 0 0; 0 0 0 0 0 -2 0 0; 0 0 0 0 0 0 -2 0; 0 0 0 0 0 0 0 -2])
  @test gram_matrix(biggest_root_sublattice(genus(N)))==A 
end