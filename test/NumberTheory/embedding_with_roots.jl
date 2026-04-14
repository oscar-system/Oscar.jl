@testset "root lattice embeddings" begin 
  #the Nikulin lattice is an overlattice of A_1^8
  N=integer_lattice(gram=matrix([-4 1 1 1 1 1 1 1;1 -2 0 0 0 0 0 0; 1 0 -2 0 0 0 0 0; 1 0 0 -2 0 0 0 0; 1 0 0 0 -2 0 0 0; 1 0 0 0 0 -2 0 0; 1 0 0 0 0 0 -2 0; 1 0 0 0 0 0 0 -2 ]))
  A=matrix([-2 0 0 0 0 0 0 0; 0 -2 0 0 0 0 0 0; 0 0 -2 0 0 0 0 0; 0 0 0 -2 0 0 0 0; 0 0 0 0 -2 0 0 0; 0 0 0 0 0 -2 0 0; 0 0 0 0 0 0 -2 0; 0 0 0 0 0 0 0 -2])
  @test gram_matrix(biggest_root_sublattice(genus(N)))==A
  
  NS=direct_sum(N, integer_lattice(gram=matrix([0 1;1 0])))[1]
  #this NS is U+N; the genus g of the orthogonal complement  of U+N in Lambda26 has many representatives 
  #g is the genus of the Kummer lattice, so the maximum number of roots is 16
  Or=embedding_in_unimodular_manyroots(NS,1,25,primitive=true,even=true, compute_overlattices=false)
  OR=Or[4]
  @test rank(root_sublattice(OR))==16
  
  #this uses the alternative option of compute_overlattices
  NS2=integer_lattice(gram=matrix([2]))
  Or2=embedding_in_unimodular_manyroots(NS2,1,1,primitive=true,even=true, compute_overlattices=true)
  OR2=Or2[4]
  @test rank(root_sublattice(OR2))==1
end

