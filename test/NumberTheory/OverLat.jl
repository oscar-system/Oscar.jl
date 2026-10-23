@testset "overlattices with no new roots" begin 
  # this lattice has 6 valid overlattices (including itself) belonging to different genera. The discriminant of L is 2^4*3^4, so we're testing different primes and powers of them.
  L=direct_sum(root_lattice(:A,5),root_lattice(:A,5),root_lattice(:A,5),root_lattice(:A,5))[1];
  ovs=Oscar.overlattices_no_new_roots(L);
  @test length(unique([genus(F) for F in ovs])) == 6

  #this lattice has 3 non-isometric overlattices of index 2, two of them belonging to the same genus
  L2=integer_lattice(gram=diagonal_matrix([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
  ov2=Oscar._overlattices_no_new_roots_index_p(L2,ZZ(2))
  @test length(ov2)==3
  @test length(unique([genus(F) for F in ov2])) == 2
end
