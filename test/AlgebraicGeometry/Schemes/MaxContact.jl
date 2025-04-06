@testset "Maximal Contact" begin
  R,(x,y,z) = polynomial_ring(QQ,3)
  J = ideal(R,[z])
  S,piS = quo(R,J)
  W = AffineScheme(S)
  I=ideal(S, [x^2-y^5])
  IS = IdealSheaf(W,I)
  WC = scheme(IS)
  inc_X = Oscar.CoveredClosedEmbedding(WC,IS)
  MC = Oscar._initialize_max_contact_object(inc_X)
  DL = Oscar._delta_list(IS)
  @test length(DL) == 2
  MC2 = Oscar._max_contact_step(MC, DL[2], 2)
  U = first(covering(MC2))
  mc_data = Oscar.maximal_contact_data(MC2)[U]
  @test Oscar.ambient_orders(mc_data) == [0,2]
  @test Oscar.dependent_variables(mc_data) == [3,1]
  @test Oscar.max_contact_minor_data(mc_data) == [(1,1),(1,1)]
  @test Oscar.prepared_jacobi_matrices(mc_data)[1] == Oscar.prepared_jacobi_matrices(maximal_contact_data(MC)[U])[1]
  R,(x,y,z) = polynomial_ring(QQ,3)
  J = ideal(R,[x-z])
  S,piS = quo(R,J)
  W = AffineScheme(S)
  I=ideal(S, [x^2+z^2-y^5])
  IS = IdealSheaf(W,I)
  WC = scheme(IS)
  inc_X = Oscar.CoveredClosedEmbedding(WC,IS)
  MC = Oscar._initialize_max_contact_object(inc_X)
  DL = Oscar._delta_list(IS)
  @test length(DL) == 2
  MC2 =_Oscar.max_contact_step(MC, DL[2], 2)
  U = first(covering(MC2))
  mc_data = Oscar.maximal_contact_data(MC2)[U]
  @test Oscar.ambient_orders(mc_data) == [0,2]
  @test Oscar.dependent_variables(mc_data) == [3,1]
  @test Oscar.max_contact_minor_data(mc_data) == [(-1,1),(-1,1)]
  @test Oscar.prepared_jacobi_matrices(mc_data)[1] == Oscar.prepared_jacobi_matrices(maximal_contact_data(MC)[U])[1]
  R,(x,y,z) = polynomial_ring(QQ,3)
  J = ideal(R,[x*y-1])
  S,piS = quo(R,J)
  W = AffineScheme(S)
  I=ideal(S, [x^2+z^2-y^5])
  IS = IdealSheaf(W,I)
  WC = scheme(IS)
  inc_X = Oscar.CoveredClosedEmbedding(WC,IS)
  MC = Oscar._initialize_max_contact_object(inc_X)
  DL = Oscar._delta_list(IS)
  @test length(DL) == 1
end