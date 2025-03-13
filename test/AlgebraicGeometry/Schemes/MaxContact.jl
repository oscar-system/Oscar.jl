@testset "Maximal Contact" begin
  R,(x,y,z) = polynomial_ring(QQ,3)
  J = ideal(R,[z])
  S,piS = quo(R,J)
  W = AffineScheme(S)
  I=ideal(S, [x^2-y^5])
  IS = IdealSheaf(W,I)
  WC = scheme(IS)
  inc_X = Oscar.CoveredClosedEmbedding(WC,IS)
  MC = _initialize_max_contact_object(inc_X)
  DL = Oscar._delta_list(IS)
  @test length(DL) == 2
  MC2 =_max_contact_step(MC, DL[2], 2)
  U = first(covering(MC2))
  mc_data = maximal_contact_data(MC2)[U]
  @test mc_data[2] == [0,2]
  amb_data = ambient_parameter_data(MC2)[U]
  @test amb_data[1] == [3,1]
  @test amb_data[2] == [(1,1),(1,1)]
  amb_data0 = ambient_parameter_data(MC)[U]
  @test amb_data[3][1] == amb_data0[3][1]
  R,(x,y,z) = polynomial_ring(QQ,3)
  J = ideal(R,[z-x])
  S,piS = quo(R,J)
  W = AffineScheme(S)
  I=ideal(S, [x^2+z^2-y^5])
  IS = IdealSheaf(W,I)
  WC = scheme(IS)
  inc_X = Oscar.CoveredClosedEmbedding(WC,IS)
  MC = _initialize_max_contact_object(inc_X)
  DL = Oscar._delta_list(IS)
  @test length(DL) == 2
  MC2 =_max_contact_step(MC, DL[2], 2)
  U = first(covering(MC2))
  mc_data = maximal_contact_data(MC2)[U]
  @test mc_data[2] == [0,2]
  amb_data = ambient_parameter_data(MC2)[U]
  @test amb_data[1] == [3,1]
  @test amb_data[2] == [(1,1),(1,1)]
  amb_data0 = ambient_parameter_data(MC)[U]
  @test amb_data[3][1] == amb_data0[3][1]
end