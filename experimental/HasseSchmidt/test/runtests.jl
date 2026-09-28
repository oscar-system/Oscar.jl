@testset "hasse_derivatives" begin
  R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);

  result_a1 = [ ([3, 0], 1),
                ([2, 0], 3*x),
                ([1, 0], 3*x^2),
                ([0, 0], x^3)]
  @test result_a1 == hasse_derivatives(x^3)

  result_a2 = [ ([0, 5], 3),
                ([0, 4], 15*y),
                ([0, 3], 30*y^2),
                ([2, 0], 5),
                ([0, 2], 30*y^3),
                ([1, 0], 10*x),
                ([0, 1], 15*y^4),
                ([0, 0], 5*x^2 + 3*y^5)]
  @test result_a2 == hasse_derivatives(5*x^2 + 3*y^5)

  result_a3 = [ ([2, 3], 1),
                ([2, 2], 3*y),
                ([1, 3], 2*x),
                ([2, 1], 3*y^2),
                ([1, 2], 6*x*y),
                ([0, 3], x^2),
                ([2, 0], y^3),
                ([1, 1], 6*x*y^2),
                ([0, 2], 3*x^2*y),
                ([1, 0], 2*x*y^3),
                ([0, 1], 3*x^2*y^2),
                ([0, 0], x^2*y^3)]
  @test result_a3 == hasse_derivatives(x^2*y^3)

  result_a4 = [ ([4, 0], 1),
                ([3, 0], 4*x),
                ([2, 0], 6*x^2),
                ([0, 2], 1),
                ([1, 0], 4*x^3),
                ([0, 1], 2*y),
                ([0, 0], x^4 + y^2)]
  @test result_a4 == hasse_derivatives(x^4 + y^2)

  result_a5 = [ ([2, 1], 1),
                ([1, 2], 1),
                ([2, 0], y),
                ([1, 1], 2*x + 2*y),
                ([0, 2], x),
                ([1, 0], 2*x*y + y^2),
                ([0, 1], x^2 + 2*x*y),
                ([0, 0], x^2*y + x*y^2)]
  @test result_a5 == hasse_derivatives(x^2*y + x*y^2)
end

@testset "hasse_derivatives finite fields" begin
  R, (x, y, z) = polynomial_ring(GF(3), ["x", "y", "z"]);

  result_b1 = [ ([2, 0, 0], 1),
                ([0, 2, 0], 1),
                ([1, 0, 0], 2*x),
                ([0, 1, 0], 2*y),
                ([0, 0, 0], x^2 + y^2)]
  @test result_b1 == hasse_derivatives(x^2 + y^2)

  result_b2 = [ ([0, 0, 6], 1),
                ([2, 1, 0], 1),
                ([0, 0, 3], 2*z^3),
                ([2, 0, 0], y),
                ([1, 1, 0], 2*x),
                ([1, 0, 0], 2*x*y),
                ([0, 1, 0], x^2),
                ([0, 0, 0], x^2*y + z^6)]
  @test result_b2 == hasse_derivatives(x^2*y + z^6)
end

@testset "_hasse_derivatives MPolyQuoRingElem" begin
  R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"]);
  I = ideal(R, [x^2 - 1]);
  RQ, _ = quo(R, I);

  result_c1 = [ ([0, 4, 0], 3),
                ([0, 3, 0], 12*y),
                ([0, 2, 0], 18*y^2),
                ([0, 1, 0], 12*y^3),
                ([0, 0, 0], 3*y^4)]
  @test result_c1 == Oscar._hasse_derivatives(RQ(3y^4))
end

@testset "_hasse_derivatives Oscar.MPolyLocRingElem" begin
  R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"]); 
  m = ideal(R, [x, y, z]); # max ideal
  U = complement_of_prime_ideal(m);
  RL, _ = localization(R, U);

  result_d1 = [ ([3, 0, 0], 5),
                ([2, 0, 0], 15*x),
                ([1, 0, 0], 15*x^2),
                ([0, 0, 0], 5*x^3)]
  @test result_d1 == Oscar._hasse_derivatives(RL(5x^3))
end

@testset "_hasse_derivatives Oscar.MPolyQuoLocRingElem" begin
  R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"]); 
  I = ideal(R, [x^2 - 1]);
  RQ, _ = quo(R, I);
  m = ideal(R, [x, y, z]); # max ideal
  U = complement_of_prime_ideal(m);
  RQL, _ = localization(RQ, U);

  result_e1 = [	([0, 0, 5], 2),
                ([0, 0, 4], 10*z),
                ([0, 0, 3], 20*z^2),
                ([0, 0, 2], 20*z^3),
                ([0, 0, 1], 10*z^4),
                ([0, 0, 0], 2*z^5)]
  @test result_e1 == Oscar._hasse_derivatives(RQL(2z^5))
end

