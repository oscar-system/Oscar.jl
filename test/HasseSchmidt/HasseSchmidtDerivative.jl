@testset "hasse_derivatives" begin
  R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);
  @test [x^3, 3x^2, 3x, R(1)] == hasse_derivatives(R(x^3))
  @test [5x^2 + 3y^5, 15y^4, 30y^3, 30y^2, 15y, R(3), 10x, R(5)] == hasse_derivatives(R(5*x^2 + 3*y^5))
  @test [x^2*y^3, 2*x*y^3, y^3, 3*x^2*y^2, 6*x*y^2, 3*y^2, 3*x^2*y, 6*x*y, 3*y, x^2, 2*x, R(1)] == hasse_derivatives(R(x^2*y^3))
  @test [y^2 + z^4, 4*z^3, 6*z^2, 4*z, R(1), 2*y, R(1)] == hasse_derivatives(R(z^4 + y^2))
end

@testset "hasse_derivatives list" begin
  R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);
  @test [[3*x^2, 6*x, R(3)], [6*y^3, 18*y^2, 18*y, R(6)]] == hasse_derivatives([R(3*x^2), R(6*y^3)])
end

@testset "_hasse_derivatives MPolyQuoRingElem" begin
  R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"]);
  I = ideal(R, [x^2 - 1]);
  RQ, _ = quo(R, I);
  @test [3*y^4, 12*y^3, 18*y^2, 12*y, RQ(3)] == _hasse_derivatives(RQ(3y^4))
end

@testset "_hasse_derivatives Oscar.MPolyLocRingElem" begin
  R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"]); 
  m = ideal(R, [x, y, z]); # max ideal
  U = complement_of_prime_ideal(m);
  RL, _ = localization(R, U);
  f1 = RL(5x^3);
  @test [5*x^3, 15*x^2, 15*x, RL(5)] == _hasse_derivatives(RL(5x^3))
end

@testset "_hasse_derivatives Oscar.MPolyQuoLocRingElem" begin
  R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"]); 
  m = ideal(R, [x, y, z]); # max ideal
  U = complement_of_prime_ideal(m);
  RL, _ = localization(R, U); # loc at m
  I = ideal(R, [x^2 - 1]);
  RQ, _ = quo(R, I);
  RQL, _ = localization(RQ, U);
  @test [2*z^5, 10*z^4, 20*z^3, 20*z^2, 10*z, RQL(2)] == _hasse_derivatives(RQL(2z^5))
end

@testset "_hasse_derivatives list" begin
  R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"]); 
  m = ideal(R, [x, y, z]); # max ideal
  U = complement_of_prime_ideal(m);
  RL, _ = localization(R, U); # loc at m
  I = ideal(R, [x^2 - 1]);
  RQ, _ = quo(R, I);
  RQL, _ = localization(RQ, U);
  @test [[5*x^3, 15*x^2, 15*x, 5], [3*y^4, 12*y^3, 18*y^2, 12*y, 3], [2*z^5, 10*z^4, 20*z^3, 20*z^2, 10*z, 2]] == _hasse_derivatives([RL(5x^3), RQ(3y^4), RQL(2z^5)])
end