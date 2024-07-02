@testset "hasse_derivatives" begin
  R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);
  @test [x^3, 3x^2, 3x, R(1)] == hasse_derivatives(R(x^3))
  @test [5x^2 + 3y^5, 15y^4, 30y^3, 30y^2, 15y, R(3), 10x, R(5)] == hasse_derivatives(R(5*x^2 + 3*y^5))
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

end


R, x = polynomial_ring(QQ, 4, "x");
m = ideal(R, [x, y, z]);
I = ideal(R, [x[1]^3 - 1]);
RQ, _ = quo(R, I);
p = ideal(R, [x[2]]);
U = complement_of_prime_ideal(p);
RQL, _ = localization(RQ, U);
f1 = RQ(4*x[3]^3);
f2 = RQL(3*x[2]^2);
_hasse_derivatives([f1, f2])

f1 = RL(5x^3);
f2 = RQ(3y^4);
f3 = RQL(2z^5);
_hasse_derivatives([f1, f2, f3])