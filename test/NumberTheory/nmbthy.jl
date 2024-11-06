using Oscar
using Test

function evalu(x::Fac)
  return x.unit * prod(p*k for (p,k) = x.fac)
end

@testset "factorizations" begin
  k, a = quadratic_field(-5)
  zk = maximal_order(k)
  f = factorizations(zk(6))
  @test length(f) == 2
  @test all(x -> evalu(x) == 6, f)
end

@testset "norm_equation.absolute" begin
  k, a = wildanger_field(3, 13)
  zk = maximal_order(k)
  na = norm(rand(zk, 1:10))
  s = norm_equation(zk, na)
  @test all(x->norm(x) == na, s)
end

@testset "norm_equation.relative" begin
  L, _ = quadratic_field(-1)
  _, x = L[:x]
  f = x^12 - 34734*x^11 + 401000259*x^10 - 1456627492885*x^9 - 2537142937228035*x^8 + 18762072755679375516*x^7 - 812368636358864062944*x^6 - 70132863629758257512231931*x^5 + 25834472514893102332821062085*x^4 + 76623280610352450247247939584745*x^3 - 45080885015422662132515763499758450*x^2 - 2070499552240812214288316981071818900*x - 550505759097778545485364826246753544;
  k, _ = number_field(f)

  zk = maximal_order(k)
  na = norm(rand(zk, 10))
  s = norm_equation(zk, na)
  @test all(x->norm(x) == na, s)
end

@testset "DiscreteLog" begin

  F = GF(3,4);
  a = gen(F)^21;
  @test disc_log(gen(F), a) == 21
  @test_throws "disc_log failed" disc_log(one(F), a)

  F2 = GF(3,5);
  @test_throws AssertionError disc_log(gen(F), gen(F2))
end

begin
  Qx, x = QQ[:x]
  k, a = number_field(x^2 - 18, "a")
  kt, t = k[:t];
  K, b = number_field(t^4 + (a + 6)*t^2 + 2a + 9, "b")
  G, m = automorphism_group(PermGroup, K)
  h = m(one(G))
  @test h(b) == b && h(K(a)) == K(a)
  @test order(G) == 4 && is_cyclic(G)
end

@testset "norm equation non-max" begin
  k, a = wildanger_field(3, 13)
  zk = maximal_order(k)
  b = zk(8*a^2 - 24*a - 1)
  o = Order(k, 8 .* basis(zk))
  @test length(norm_equation(zk, norm(b))) == 4
  @test length(norm_equation(o, norm(b))) == 3

  @test length(norm_equation(o, 8)) == 2
  @test length(norm_equation(zk, 8)) == 1
end
