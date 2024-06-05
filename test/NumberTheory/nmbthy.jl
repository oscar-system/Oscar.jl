using Oscar
using Test

function evalu(x::Fac)
  return x.unit * prod(p*k for (p,k) = x.fac)
end
@testset "Polymake.factorizations" begin
  k, a = quadratic_field(-5)
  zk = maximal_order(k)
  f = factorizations(zk(6))
  @test length(f) == 2
  @test all(x -> evalu(x) == 6, f)
end

@testset "Polymake.norm_equation" begin
  k, a = wildanger_field(3, 13)
  zk = maximal_order(k)
  na = norm(rand(zk, 1:10))
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
  Qx, x = QQ["x"]
  k, a = number_field(x^2 - 18, "a")
  kt, t = k["t"];
  K, b = number_field(t^4 + (a + 6)*t^2 + 2a + 9, "b")
  G, m = automorphism_group(PermGroup, K)
  h = m(one(G))
  @test h(b) == b && h(K(a)) == K(a)
  @test order(G) == 4 && is_cyclic(G)
end
