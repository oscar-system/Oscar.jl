function evalu(x::Fac)
  return x.unit * prod(p*k for (p,k) = x.fac)
end
@testset "Polymake.factorisations" begin
  k, a = quadratic_field(-5)
  zk = maximal_order(k)
  f = factorisations(zk(6))
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
  @test_throws "disc_log failed" Oscar.DiscreteLog.disc_log(one(F), a)

  F2 = GF(3,5);
  @test_throws AssertionError disc_log(gen(F), gen(F2))
end
