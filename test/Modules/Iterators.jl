@testset "iterators over monomials in modules" begin
  IP1 = projective_space(NormalToricVariety, 1)
  X = IP1*IP1
  S = cox_ring(X)
  G = grading_group(S)
  F = graded_free_module(S, [zero(G), G[1], 3*G[2]])
  @test length(collect(Oscar.all_monomials(F, 7*G[2]+3*G[1]))) == 76
  @test length(collect(Oscar.all_exponents(F, 7*G[2]+3*G[1]))) == 76
end
