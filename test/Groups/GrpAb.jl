@testset "Groups/GrpAb.jl" begin
  A = abelian_group([4])
  f = hom(A, A, [2*gens(A)[1]])
  g = Oscar.restrict_codomain(f)
  @test issurjective(g)
  @test order(codomain(g)) == 2
  
  A = abelian_group([0, 2])
  f = hom(A, A, [gens(A)[2], zero(A)])
  g = Oscar.restrict_codomain(f)
  @test issurjective(g)
  @test order(codomain(g)) == 2
end
