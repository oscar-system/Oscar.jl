@testset "Groups/GrpAb.jl" begin
  A = abelian_group([4])
  f = hom(A, A, [2*gens(A)[1]])
  g = Oscar.restrict_codomain(f)
  @test is_surjective(g)
  @test order(codomain(g)) == 2
  
  A = abelian_group([0, 2])
  f = hom(A, A, [gens(A)[2], zero(A)])
  g = Oscar.restrict_codomain(f)
  @test is_surjective(g)
  @test order(codomain(g)) == 2
end

@testset "describe for GrpAbFinGen" begin
  @test describe(abelian_group(GrpAbFinGen, Int[])) == "0"
  @test describe(abelian_group(GrpAbFinGen, Int[0])) == "Z"
  @test describe(abelian_group(GrpAbFinGen, Int[0, 0])) == "Z^2"
  @test describe(abelian_group(GrpAbFinGen, Int[2])) == "Z/2"
  @test describe(abelian_group(GrpAbFinGen, Int[0, 2])) == "Z/2 + Z"
  @test describe(abelian_group(GrpAbFinGen, Int[0, 0, 2])) == "Z/2 + Z^2"
  @test describe(abelian_group(GrpAbFinGen, Int[2, 4])) == "Z/2 + Z/4"
  @test describe(abelian_group(GrpAbFinGen, Int[0, 2, 4])) == "Z/2 + Z/4 + Z"
  @test describe(abelian_group(GrpAbFinGen, Int[0, 0, 2, 4])) == "Z/2 + Z/4 + Z^2"
end
