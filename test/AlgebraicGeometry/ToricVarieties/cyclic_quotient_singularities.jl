@testset "Cyclic Quotient Singularities" begin
  cyc = cyclic_quotient_singularity(2, 1)

  @testset "Basic properties" begin
    @test is_affine(cyc) == true
    @test is_smooth(cyc) == false
    @test is_simplicial(cyc) == true
    @test is_projective_space(cyc) == false
  end

  @testset "Basic attributes" begin
    @test continued_fraction_hirzebruch_jung(cyc)[1] == 2
    @test dual_continued_fraction_hirzebruch_jung(cyc)[1] == 2
  end
end
