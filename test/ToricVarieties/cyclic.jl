@testset "Cyclic Quotient Singularities" begin
    cyc = CyclicQuotientSingularity(2, 1)

    @testset "Basic properties" begin
        @test is_affine(cyc) == true
        @test is_smooth( cyc ) == false
        @test is_simplicial( cyc ) == true
        @test continued_fraction_hirzebruch_jung( cyc )[1] == 2
        @test dual_continued_fraction_hirzebruch_jung(cyc)[1] == 2
        @test is_projective_space(cyc) == false
    end
end

