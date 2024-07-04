@testset "Torus-invariant divisor classes" begin

  F5 = hirzebruch_surface(NormalToricVariety, 5)
  dP3 = del_pezzo_surface(NormalToricVariety, 3)
  P2 = projective_space(NormalToricVariety, 2)
    
  DC = toric_divisor_class(F5, [ZZRingElem(0), ZZRingElem(0)])
  DC2 = toric_divisor_class(F5, [1, 2])
  DC3 = toric_divisor_class(dP3, [4, 3, 2, 1])
  DC4 = canonical_divisor_class(dP3)
  DC5 = anticanonical_divisor_class(dP3)
  DC6 = trivial_divisor_class(dP3)
  DC7 = toric_divisor_class(P2, [1])
  DC8 = toric_divisor_class(P2, [-1])
  
  @testset "Basic properties" begin
    @test is_trivial(toric_divisor(DC2)) == false
    @test is_effective(DC7) == true
    @test is_effective(DC8) == false
  end
  
  @testset "Basic attributes" begin
    @test torsion_free_rank(parent(divisor_class(DC2))) == 2
    @test dim(toric_variety(DC2)) == 2
  end
  
  @testset "Arithmetic" begin
    @test is_trivial(ZZRingElem(2)*DC+DC2) == false
    @test is_trivial(2*DC-DC2) == false
    @test (DC == DC2) == false
    @test (DC4 - DC5 == DC6) == false
    @test (DC == DC3) == false
  end
end
