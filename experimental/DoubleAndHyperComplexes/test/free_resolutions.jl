@testset "free resolutions of ideals" begin
  X = veronese();
  I = defining_ideal(X);
  Pn = base_ring(I)
  FI = free_resolution(Oscar.SimpleFreeResolution, I);
  F = graded_free_module(Pn, 1)
  dualFIC = hom(FI, F);
  # Duals invert the range from 0:10 to 0:-1:-10 etc. 
  @test dualFIC[0] isa FreeMod
  @test dualFIC[-1] isa FreeMod
  @test dualFIC[-10] isa FreeMod
  @test is_graded(dualFIC[0])
  @test is_graded(dualFIC[-1])
  @test is_graded(dualFIC[-10])
  @test !Oscar.can_compute_index(dualFIC, 1)
end
