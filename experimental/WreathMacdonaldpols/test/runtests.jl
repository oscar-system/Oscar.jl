@testset "wreath macdonald polynomials" begin

  result_p1 = [t^2 t 1; q t 1; q q^2 1]
  @test result_p1 == wreath_macs(1,3,@perm r (3,1,2), [1,-1,0])
end


