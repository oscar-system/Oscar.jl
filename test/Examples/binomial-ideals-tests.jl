Oscar.include(joinpath(Oscar.oscardir, "examples", "TranslateBinomial2.jl"))

@testset "Binomial Ideals" begin
  Qxy, (x, y, z, t) = PolynomialRing(FlintQQ, 4)
  I = ideal([x*y, z*t^2-t^3, z^2-y^2])
  @test Oscar.isbinomial(I)
  @test Oscar.isunital(I)
  J = ideal([x*y - z*t^2 + t^3, z*t^2-t^3])
  @test Oscar.isbinomial(J)
  J1 = ideal(Qxy, x^2+y^2+z^2)
  @test !Oscar.isbinomial(J1)
  @test !Oscar.iscellular(I)[1]
  lI = Oscar.cellular_decomposition(I)
  lI2 = Oscar.cellular_decomposition_macaulay(I)
  @test length(lI) == length(lI2)
  for x in lI
    @test x in lI2
  end
  lP = Oscar.primary_decomposition(I)
  lP1 = Oscar.binomial_primary_decomposition(I)
  @test length(lP) == length(lP1)
  for x in lP
    @test x in lP1
  end

end
