@testset "Modules over MPolyQuos" begin
  P, (x, y) = QQ["x", "y"]
  I = ideal(P, x)
  R, _ = quo(P, I)
  FR = FreeMod(R, 2, cached=false)
  M = oscar._as_poly_module(FR)
  @test M isa SubQuo{<:MPolyElem}
  @test FR isa FreeMod{<:MPolyQuoElem}
  FP = oscar._poly_module(FR)
  @test FP isa FreeMod{<:MPolyElem}

  M = FreeMod(R, 1, cached=false)
  psi = hom(FR, M, [2*M[1], M[1]])

  K, inc = kernel(psi)

  v = FR[1] - 2*FR[2] 
  @test v in K

  N = oscar.SubModuleOfFreeModule(FR, [R(y-1)*FR[1], R(y)*FR[2]])
  c = coordinates(y*(y-1-x)*FR[1] + y^5*(1+x)*FR[2], N)
  @test c[1] == R(y)
  @test c[2] == R(y^4)
end
