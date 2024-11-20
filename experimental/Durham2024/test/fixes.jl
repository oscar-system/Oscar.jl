@testset "components of cycles" begin
  IA2 = affine_space(QQ, [:x, :y])
  bl0 = blow_up(IA2, ideal(OO(IA2), gens(OO(IA2))))
  E = exceptional_divisor(bl0)
  EW = weil_divisor(E)
  c = first(components(WeilDivisor, EW))
  D = Oscar.underlying_cycle(EW)
  cc = first(components(AlgebraicCycle, D))
  @test EW[c] == D[cc]
  
  X = domain(bl0)
  FX = function_field(X)
  U = first(affine_charts(X))
  x = first(gens(OO(U)))
  f = FX(x)
  D = weil_divisor(f)
  @test 2*D == weil_divisor(f^2) 
end

