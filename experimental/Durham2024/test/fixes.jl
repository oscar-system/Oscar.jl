@testset "components of cycles" begin
  IA2 = affine_space(QQ, [:x, :y])
  bl0 = blow_up(IA2, ideal(OO(IA2), gens(OO(IA2))))
  E = exceptional_divisor(bl0)
  EW = weil_divisor(E)
  c = first(components(WeilDivisor, EW))
  D = Oscar.underlying_cycle(EW)
  cc = first(components(AlgebraicCycle, D))
  @test EW[c] == D[cc]
end

