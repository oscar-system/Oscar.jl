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
  
  (s, x, y) = ambient_coordinates(U)
  I = ideal_sheaf(E)
  IU = I(U)
  @test OO(U)(x) in IU
  @test OO(U)(y) in IU
  @test !(OO(U)(s) in IU)
  
  f = FX(x)
  W = weil_divisor(f)
  I1, I2 = components(W)
  @test I == I1 || I == I2
  @test I != I1 || I != I2
  
  D = weil_divisor(E)
  D2 = irreducible_decomposition(D - W)
  @test integral(intersect(D2, E)) == -1
  
  @test integral(intersect(Oscar.move_divisor(D), E)) == -1
end

