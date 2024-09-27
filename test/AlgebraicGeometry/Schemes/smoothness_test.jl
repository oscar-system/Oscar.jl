@testset "_is_smooth" begin
  R, (x, y, z) = QQ[:x, :y, :z]

  I1 = ideal(R, x*y)
  IA3 = spec(R)
  X1, _ = sub(IA3, I1)
  @test !Oscar._is_smooth(X1)
  @test !Oscar._is_smooth(X1; jacobian_cut_off=0)

  I2 = ideal(R, x*y^2)
  X2, _ = sub(IA3, I2)
  @test !Oscar._is_smooth(X2)
  @test !Oscar._is_smooth(X2; jacobian_cut_off=0)

  A = R[x y 0; 0 y z]
  I3 = ideal(R, minors(A, 2))
  X3, _ = sub(IA3, I3)
  @test !Oscar._is_smooth(X3)
  @test !Oscar._is_smooth(X3; jacobian_cut_off=0)

  A = R[x y 0; 0 y-1 z]
  I4 = ideal(R, minors(A, 2))
  X4, _ = sub(IA3, I4)
  @test !Oscar._is_smooth(X4)
  @test !Oscar._is_smooth(X4; jacobian_cut_off=0)

  I5 = ideal(R, [x, y]) * ideal(R, [x, y-1, z])
  X5, _ = sub(IA3, I5)
  @test Oscar._is_smooth(X5)
  @test Oscar._is_smooth(X5; jacobian_cut_off=0)

  I6 = ideal(R, [x, y^2]) * ideal(R, [x, y, z-1])
  X6, _ = sub(IA3, I6)
  @test !Oscar._is_smooth(X6)
  @test !Oscar._is_smooth(X6; jacobian_cut_off=0)

  I7 = ideal(R, [x, y]) * ideal(R, [x-1])
  X7, _ = sub(IA3, I7)
  @test Oscar._is_smooth(X7)
  @test Oscar._is_smooth(X7; jacobian_cut_off=0)

  I5 = ideal(R, [x, y]) * ideal(R, [x, y-1, z^2])
  X5, _ = sub(IA3, I5)
  @test !Oscar._is_smooth(X5)
  @test !Oscar._is_smooth(X5; jacobian_cut_off=0)
end

