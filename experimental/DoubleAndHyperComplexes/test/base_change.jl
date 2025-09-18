@testset "base change" begin
  R, (x, y) = QQ[:x, :y]
  R2 = FreeMod(R, 2)
  C = koszul_complex(Oscar.KoszulComplex, x*R2[1] + y*R2[2])

  U = complement_of_point_ideal(R, [0, 0])

  L, loc = localization(R, U)

  CC, red = change_base_ring(loc, C)
  @test !is_zero(CC[0])
  @test !is_zero(CC[1])
  @test !is_zero(CC[2])
  a = compose(map(C, 1), red[0]) 
  b = compose(red[1], map(CC, 1))
  @test all(a(g) == b(g) for g in gens(C[1]))
  a = compose(map(C, 2), red[1]) 
  b = compose(red[2], map(CC, 2))
  @test all(a(g) == b(g) for g in gens(C[2]))

  F1 = FreeMod(ZZ, 2)
  v = 15*F1[1] + 21*F1[2]
  C = koszul_complex(Oscar.KoszulComplex, v)

  CC, red = change_base_ring(QQ, C)
  @test !is_zero(CC[0])
  @test !is_zero(CC[1])
  @test !is_zero(CC[2])
  a = compose(map(C, 1), red[0]) 
  b = compose(red[1], map(CC, 1))
  @test all(a(g) == b(g) for g in gens(C[1]))
  a = compose(map(C, 2), red[1]) 
  b = compose(red[2], map(CC, 2))
  @test all(a(g) == b(g) for g in gens(C[2]))
end

