@testset "base change" begin
  R, (x, y) = QQ[:x, :y]
  R2 = FreeMod(R, 2)
  C = koszul_complex(Oscar.KoszulComplex, x*R2[1] + y*R2[2])

  U = complement_of_point_ideal(R, [0, 0])

  L, loc = localization(R, U)

  change_base_ring(loc, C)
end
