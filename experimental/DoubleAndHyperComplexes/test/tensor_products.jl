@testset "tensor products of hyper complexes" begin
  S, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])

  # Hilbert-Burch example of an affine twisted cubic
  A = S[x y z; y-1 z-1 x-1]
  I = ideal(S, minors(A, 2))
  F = FreeMod(S, 1)
  IF, _ = I*F
  M, _ = quo(F, IF)
  res, _ = free_resolution(Oscar.SimpleFreeResolution, M)

  c = tensor_product(res, res)
  @test dim(c) == 2
  @test lower_bound(c, 1) == lower_bound(c, 2) == 0
  @test Oscar.has_upper_bound(c, 1) 
  @test Oscar.has_upper_bound(c, 2)
  @test rank(c[0, 0]) == 1
  @test rank(c[1, 0]) == rank(c[0, 1]) == 3
  @test rank(c[1, 1]) == rank(c[1, 1]) == 9
  @test rank(c[2, 1]) == rank(c[1, 2]) == 6
  @test rank(c[2, 2]) == 4
  @test iszero(c[10, 10])
  @test !Oscar.can_compute_index(c, (-1, -2))
  @test Oscar.can_compute_map(c, 2, (0, 1))
  @test map(c, 1, (1, 0)) isa FreeModuleHom
  @test map(c, 2, (0, 1)) isa FreeModuleHom

  # Check that all squares commute
  @test compose(map(c, 2, (1, 1)), map(c, 1, (1, 0))) == compose(map(c, 1, (1, 1)), map(c, 2, (0, 1)))
  @test compose(map(c, 2, (1, 2)), map(c, 1, (1, 1))) == compose(map(c, 1, (1, 2)), map(c, 2, (0, 2)))
  @test compose(map(c, 2, (2, 2)), map(c, 1, (2, 1))) == compose(map(c, 1, (2, 2)), map(c, 2, (1, 2)))
  @test compose(map(c, 2, (2, 1)), map(c, 1, (2, 0))) == compose(map(c, 1, (2, 1)), map(c, 2, (1, 1)))

  c2 = tensor_product(c, res)
  @test dim(c2) == 3
  @test lower_bound(c2, 1) == lower_bound(c2, 2) == lower_bound(c2, 3) == 0
  @test Oscar.has_upper_bound(c2, 1)
  @test Oscar.has_upper_bound(c2, 2) 
  @test Oscar.has_upper_bound(c2, 3) 
  @test rank(c2[0, 0, 0]) == 1
  @test rank(c2[0, 1, 0]) == rank(c2[0, 1, 0]) == 3
  @test rank(c2[0, 1, 1]) == rank(c2[1, 1, 0]) == 9
  @test rank(c2[0, 2, 1]) == rank(c2[1, 2, 0]) == 6
  @test rank(c2[0, 2, 2]) == 4

  @test rank(c2[1, 0, 0]) == 1*3
  @test rank(c2[1, 1, 0]) == rank(c2[0, 1, 1]) == 3*3
  @test rank(c2[1, 1, 1]) == rank(c2[1, 1, 1]) == 9*3
  @test rank(c2[1, 2, 1]) == rank(c2[1, 2, 1]) == 6*3
  @test rank(c2[1, 2, 2]) == 4*3

  @test rank(c2[2, 0, 0]) == 1*2
  @test rank(c2[2, 1, 0]) == rank(c2[0, 1, 2]) == 3*2
  @test rank(c2[2, 1, 1]) == rank(c2[1, 1, 2]) == 9*2
  @test rank(c2[2, 2, 1]) == rank(c2[1, 2, 2]) == 6*2
  @test rank(c2[2, 2, 2]) == 4*2

  # Check that some squares commute
  @test compose(map(c, 2, (1, 1, 0)), map(c, 1, (1, 0, 0))) == compose(map(c, 1, (1, 1, 0)), map(c, 2, (0, 1, 0)))
  @test compose(map(c, 2, (1, 2, 0)), map(c, 1, (1, 1, 0))) == compose(map(c, 1, (1, 2, 0)), map(c, 2, (0, 2, 0)))
  @test compose(map(c, 2, (2, 2, 0)), map(c, 1, (2, 1, 0))) == compose(map(c, 1, (2, 2, 0)), map(c, 2, (1, 2, 0)))
  @test compose(map(c, 2, (2, 1, 0)), map(c, 1, (2, 0, 0))) == compose(map(c, 1, (2, 1, 0)), map(c, 2, (1, 1, 0)))

  @test compose(map(c, 2, (1, 1, 1)), map(c, 1, (1, 0, 1))) == compose(map(c, 1, (1, 1, 1)), map(c, 2, (0, 1, 1)))
  @test compose(map(c, 2, (1, 2, 1)), map(c, 1, (1, 1, 1))) == compose(map(c, 1, (1, 2, 1)), map(c, 2, (0, 2, 1)))
  @test compose(map(c, 2, (2, 2, 1)), map(c, 1, (2, 1, 1))) == compose(map(c, 1, (2, 2, 1)), map(c, 2, (1, 2, 1)))
  @test compose(map(c, 2, (2, 1, 1)), map(c, 1, (2, 0, 1))) == compose(map(c, 1, (2, 1, 1)), map(c, 2, (1, 1, 1)))

  @test compose(map(c, 2, (1, 1, 2)), map(c, 1, (1, 0, 2))) == compose(map(c, 1, (1, 1, 2)), map(c, 2, (0, 1, 2)))
  @test compose(map(c, 2, (1, 2, 2)), map(c, 1, (1, 1, 2))) == compose(map(c, 1, (1, 2, 2)), map(c, 2, (0, 2, 2)))
  @test compose(map(c, 2, (2, 2, 2)), map(c, 1, (2, 1, 2))) == compose(map(c, 1, (2, 2, 2)), map(c, 2, (1, 2, 2)))
  @test compose(map(c, 2, (2, 1, 2)), map(c, 1, (2, 0, 2))) == compose(map(c, 1, (2, 1, 2)), map(c, 2, (1, 1, 2)))

  c3 = tensor_product(res, c)
  @test dim(c3) == 3
  @test lower_bound(c3, 1) == lower_bound(c3, 2) == lower_bound(c3, 3) == 0
  @test Oscar.has_upper_bound(c3, 1) 
  @test Oscar.has_upper_bound(c3, 2)
  @test Oscar.has_upper_bound(c3, 3)
  @test rank(c3[0, 0, 0]) == 1
  @test rank(c3[1, 0, 0]) == rank(c3[0, 1, 0]) == 3
  @test rank(c3[1, 1, 0]) == rank(c3[1, 1, 0]) == 9
  @test rank(c3[2, 1, 0]) == rank(c3[1, 2, 0]) == 6
  @test rank(c3[2, 2, 0]) == 4

  @test rank(c3[0, 0, 1]) == 1*3
  @test rank(c3[1, 0, 1]) == rank(c3[0, 1, 1]) == 3*3
  @test rank(c3[1, 1, 1]) == rank(c3[1, 1, 1]) == 9*3
  @test rank(c3[2, 1, 1]) == rank(c3[1, 2, 1]) == 6*3
  @test rank(c3[2, 2, 1]) == 4*3

  @test rank(c3[0, 0, 2]) == 1*2
  @test rank(c3[1, 0, 2]) == rank(c3[0, 1, 2]) == 3*2
  @test rank(c3[1, 1, 2]) == rank(c3[1, 1, 2]) == 9*2
  @test rank(c3[2, 1, 2]) == rank(c3[1, 2, 2]) == 6*2
  @test rank(c3[2, 2, 2]) == 4*2

  # The induced maps should be independent of the order of creation of the 
  # tensor products
  for i1 in 1:2
    for i2 in 1:2
      for i3 in 1:2
        @test matrix(map(c2, 1, (i1, i2, i3))) == matrix(map(c3, 1, (i1, i2, i3)))
        @test matrix(map(c2, 2, (i1, i2, i3))) == matrix(map(c3, 2, (i1, i2, i3)))
        @test matrix(map(c2, 3, (i1, i2, i3))) == matrix(map(c3, 3, (i1, i2, i3)))
      end
    end
  end


  c4 = tensor_product(res, res, res)
  @test dim(c4) == 3
  @test map(c4, 2, (1,1,1,1)) isa FreeModuleHom
  c5 = tensor_product(c4, c2)
  @test dim(c5) == 6
  @test map(c5, 3, (1,1,1,1,1,1)) isa FreeModuleHom
  @test map(c5, 5, (1,1,1,1,1,1)) isa FreeModuleHom
  @test map(c5, 1, (1,1,1,1,1,1)) isa FreeModuleHom
  @test map(c5, 6, (1,1,1,1,1,1)) isa FreeModuleHom
end
