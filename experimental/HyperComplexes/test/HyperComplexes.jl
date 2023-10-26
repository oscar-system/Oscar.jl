@testset "hypercomplexes and tensor products" begin
  R, (x, y, z) = QQ[:x, :y, :z]

  R1 = FreeMod(R, 1)
  mx = hom(R1, R1, [x*R1[1]])
  my = hom(R1, R1, [y*R1[1]])
  mz = hom(R1, R1, [z*R1[1]])

  R0 = FreeMod(R, 0)
  inc0 = hom(R0, R1, elem_type(R1)[])
  pr0 = hom(R1, R0, [zero(R0)])

  Cx = ComplexOfMorphisms([inc0, mx, pr0], seed = -1)
  Cy = ComplexOfMorphisms([inc0, my, pr0], seed = -1)
  Cz = ComplexOfMorphisms([inc0, mz, pr0], seed = -1)

  HCx = Oscar.hyper_complex(Cx)
  HCy = Oscar.hyper_complex(Cy)
  HCz = Oscar.hyper_complex(Cz)

  @test HCx[(0,)] === R1
  @test HCx[(1,)] === R1
  @test map(HCx, 1, (1,))(R1[1]) == x*R1[1]
  @test iszero(HCx[(10,)])
  @test iszero(HCx[(-10,)])
  @test_throws AssertionError map(HCx, 2, (1,))
  @test_throws AssertionError map(HCx, 1, (1, 10))
  @test_throws AssertionError HCx[(1, 10)]

  HC = tensor_product([Cx, Cy, Cz])

  @test HC[(0, 0, 0)] isa FreeMod
  @test iszero(HC[(-1, -1, -1)])
  @test iszero(map(HC, 1, (0, 0, 0)))
  @test !iszero(map(HC, 1, (1, 1, 1)))
  @test !iszero(map(HC, 2, (1, 1, 1)))
  @test !iszero(map(HC, 3, (1, 1, 1)))
end
