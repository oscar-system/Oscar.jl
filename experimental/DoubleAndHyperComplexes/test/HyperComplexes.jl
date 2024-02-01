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

  HCx = Oscar.hyper_complex(Cx, auto_extend=true)
  HCy = Oscar.hyper_complex(Cy, auto_extend=true)
  HCz = Oscar.hyper_complex(Cz, auto_extend=true)

  @test HCx[(0,)] === R1
  @test HCx[(1,)] === R1
  @test map(HCx, 1, (1,))(R1[1]) == x*R1[1]
  @test iszero(HCx[(10,)])
  @test iszero(HCx[(-10,)])
  @test_throws AssertionError map(HCx, 2, (1,))
  @test_throws AssertionError map(HCx, 1, (1, 10))
  @test_throws AssertionError HCx[(1, 10)]

  HCx = Oscar.SimpleComplexWrapper(Cx, auto_extend=true)
  HCy = Oscar.SimpleComplexWrapper(Cy, auto_extend=true)
  HCz = Oscar.SimpleComplexWrapper(Cz, auto_extend=true)

  @test range(HCx) == range(Cx)
  @test map_range(HCx) == map_range(Cx)
  @test HCx[0] === R1
  @test HCx[1] === R1
  @test map(HCx, 1)(R1[1]) == x*R1[1]
  @test iszero(HCx[10])
  @test iszero(HCx[-10])
  @test_throws AssertionError map(HCx, 2, (1,))
  @test_throws AssertionError map(HCx, 1, (1, 10))
  @test_throws AssertionError HCx[(1, 10)]
  @test Oscar.direction(HCx) == :chain

  HCx = Oscar.hyper_complex(Cx, auto_extend=false)
  HCy = Oscar.hyper_complex(Cy, auto_extend=false)
  HCz = Oscar.hyper_complex(Cz, auto_extend=false)

  @test HCx[(0,)] === R1
  @test HCx[(1,)] === R1
  @test map(HCx, 1, (1,))(R1[1]) == x*R1[1]
  @test_throws ErrorException iszero(HCx[(10,)])
  @test_throws ErrorException iszero(HCx[(-10,)])
  @test_throws AssertionError map(HCx, 2, (1,))
  @test_throws AssertionError map(HCx, 1, (1, 10))
  @test_throws AssertionError HCx[(1, 10)]

  HCx = Oscar.SimpleComplexWrapper(Cx, auto_extend=false)
  HCy = Oscar.SimpleComplexWrapper(Cy, auto_extend=false)
  HCz = Oscar.SimpleComplexWrapper(Cz, auto_extend=false)

  @test HCx[0] === R1
  @test HCx[1] === R1
  @test map(HCx, 1)(R1[1]) == x*R1[1]
  @test_throws ErrorException map(HCx, 20)
  @test_throws ErrorException iszero(HCx[10])
  @test_throws ErrorException iszero(HCx[-10])
  @test_throws AssertionError map(HCx, 2, (1,))
  @test_throws AssertionError map(HCx, 1, (1, 10))
  @test_throws AssertionError HCx[(1, 10)]
  @test Oscar.direction(HCx) == :chain


  HC = tensor_product([Cx, Cy, Cz])

  @test HC[(0, 0, 0)] isa FreeMod
  @test iszero(HC[(-1, -1, -1)])
  @test iszero(map(HC, 1, (0, 0, 0)))
  @test !iszero(map(HC, 1, (1, 1, 1)))
  @test !iszero(map(HC, 2, (1, 1, 1)))
  @test !iszero(map(HC, 3, (1, 1, 1)))
end
  
@testset "Double complex wrapper" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  Z = FreeMod(R, 0)
  R1 = FreeMod(R, 1)
  Kx = ComplexOfMorphisms([hom(Z, R1, elem_type(R1)[]), hom(R1, R1, [x*R1[1]]), hom(R1, Z, [zero(Z)])], seed = -1)
  Ky = ComplexOfMorphisms([hom(Z, R1, elem_type(R1)[]), hom(R1, R1, [y*R1[1]]), hom(R1, Z, [zero(Z)])], seed = -1)
  Kxy = tensor_product(Kx, Ky)
  @test Kxy isa Oscar.DoubleComplexOfMorphisms

  K = Oscar.DoubleComplexWrapper(Kxy)
  @test K[0, 0] isa FreeMod
  @test can_compute_index(K, 1, 0)
  @test !has_index(K, 1, 1)
  @test vertical_map(K, 1, 0) isa ModuleFPHom
  @test !has_vertical_map(K, 1, 1)
  @test can_compute_vertical_map(K, 1, 1)
  @test !has_horizontal_map(K, 1, 1)
  @test can_compute_horizontal_map(K, 1, 1)
  @test vertical_direction(K) == :chain
  @test horizontal_direction(K) == :chain
  @test has_upper_bound(K)
  @test has_lower_bound(K)
  @test has_right_bound(K)
  @test has_left_bound(K)
  
  Kxy = tensor_product([Kx, Ky])
  @test Kxy isa Oscar.HyperComplex

  K = Oscar.DoubleComplexWrapper(Kxy)
  @test K[0, 0] isa FreeMod
  @test can_compute_index(K, 1, 0)
  @test !has_index(K, 1, 1)
  @test vertical_map(K, 1, 0) isa ModuleFPHom
  @test !has_vertical_map(K, 1, 1)
  @test can_compute_vertical_map(K, 1, 1)
  @test !has_horizontal_map(K, 1, 1)
  @test can_compute_horizontal_map(K, 1, 1)
  @test vertical_direction(K) == :chain
  @test horizontal_direction(K) == :chain
  @test has_upper_bound(K)
  @test has_lower_bound(K)
  @test has_right_bound(K)
  @test has_left_bound(K)
end
