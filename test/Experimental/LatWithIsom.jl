using Oscar.LWI

@testset "Constructors and accessors" begin
  A4 = root_lattice(:A, 4)
  agg = = automorphism_group_generators(A4, ambient_representation = false)
  agg_ambient = automorphism_group_generators(A4, ambient_representation = true)
  f = rand(agg)
  g_ambient = rand(agg_ambient)

  L = @inferred lattice_with_isometry(A4)
  @test isone(isometry(L))
  @test isone(ambient_isometry(L))
  @test isone(order_of_isometry(L))

  for func in [rank, charpoly, minpoly, genus, ambient_space, basis_matrix,
               gram_matrix, rational_span, det, scale, norm, is_integral,
               degree, is_even, discriminant, signature_tuple]
    out = @inferred func(L)
  end

  nf = multiplicative_order(f)
  @test_throws ArgumentError lattice_with_isometry(A4, f, nf+1)
  @test_throws ArgumentError lattice_with_isometry(A4, zero_matrix(QQ, 0, 0), -1)

  L2 = @inferred lattice_with_isometry(A4, f, ambient_representation = false)
  @test order_of_isometry(L2) == nf
  L2v = @inferred dual(L2)
  @test order_of_isometry(L2v) == nf
  @test ambient_isometry(L2v) == ambient_isometry(L2)
  
  L3 = @inferred lattice_with_isometry(A4, g_ambient, ambient_representation = true)
  @test order_of_isometry(L2) == multiplicative_order(g_ambient)

  L4 = @inferred rescale(L3, QQ(1//4))
  @test !is_integral(L4)
  @test order_of_isometry(L4) == order_of_isometry(L3)
  @test_throws ArgumentError dual(L4)
  @test ambient_isometry(lll(L4)) == ambient_isometry(L4)
end
