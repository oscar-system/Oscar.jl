@testset "affine schemes" begin
  An = affine_space( QQ, 4 )

  R = ambient_ring(An)
  # test that we are caching.
  @test R == ambient_ring(An)
  x = gens(R)
  I = ideal([x[1] + x[2]^2 - x[3]^3, x[1]^2 - 1])

  A = Spec(R, I)
  f1 = x[1]^2 + x[2]^2 - 1
  f2 = x[3]
  D1 = localize(A, f1)
  D2 = localize(A, f2)
  ambient_ring(D2)
  k = base_ring(An)
  D12 = localize(D1, f2)
  D21 = localize(D2, f1)
  @test base_ring(D1) == k
  @test base_ring(D2) == k
  @test base_ring(D12) == k
  @test root(D1) == A
  @test root(D12) == A
  @test root(D21) == A
  @test D12 != D21
  I = defining_ideal(D12)

  # Check whether the fractional representation of homomorphisms work
  IA2 = affine_space( QQ, 2 )
  R = ambient_ring( IA2 )
  vars = gens(R)
  x = vars[1]
  y = vars[2]
  IA2x = localize( IA2, x^3 )
  IA2y = localize( IA2, y )
  IA2xy = localize( IA2x, y^4 )
  Q = FractionField(R)
  h = AffSchMorphism( IA2xy, IA2y, [Q(x), Q(y)] )
  @test pullback(h)
  
  # Check whether the localization at elements which are already units 
  # also works.
  IA2 = affine_space( QQ, 2 )
  R = ambient_ring( IA2 )
  vars = gens(R)
  x = vars[1]
  y = vars[2]
  IA2xy = localize( IA2, x*y )
  @test U = localize( IA2xy, x )
  @test defining_ideal(U)
end

@testset "Misc routines" begin
  using Oscar
  using Oscar.Misc
  R,vars = PolynomialRing( QQ, [:x,:y] )
  x = vars[1]
  y = vars[2]
  f = (x^2+y^2+1)
  I = ideal( R, [f^6] )
  P, projection = quo( R, I )
  Q, phi, t = add_variables( P, ["t"] )
  g2 = projection(y)
  g = projection(x)
  h = lift( g )
  @test radical_membership( g*projection(y), g^4*projection(y) )
end

