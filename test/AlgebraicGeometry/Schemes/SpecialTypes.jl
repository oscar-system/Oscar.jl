@testset "tests for SpecialTypes.jl" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  X = spec(R)
  h = x^2 - y^2 + z^2 - 1
  U = PrincipalOpenSubset(X, h)
  @test complement_equations(U)[1] == h
  @test U[1] == U

  inc = inclusion_morphism(U, check=true)
  @test inc === inclusion_morphism(U) # test caching
  @test inc isa Oscar.PrincipalOpenEmbedding

  f = OO(U)(x)/OO(U)(h^2)
  @test OO(U)(Oscar.generic_fraction(f, U)) == f

  @test domain(Oscar.underlying_morphism(inc)) === U
  @test codomain(Oscar.underlying_morphism(inc)) === X
  @test pullback(Oscar.underlying_morphism(inc)) == pullback(inc)
  @test h in ideal(OO(X), complement_equations(inc))

  V = hypersurface_complement(X, x)
  Y = subscheme(V, h)
  inc = inclusion_morphism(Y, V)
  clemb = ClosedEmbedding(inc, ideal(OO(V), h))
  clemb2 = ClosedEmbedding(V, ideal(OO(V), x^7*h))

  @test clemb == clemb2
  @test !(clemb === clemb2)
  @test !(domain(clemb) === domain(clemb2))
  @test complement(clemb) == complement(clemb2)

  @test ideal_type(OO(V)) == typeof(ideal(OO(V), one(OO(V))))
  @test ideal_type(OO(X)) == typeof(ideal(OO(X), one(OO(X))))
  @test ideal_type(OO(Y)) == typeof(ideal(OO(Y), one(OO(Y))))
  Z = subscheme(X, h)
  @test ideal_type(OO(Z)) == typeof(ideal(OO(Z), one(OO(Z))))
end
