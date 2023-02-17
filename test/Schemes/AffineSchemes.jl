@testset "affine schemes" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  A3 = Spec(R)
  deepcopy(A3)
  set_name!(A3, "𝔸³")
  @test iszero(Oscar.ambient_closure_ideal(A3))
  f = x*y-z^2
  I = ideal(R, f)
  J = ideal(R, [f, x])
  A3empty = subscheme(A3,ideal(R,R(1)))
  absempty = EmptyScheme(QQ)
  @test (A3empty==absempty)
  @test (absempty==A3empty)
  X = subscheme(A3, I)
  @test_broken !is_non_zero_divisor(f,X)
  @test is_non_zero_divisor(f,A3)
  Xsub = subscheme(A3,J)
  @test !is_open_embedding(X,A3)
  @test issubset(Xsub,X)
  @test !issubset(X, Xsub)
  set_name!(X, "X")
  @test iszero(OO(X)(f))
  U = hypersurface_complement(A3, x)
  UX = hypersurface_complement(X,x)
  @test is_non_zero_divisor(f,U)
  @test !is_non_zero_divisor(f,UX)
  @test !is_open_embedding(UX,U)
  @test_broken is_closed_embedding(UX,U)
  @test !is_open_embedding(UX,A3)
  @test issubset(UX,U)
  @test issubset(UX,X)
  @test issubset(U,A3)
  @test !issubset(A3,U)
  line = subscheme(A3, x)
  line_complement = hypersurface_complement(A3,x)
  @test !issubset(line,line_complement)
  @test !issubset(line_complement,line)
  @test !issubset(line, X)
  @test issubset(Xsub, line)
  U1 = hypersurface_complement(A3, [x])
  @test ambient_coordinate_ring(U) === R
  set_name!(U, "U")
  UX = intersect(X, U)
  set_name!(UX, "U ∩ X")
  @test issubset(UX, X)
  @test issubset(UX, U)
  @test name(UX) == "U ∩ X"
  @test X == closure(UX, A3)
  @test is_open_embedding(UX, X)
  @test is_closed_embedding(X, A3)
  UZ = subscheme(UX, y^2)
  subscheme(UX, [y^2])
  Z = subscheme(X, y^2)
  @test closure(UZ, X)==Z
  
  S, (u,v) = QQ["u", "v"]
  A2 = Spec(S)
  set_name!(A2, "𝔸²")
  @test OO(UX)(y//z) == OO(UX)(z//x)
  phi = SpecMor(UX, A2, [y//z, z])
  L = subscheme(A2, u-v)
  phi_L = preimage(phi, L)
  @test OO(phi_L)(y//z) == OO(phi_L)(z)
  psi = restrict(phi, phi_L, L)
  Gamma_psi, p, q = graph(psi)
  @test iszero(pullback(p)(OO(phi_L)(y//z)) - pullback(q)(OO(L)(v)))
  
  Xstd = Oscar.standard_spec(X)
  mirr = SpecMor(Xstd, Xstd, [y, x, z])
  @test is_isomorphism(mirr)
  @test pullback(compose(inverse(mirr), mirr))(OO(Xstd)(x^2-34*z)) == OO(Xstd)(x^2-34*z+ f^2)
  @test is_empty(EmptyScheme(QQ))
  @test issubset(EmptyScheme(QQ),A3)
  @test issubset(EmptyScheme(QQ),U)
  @test !issubset(U,EmptyScheme(QQ))
  @test issubset(X,A3)
  @test !issubset(A3, X)
  @test issubset(A3,A3)
  @test issubset(intersect(A3,A3), A3)
  @test dim(A3) == 3
  @test dim(U) == 3
  @test dim(X) == 2
  @test codim(A3) == 0
  @test codim(X) == 1
end

@testset "smoothness tests" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  M = R[x y+1 z-2 x+y; y z-3 x-y+5 y-z; z x-y z-2 x+y+z]
  I = ideal(R, minors(M, 3))
  Q, _ = quo(R, I)
  X = Spec(Q)
  @test is_smooth(X)
end

@testset "AbsSpec interface" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  A3 = Spec(R)
  set_name!(A3, "𝔸³")
  f = x*y-z^2
  I = ideal(R, f)
  S = powers_of_element(x)
  Spec(R,S)
  S = complement_of_prime_ideal(I)
  Spec(R,I)
  X = subscheme(A3, I)
  set_name!(X, "X")
  @test iszero(OO(X)(f))
  U = hypersurface_complement(A3, x)
  @test issubset(intersect(A3,U),intersect(U,A3))
  V  = PrincipalOpenSubset(U)
  @test dim(V) == 3
  @test issubset(U,V)
  @test issubset(V,V)

  @test_broken issubset(intersect(A3,V),intersect(V,A3))
  @test issubset(intersect(U,A3),V)
  @test_broken issubset(intersect(A3,V),A3)
  @test ambient_coordinate_ring(V)===R
  @test ambient_coordinate_ring(U) === R
  @test ring_type(V) == typeof(OO(V))
  @test base_ring_type(typeof(V)) == typeof(QQ)
  @test base_ring_elem_type(V) == fmpq
  @test base_ring(V) == QQ
  @test issubset(V,U)
  @test issubset(U,V)
  @test issubset(V,A3)
  @test !issubset(A3, V)
  h = gens(OO(V))[1]
  hypersurface_complement(V, h)
  hypersurface_complement(V, [h])
  inclusion_morphism(X,A3)
  inclusion_morphism(V,A3)
  inclusion_morphism(U,A3)
  inclusion_morphism(A3,A3)
end

@testset "Spec ZZ" begin
  Spec(ZZ)
  Spec(QQ)
end

@testset "fiber product" begin
  R, _ = QQ["x","t"]
  S, _ = QQ["y","t"]
  T, _ = PolynomialRing(QQ,["t"])
  X = Oscar.standard_spec(Spec(R))
  Y = Oscar.standard_spec(Spec(S))
  B = Oscar.standard_spec(Spec(T))
  phi1 = hom(OO(B), OO(X), [gens(OO(X))[2]])
  phi2 = hom(OO(B), OO(Y), [gens(OO(Y))[2]])
  Phi1 = SpecMor(X, B, phi1)
  Phi2 = SpecMor(Y, B, phi2)
  Z = fiber_product(Phi1, Phi2)[1]
  A = ambient_coordinate_ring(Z)
  a = gens(A)
  fib = subscheme(Spec(A), ideal(A, [a[2]-a[4]]))
  @test fib==Z
end

@testset "ClosedEmbedding" begin
  R, (x, y) = QQ["x", "y"]
  X = Spec(R)
  h = x^2 + y^2 -1
  I = ideal(R, [h])
  Y, inc = sub(X, I)
  @test inc isa ClosedEmbedding
  @test Y == subscheme(X, I)
  @test pullback(inc)(x) == OO(Y)(x)
  @test pullback(inc)(y) == OO(Y)(y)
  @test image_ideal(inc) == modulus(OO(Y))
  @test complement(inc) == SpecOpen(X, [h])
end
