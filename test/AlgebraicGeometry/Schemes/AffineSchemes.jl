@testset "affine schemes" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  A3 = spec(R)
  deepcopy(A3)
  set_name!(A3, "𝔸³")
  @test iszero(Oscar.saturated_ideal(defining_ideal(A3)))
  @test iszero(defining_ideal(A3))
  f = x*y-z^2
  I = ideal(R, f)
  J = ideal(R, [f, x])
  A3empty = subscheme(A3,ideal(R,R(1)))
  absempty = EmptyScheme(QQ)
  @test (A3empty==absempty)
  @test (absempty==A3empty)
  X = subscheme(A3, I)
  @test defining_ideal(X) isa Oscar.MPolyIdeal
  @test_broken !is_non_zero_divisor(f,X)
  @test is_non_zero_divisor(f,A3)
  Xsub = subscheme(A3,J)
  @test !is_open_embedding(X,A3)
  @test is_subscheme(Xsub,X)
  @test !is_subscheme(X, Xsub)
  set_name!(X, "X")
  @test iszero(OO(X)(f))
  U = hypersurface_complement(A3, x)
  @test defining_ideal(U) isa Oscar.MPolyLocalizedIdeal
  UX = hypersurface_complement(X,x)
  @test is_non_zero_divisor(f,U)
  @test !is_non_zero_divisor(f,UX)
  @test !is_open_embedding(UX,U)
  @test_broken is_closed_embedding(UX,U)
  @test !is_open_embedding(UX,A3)
  @test is_subscheme(UX,U)
  @test is_subscheme(UX,X)
  @test is_subscheme(U,A3)
  @test !is_subscheme(A3,U)
  line = subscheme(A3, x)
  line_complement = hypersurface_complement(A3,x)
  @test !is_subscheme(line,line_complement)
  @test !is_subscheme(line_complement,line)
  @test !is_subscheme(line, X)
  @test is_subscheme(Xsub, line)
  U1 = hypersurface_complement(A3, [x])
  @test ambient_coordinate_ring(U) === R
  set_name!(U, "U")
  UX = intersect(X, U)
  @test defining_ideal(UX) isa Oscar.MPolyLocalizedIdeal
  set_name!(UX, "U ∩ X")
  @test is_subscheme(UX, X)
  @test is_subscheme(UX, U)
  @test name(UX) == "U ∩ X"
  @test X == closure(UX, A3)
  @test is_open_embedding(UX, X)
  @test is_closed_embedding(X, A3)
  UZ = subscheme(UX, y^2)
  subscheme(UX, [y^2])
  Z = subscheme(X, y^2)
  @test closure(UZ, X)==Z

  S, (u,v) = QQ[:u, :v]
  A2 = spec(S)
  set_name!(A2, "𝔸²")
  @test OO(UX)(y//z) == OO(UX)(z//x)
  phi = morphism(UX, A2, [y//z, z])
  L = subscheme(A2, u-v)
  phi_L = preimage(phi, L)
  @test OO(phi_L)(y//z) == OO(phi_L)(z)
  psi = restrict(phi, phi_L, L)
  Gamma_psi, p, q = graph(psi)
  @test iszero(pullback(p)(OO(phi_L)(y//z)) - pullback(q)(OO(L)(v)))

  Xstd = Oscar.standard_spec(X)
  mirr = morphism(Xstd, Xstd, [y, x, z])
  @test is_isomorphism(mirr)
  @test pullback(compose(inverse(mirr), mirr))(OO(Xstd)(x^2-34*z)) == OO(Xstd)(x^2-34*z+ f^2)
  @test is_empty(EmptyScheme(QQ))
  @test is_subscheme(EmptyScheme(QQ),A3)
  @test is_subscheme(EmptyScheme(QQ),U)
  @test !is_subscheme(U,EmptyScheme(QQ))
  @test is_subscheme(X,A3)
  @test !is_subscheme(A3, X)
  @test is_subscheme(A3,A3)
  @test is_subscheme(intersect(A3,A3), A3)
end

# Tests for dimension when localizing with respect to either a prime
# ideal or powers of an element
@testset "dimensions of affine schemes" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  A3 = spec(R)
  X = subscheme(A3, x*y)
  U = hypersurface_complement(A3, z)
  @test dim(A3) == 3
  @test dim(U) == 3
  @test dim(X) == 2
  @test codim(A3) == 0
  @test codim(X) == 1
  disjoint_plane_and_line = subscheme(A3, [x*(x - 1), x*y])
  line = hypersurface_complement(disjoint_plane_and_line, x)
  plane = hypersurface_complement(disjoint_plane_and_line, y)
  @test dim(line) == 1
  @test dim(plane) == 2
  A3_localized_along_line = spec(localization(R, complement_of_prime_ideal(ideal(R, [x, y])))[1])
  @test dim(A3_localized_along_line) == 2
  @test dim(Oscar.standard_spec(A3_localized_along_line)) == 2

  S = complement_of_point_ideal(R, [1, 1, 1])
  I = ideal(R, [x-1, y-1])*ideal(R, z)
  L, _ = localization(R, S)
  W, _ = quo(L, L(I))
  @test dim(spec(W)) == 1
end

@testset "dimensions of affine schemes over the integers" begin
  R, (x,y,z) = ZZ[:x, :y, :z]
  A3 = spec(R)
  X = subscheme(A3, x*y)
  U = hypersurface_complement(A3, z)
  @test dim(A3) == 4
  @test dim(U) == 4
  @test codim(A3) == 0
  @test codim(X) == 1
  disjoint_plane_and_line = subscheme(A3, [x*(x - 1), x*y])
  line = hypersurface_complement(disjoint_plane_and_line, x)
  plane = hypersurface_complement(disjoint_plane_and_line, y)
  @test dim(line) == 2
  @test dim(plane) == 3
  A3_localized_along_line = spec(localization(R, complement_of_prime_ideal(ideal(R, [x, y])))[1])
  @test dim(A3_localized_along_line) == 2
  @test dim(Oscar.standard_spec(A3_localized_along_line)) == 2
  P = ideal(R, R(5))
  @test is_prime(P)
  S = complement_of_prime_ideal(P)
  Y = spec(R, S)
  @test dim(Y) == 1
  Q = ideal(R, R.([7, x, y, z]))
  S = complement_of_prime_ideal(Q)
  Z = spec(R, S)
  @test dim(Z) == 4

  S = complement_of_prime_ideal(ideal(R, [x-1, y-1, z-1]))
  I = ideal(R, [x-1, y-1])*ideal(R, z)
  L, _ = localization(R, S)
  W, _ = quo(L, L(I))
  @test dim(spec(W)) == 1
end

@testset "dimensions of affine schemes over quotients of the integers" begin
  kk, _ = quo(ZZ, 4)
  R, (x,y,z) = kk[:x, :y, :z]
  A3 = spec(R)
  X = subscheme(A3, x*y)
  U = hypersurface_complement(A3, z)
  @test dim(A3) == 3
  @test dim(U) == 3
  @test codim(A3) == 0
  @test codim(X) == 1
  disjoint_plane_and_line = subscheme(A3, [x*(x - 1), x*y])
  line = hypersurface_complement(disjoint_plane_and_line, x)
  plane = hypersurface_complement(disjoint_plane_and_line, y)
  @test dim(line) == 1
  @test dim(plane) == 2
  # The other tests from above do not run, because the singular side does not digest the rings.
end


@testset "smoothness tests" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  M = R[x y+1 z-2 x+y; y z-3 x-y+5 y-z; z x-y z-2 x+y+z]
  I = ideal(R, minors(M, 3))
  Q, _ = quo(R, I)
  X = spec(Q)
  @test is_smooth(X)
end

@testset "AbsAffineScheme interface" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  A3 = spec(R)
  set_name!(A3, "𝔸³")
  f = x*y-z^2
  I = ideal(R, f)
  S = powers_of_element(x)
  spec(R,S)
  S = complement_of_prime_ideal(I)
  spec(R,I)
  X = subscheme(A3, I)
  set_name!(X, "X")
  @test iszero(OO(X)(f))
  U = hypersurface_complement(A3, x)
  @test is_subscheme(intersect(A3,U),intersect(U,A3))
  V  = PrincipalOpenSubset(U)
  @test dim(V) == 3
  @test is_subscheme(U,V)
  @test is_subscheme(V,V)

  @test_broken is_subscheme(intersect(A3,V),intersect(V,A3))
  @test is_subscheme(intersect(U,A3),V)
  @test_broken is_subscheme(intersect(A3,V),A3)
  @test ambient_coordinate_ring(V)===R
  @test ambient_coordinate_ring(U) === R
  @test Oscar.ring_type(V) == typeof(OO(V))
  @test Oscar.base_ring_type(typeof(V)) == typeof(QQ)
  @test base_ring(V) == QQ
  @test is_subscheme(V,U)
  @test is_subscheme(U,V)
  @test is_subscheme(V,A3)
  @test !is_subscheme(A3, V)
  h = gens(OO(V))[1]
  hypersurface_complement(V, h)
  hypersurface_complement(V, [h])
  inclusion_morphism(X,A3)
  inclusion_morphism(V,A3)
  inclusion_morphism(U,A3)
  inclusion_morphism(A3,A3)
end

@testset "spec ZZ" begin
  spec(ZZ)
  spec(QQ)
end

@testset "fiber product" begin
  R, _ = QQ[:x, :t]
  S, _ = QQ[:y, :t]
  T, _ = polynomial_ring(QQ,[:t])
  X = Oscar.standard_spec(spec(R))
  Y = Oscar.standard_spec(spec(S))
  B = Oscar.standard_spec(spec(T))
  phi1 = hom(OO(B), OO(X), [gens(OO(X))[2]])
  phi2 = hom(OO(B), OO(Y), [gens(OO(Y))[2]])
  Phi1 = morphism(X, B, phi1)
  Phi2 = morphism(Y, B, phi2)
  Z = fiber_product(Phi1, Phi2)[1]
  A = ambient_coordinate_ring(Z)
  a = gens(A)
  fib = subscheme(spec(A), ideal(A, [a[2]-a[4]]))
  @test fib==Z
end

@testset "ClosedEmbedding" begin
  R, (x, y) = QQ[:x, :y]
  X = spec(R)
  h = x^2 + y^2 -1
  I = ideal(R, [h])
  Y, inc = sub(X, I)
  @test inc isa ClosedEmbedding
  @test Y == subscheme(X, I)
  @test pullback(inc)(x) == OO(Y)(x)
  @test pullback(inc)(y) == OO(Y)(y)
  @test image_ideal(inc) == modulus(OO(Y))
  @test complement(inc) == AffineSchemeOpenSubscheme(X, [h])
end

@testset "fix for is_smooth" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  I1 = ideal(R, [x,y]);
  I2 = ideal(R, [z]);
  I3 = ideal(R, [x, z-1])
  I = intersect(I1, I2, I3)
  X = spec(R, I)
  @test !is_smooth(X)
end

@testset "reduction mod p" begin
  IA2 = affine_space(ZZ, 2)
  P = OO(IA2)
  x, y = gens(P)

  X, inc_X = sub(IA2, x^2 + y^2)

  U = hypersurface_complement(IA2, x)
  inc_U = inclusion_morphism(U, IA2)

  V = hypersurface_complement(X, x)
  inc_V = inclusion_morphism(V, X)

  kk, pr = quo(ZZ, 5)
  IA2_red, phi1 = base_change(pr, IA2)
  red_X, phi2, _ = base_change(pr, inc_X; codomain_map=phi1)
  X_red = domain(red_X)
  @test ambient_coordinate_ring(IA2_red) === ambient_coordinate_ring(X_red)
  red_U, phi3, _ = base_change(pr, inc_U; codomain_map=phi1)
  U_red = domain(red_U)
  @test ambient_coordinate_ring(IA2_red) === ambient_coordinate_ring(U_red)
  red_V, phi4, _ = base_change(pr, inc_V; codomain_map=red_X)
  V_red = domain(red_V)
  @test ambient_coordinate_ring(IA2_red) === ambient_coordinate_ring(V_red)

  m1 = compose(inclusion_morphism(V_red, IA2_red), phi1);
  m2 = compose(red_V, inclusion_morphism(V, IA2));
  @test m1 == m2

  # Testing morphisms
  inc_U = inclusion_morphism(V, X)
  (a, b, c) = base_change(pr, inc_U)
  @test compose(a, inc_U) == compose(b, c)

  # Behavior for special types
  U = PrincipalOpenSubset(IA2, y)
  UU, f = base_change(pr, U)
  @test UU isa PrincipalOpenSubset

  W = AffineSchemeOpenSubscheme(IA2, [x, y])
  WW, f = base_change(pr, W)
  @test WW isa AffineSchemeOpenSubscheme
end

@testset "principal open embeddings" begin
  IA3 = affine_space(QQ, [:x, :y, :z])
  (x, y, z) = gens(OO(IA3))
  U, inc = complement(IA3, y)
  @test [y] == complement_equations(inc)
  phi = inverse_on_image(inc)
  X, X_to_IA3 = sub(IA3, [x-y])
  X_simp = simplify(X)
  X_simp_to_X, X_to_X_simp = Oscar.identification_maps(X_simp)
  V, inc_V = complement(X_simp, OO(X_simp)[1])
  g = compose(inc_V, X_simp_to_X)
  there = Oscar.PrincipalOpenEmbedding(g, OO(X).([x]))
  back = inverse_on_image(there)
end

@testset "fiber products with principal open embeddings" begin
  IA3 = affine_space(QQ, [:x, :y, :z])
  (x, y, z) = gens(OO(IA3))
  U, inc = complement(IA3, y)
  @test [y] == complement_equations(inc)
  phi = inverse_on_image(inc)
  X, X_to_IA3 = sub(IA3, [x-y])
  X_simp = simplify(X)
  X_simp_to_X, X_to_X_simp = Oscar.identification_maps(X_simp)
  V, inc_V = complement(X_simp, OO(X_simp)[1])
  g = compose(inc_V, X_simp_to_X)
  there = Oscar.PrincipalOpenEmbedding(g, OO(X).([x]))

  # The classical case of two maps
  prod1, p1, p2 = fiber_product(there, there)
  id_V = id_hom(V)
  h = Oscar.induced_map_to_fiber_product(id_V, id_V, there, there, fiber_product=(prod1, p1, p2))
  @test codomain(h) === prod1
  @test is_isomorphism(h)

  # one open inclusion and a closed embedding
  Y, inc_Y = sub(X, z)
  pro, p1, p2 = fiber_product(there, inc_Y)
  h = Oscar.induced_map_to_fiber_product(p1, p2, there, inc_Y, fiber_product=(pro, p1, p2))
  @test codomain(h) === pro
  @test is_isomorphism(h)
  @test compose(p1, there) == compose(p2, inc_Y)

  # the same the other way around
  pro, p1, p2 = fiber_product(inc_Y, there)
  h = Oscar.induced_map_to_fiber_product(p1, p2, inc_Y, there, fiber_product=(pro, p1, p2))
  @test codomain(h) === pro
  @test is_isomorphism(h)
  @test codomain(p1) === domain(inc_Y)
  @test codomain(p2) === domain(there)
  @test compose(p1, inc_Y) == compose(p2, there)

  # two open inclusions
  pro, p1, p2 = fiber_product(there, there)
  h = Oscar.induced_map_to_fiber_product(p1, p2, there, there, fiber_product=(pro, p1, p2))
  @test codomain(h) === pro
  @test is_isomorphism(h)
  @test compose(p1, there) == compose(p2, there)
end

@testset "degree of zero-dimensional affine schemes" begin
  R, (x,) = polynomial_ring(QQ, [:x])
  I = ideal(R, x^2 - 1)
  X = spec(R)
  Q, _ = quo(R, I)
  Y = spec(Q)
  @test dim(Y) == 0
  @test degree(Y) == 2

  R, (x,y) = QQ[:x,:y]
  I = ideal(R, [x^2 + y^2 - 1, 2x^2 + (1//2)y^2 - 1])
  X = spec(R)
  Y,_ = sub(X, I)
  @test dim(Y) == 0
  @test degree(Y) == 4
end

@testset "normalization" begin
  R, (x,y) = polynomial_ring(QQ,[:x,:y])
  # non-normal
  I = ideal(R, x^3-y^2)
  X1 = spec(R,I)
  @test !is_normal(X1; check=false)
  N1 = normalization(X1)
  @test length(N1) == 1
  @test is_normal(N1[1][1]; check=false)
  @test is_smooth(N1[1][1])
  # normal
  J = ideal(R, x)
  X2 = spec(R, J)
  N2 = normalization(X2)
  @test length(N2) == 1
  @test is_smooth(N2[1][1])
  # reducible
  K = ideal(R, x^4-y^2*x)
  X3 = spec(R,K)
  N3 = normalization(X3)
  @test all(is_smooth(i[1]) for i in N3)
  @test length(N3) == 2

  # normalization in positive characteristic
  R, (x,y) = polynomial_ring(GF(5),[:x,:y])
  # non-normal
  I = ideal(R, x^3-y^2)
  X1 = spec(R,I)
  N1 = normalization(X1)
  @test length(N1) == 1
  @test is_smooth(N1[1][1])
  # normal
  J = ideal(R, x)
  X2 = spec(R, J)
  N2 = normalization(X2)
  @test length(N2) == 1
  @test is_smooth(N2[1][1])
  # reducible
  K = ideal(R, x^4-y^2*x)
  X3 = spec(R,K)
  N3 = normalization(X3)
  @test all(is_smooth(i[1]) for i in N3)
  @test length(N3) == 2

  # data corruption bug, throws an error if bug reappears
  P,(a,b) = polynomial_ring(QQ,[:a,:b])
  I = ideal(P,a^2-b)
  Q = normalization(quo(P,I)[1])[1][1]
  R = base_ring(Q)
  J = ideal(R, R[1])
  Sat = saturation(modulus(Q),J)
end
  
@testset "irreducible components" begin
  A = affine_space(QQ,2)
  (x,y) = coordinates(A)
  Y1 = subscheme(A, [x*y])
  @test length(irreducible_components(Y1))==2
  Y2 = hypersurface_complement(Y1,x)
  @test length(irreducible_components(Y2))==1
  P = ideal([x,y])
  l,_ = localization(OO(Y1), complement_of_prime_ideal(P))
  Y3 = spec(l)
  @test length(irreducible_components(Y3))==2
end 
