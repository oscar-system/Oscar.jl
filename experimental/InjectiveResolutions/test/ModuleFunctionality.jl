@testset "module functionality for MonoidAlgebras" begin
  # get MonoidAlgebra
  kQ = monoid_algebra([[0, 1], [1, 1], [2, 1]],QQ)

  F = FreeMod(kQ, 2)
  (x, y, z) = gens(kQ)
  I, inc_I = sub(F, [x*F[1], y*F[2]])
  u, v = gens(I)
  @test repres(v) in I
  c = coordinates(repres(u + v), I)
  @test c[1] == 1 && c[2] == 1

  phi = hom(F, F, [F[1], F[1]])
  @test F[1] - F[2] in kernel(phi)[1]
end

@testset "graded modules over MonoidAlgebras and their ideals" begin
  # get MonoidAlgebra
  kQ = monoid_algebra([[0, 1], [1, 1], [2, 1]],QQ)

  @test is_graded(kQ)
  G = grading_group(kQ)
  (x, y, z) = gens(kQ)
  degree.(gens(kQ))

  F = graded_free_module(kQ, [zero(G)])
  I, inc = sub(F, [x*F[1] for x in gens(kQ)])
  J = ideal(kQ, gens(kQ))
  Oscar._saturation(I.sub, J)

  M = quotient_ring_as_module(J)
  res = free_resolution(M; length=5)
  prune_with_map(M)
end

@testset "quotient_ring_as_module" begin
  # get MonoidAlgebra
  kQ = monoid_algebra([[1, 0], [0, 1]],QQ)
  x, y = gens(kQ)

  # define ideal over monoid algebra
  I = ideal(kQ, [x^4, x^2*y^2, y^4])

  M = quotient_ring_as_module(I)
  # res = free_resolution(M)
  prune_with_map(M)
end

@testset "monomial bases" begin
  # get a MonoidAlgebra with a quotient ring as internal
  A = monoid_algebra([[0, 1], [1, 1], [2, 1]],QQ)
  GA = grading_group(A)
  @test all(parent(a) === base_ring(A.algebra) for a in monomial_basis(A, GA[1]+5*GA[2]))
  A.algebra[1]*A[1] # Test promotion
  @test is_normal(A)

  # get MonoidAlgebra with a polynomial ring as internal
  B = monoid_algebra([[1, 0], [0, 1]],QQ)
  GB = grading_group(B)
  @test all(parent(a) === B.algebra for a in monomial_basis(B, GB[1]+5*GB[2]))
  B.algebra[1] *B[1] # Test promotion
  @test is_normal(B)
end

@testset "coefficients for SubquoModules" begin
  # get a MonoidAlgebra with a quotient ring as internal
  A = monoid_algebra([[0, 1], [1, 1], [2, 1]],QQ)
  GA = grading_group(A)
  monomial_basis(A, GA[1]+5*GA[2])
  A.algebra[1]*A[1] # Test promotion

  F = graded_free_module(A, [zero(GA)])
  I, inc = sub(F, [x*F[1] for x in gens(A)])
  mon_base = monomial_basis(A, GA[2])
  [x*I[1] for x in mon_base]
  [x*one(A) for x in mon_base]
  [x*F[1] for x in mon_base]

  for f in Oscar.faces(A)
    coefficients(I, f)
  end
end

@testset "workaround for issue 4833" begin
  # get a MonoidAlgebra with a quotient ring as internal
  A = monoid_algebra([[1,0],[0,1]],QQ)
  G = grading_group(A)

  F = graded_free_module(A, [zero(G)])
  I, inc = sub(F, push!([x*F[1] for x in gens(A)], zero(F)))
  J, inc = sub(F, [x*F[1] for x in gens(A)])

  coordinates(repres(J[1]), J)
  coordinates(repres(I[1]), I)
end

