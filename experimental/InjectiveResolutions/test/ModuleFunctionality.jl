@testset "module functionality for MonoidAlgebras" begin
  # definition of monoid algebra as quotient of polynomial ring
  S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
  J = ideal(S, [x*z-y^2])
  R_Q, phi = quo(S, J)

  # get MonoidAlgebra
  kQ = Oscar.MonoidAlgebra(R_Q)

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
# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x*z-y^2])
R_Q, phi = quo(S, J)

# get MonoidAlgebra
kQ = Oscar.MonoidAlgebra(R_Q)

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
end

@testset "quotient_ring_as_module" begin
# definition of polynomial ring k[x,y]
R_Q, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]; weights=[[1, 0], [0, 1]])

# get MonoidAlgebra
kQ = Oscar.MonoidAlgebra(R_Q)

# define ideal over monoid algebra
I = ideal(kQ, [x^4, x^2*y^2, y^4])

M = quotient_ring_as_module(I)
res = free_resolution(M)
end
