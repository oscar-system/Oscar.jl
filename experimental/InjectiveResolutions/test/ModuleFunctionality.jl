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
prune_with_map(M)

end

@testset "quotient_ring_as_module" begin
# definition of polynomial ring k[x,y]
R_Q, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]; weights=[[1, 0], [0, 1]])

# get MonoidAlgebra
kQ = Oscar.MonoidAlgebra(R_Q)
x, y = gens(kQ)

# define ideal over monoid algebra
I = ideal(kQ, [x^4, x^2*y^2, y^4])

M = quotient_ring_as_module(I)
res = free_resolution(M)
prune_with_map(M)

end

@testset "monomial bases" begin
# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x*z-y^2])
R_Q, phi = quo(S, J)

# get a MonoidAlgebra with a quotient ring as internal
A = Oscar.MonoidAlgebra(R_Q)
GA = grading_group(A)
@test all(parent(a) === A for a in monomial_basis(A, GA[1]+5*GA[2]))
A.algebra[1]*A[1] # Test promotion
@test is_normal(A)

# definition of polynomial ring k[x,y]
R_Q, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]; weights=[[1, 0], [0, 1]])

# get MonoidAlgebra with a polynomial ring as internal
B = Oscar.MonoidAlgebra(R_Q)
GB = grading_group(B)
@test all(parent(a) === B for a in monomial_basis(B, GB[1]+5*GB[2]))
B.algebra[1] *B[1] # Test promotion
@test is_normal(B)

end

@testset "coefficients for SubquoModules" begin
  # definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x*z-y^2])
R_Q, phi = quo(S, J)

# get a MonoidAlgebra with a quotient ring as internal
A = Oscar.MonoidAlgebra(R_Q)
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
  coefficients(I, f, A)
end

end

