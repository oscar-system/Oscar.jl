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

