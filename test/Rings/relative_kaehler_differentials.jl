@testset "relative kaehler differentials" begin
  R, (x, y, z, t) = QQ[:x, :y, :z, :t]
  A = R[x z-t; z+t y]
  I = ideal(R, det(A))
  OX, pr_OX = quo(R, I)
  S, (u, v) = graded_polynomial_ring(R, [:u, :v])

  Aext = hcat(map_entries(S, A), S[u; v])
  J = ideal(S, minors(Aext, 2))
  SY, pr_SY = quo(S, J)

  X = spec(OX)
  Y = projective_scheme(SY)

  MX1 = Oscar.relative_kaehler_differentials(OX; base_variables=[t])
  MX1 = Oscar.relative_kaehler_differentials(OX; base_variables=[t, y])

  U = hypersurface_complement(X, t)
  MU1 = Oscar.relative_kaehler_differentials(OO(U); base_variables=[t])
  MU1 = Oscar.relative_kaehler_differentials(OO(U); base_variables=[t, y])
  
  Oscar._relative_kaehler_differentials(S, [4])
  Oscar._relative_kaehler_differentials(graded_polynomial_ring(OX, [:u, :v])[1], [4])
  Oscar._relative_kaehler_differentials(SY, [4])
end

