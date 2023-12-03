@testset "tor for modules over localizations at prime ideals" begin
  R, (x0, x1, x2, x3) = polynomial_ring(QQ, ["x0", "x1", "x2", "x3"])
  p = ideal(R, [x0-x1,x2,x3])
  Rp, _ = localization(R, complement_of_prime_ideal(p))
  I = ideal(Rp, [x0-x1,x2-x3])
  Ym = ideal(Rp, [x0^3+(x2^2-x3^2)*x0-x1^3+x2^3+x3^3])
  Rp1 = FreeMod(Rp, 1)
  IRp1, inc_I = I*Rp1
  YmRp1, inc_Ym = Ym*Rp1
  M1 = cokernel(inc_I)
  M2 = cokernel(inc_Ym)
  H0, _ = tor(M1, M2, 0)
  @test length(H0) == 3
end

