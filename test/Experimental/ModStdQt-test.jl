@testset "Examples.ModStdQt" begin
  Qt, t = PolynomialRing(QQ, :t=>1:2)
  Qtx, x = PolynomialRing(FractionField(Qt), :x => 1:3)

  f1 = x[1]^2*x[2]^3*x[3] + 2*t[1]*x[1]*x[2]*x[3]^2 + 7*x[2]^3
  f2 = x[1]^2*x[2]^4*x[3] + (t[1]-7*t[2])*x[1]^2*x[2]*x[3]^2 - x[1]*x[2]^2*x[3]^2+2*x[1]^2*x[2]*x[3] - 12*x[1] + t[2]*x[2]
  f3 = (t[1]^2+t[2]-2)*x[2]^5*x[3] + (t[1]+5*t[2])*x[1]^2*x[2]^2*x[3] - t[2]*x[1]*x[2]^3*x[3] - t[2]*x[1]*x[2]^3*x[3] - x[1]*x[2]^3+x[2]^4+2*t[1]^2*x[2]^2*x[3]
  f4 = t[1]*x[2]^2*x[2]^2*x[3] - x[1]*x[2]^3 *x[3] + (-t[1]+4)*x[2]^3*x[3]^2 + 3*t[1]*x[1]*x[2]*x[3]^3 + 4*x[3]^2 - t[2]*x[1]

  #z = groebner_basis(ideal(Qtx, [f1, f2, f3, f4]))
  #@test length(z) == 4
  #z = groebner_basis(ideal(Qtx, [f2, f2, f3, f4]), ord = :lex)
  #@test length(z) == 3
end


@testset "Examples.factor_absolute" begin
  S, (a, b) = PolynomialRing(QQ, ["a", "b"]);
  S = parent(a//b)
  St, (u, v) = PolynomialRing(S, ["u", "v"])
  f = factor_absolute(u^2+v^2*a)
  @test length(f) == 2
end
