@testset "Modules: orderings" begin
  R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
  F = FreeMod(R,2)
  a = (1 + 2*w + 3*x + 4*y + 5*z)*(F[1] + F[2])

  @test length(string(terms(a))) > 2

  o = lex(R)*lex(F)

  @test collect(coefficients(a; ordering = o)) ==
    [2, 2, 3, 3, 4, 4, 5, 5, 1, 1]

  @test collect(coefficients_and_exponents(a; ordering = o)) ==
    [(2, ([1, 0, 0, 0], 1)), (2, ([1, 0, 0, 0], 2)),
     (3, ([0, 1, 0, 0], 1)), (3, ([0, 1, 0, 0], 2)),
     (4, ([0, 0, 1, 0], 1)), (4, ([0, 0, 1, 0], 2)),
     (5, ([0, 0, 0, 1], 1)), (5, ([0, 0, 0, 1], 2)),
     (1, ([0, 0, 0, 0], 1)), (1, ([0, 0, 0, 0], 2))]

  @test collect(exponents(a; ordering = o)) == 
    [([1, 0, 0, 0], 1), ([1, 0, 0, 0], 2),
     ([0, 1, 0, 0], 1), ([0, 1, 0, 0], 2),
     ([0, 0, 1, 0], 1), ([0, 0, 1, 0], 2),
     ([0, 0, 0, 1], 1), ([0, 0, 0, 1], 2),
     ([0, 0, 0, 0], 1), ([0, 0, 0, 0], 2)]

  @test collect(terms(a; ordering = o)) ==
    [2*w*F[1], 2*w*F[2], 3*x*F[1], 3*x*F[2], 4*y*F[1], 4*y*F[2],
     5*z*F[1], 5*z*F[2], F[1], F[2]]

  @test collect(monomials(a; ordering = o)) ==
    [w*F[1], w*F[2], x*F[1], x*F[2], y*F[1], y*F[2],
     z*F[1], z*F[2], F[1], F[2]]

  @test leading_coefficient(a; ordering = o) == first(coefficients(a; ordering = o))
  @test leading_coefficient_and_exponent(a; ordering = o) == first(coefficients_and_exponents(a; ordering = o))
  @test leading_exponent(a; ordering = o) == first(exponents(a; ordering = o))
  @test leading_term(a; ordering = o) == first(terms(a; ordering = o))
  @test leading_monomial(a; ordering = o) == first(monomials(a; ordering = o))

  t = tail(a)
  @test a == leading_term(a) + t

  @test collect(terms(a; ordering = lex([w,x])*revlex(F)*lex([y,z]))) ==
    [2*w*F[2], 2*w*F[1], 3*x*F[2], 3*x*F[1],
     4*y*F[2], 5*z*F[2], F[2], 4*y*F[1], 5*z*F[1], F[1]]

  @test_throws ErrorException induced_ring_ordering(revlex(F))
  @test induced_ring_ordering(lex([w,x])*revlex(F)*lex([y,z])) == lex([w,x,y,z])
end
