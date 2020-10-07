@testset "Polynomial Orderings" begin
  R, (x, y, z) = PolynomialRing(QQ, 3)
  f = x*y + 5*z^3

  @test collect(monomials(f, :lex)) == [ x*y, z^3 ]
  @test collect(monomials(f, :revlex)) == [ z^3, x*y ]
  @test collect(monomials(f, :deglex)) == [ z^3, x*y ]
  @test collect(monomials(f, :degrevlex)) == [ z^3, x*y ]
  @test collect(monomials(f, :neglex)) == [ z^3, x*y ]
  @test collect(monomials(f, :negrevlex)) == [ x*y, z^3 ]
  @test collect(monomials(f, :negdeglex)) == [ x*y, z^3 ]
  @test collect(monomials(f, :negdegrevlex)) == [ x*y, z^3 ]

  w = [ 1, 2, 1 ]
  @test collect(monomials(f, :weightlex, w)) == [ x*y, z^3 ]
  @test collect(monomials(f, :weightrevlex, w)) == [ x*y, z^3 ]

  M = [ 1 1 1; 1 0 0; 0 1 0 ]
  @test collect(monomials(f, M)) == collect(monomials(f, :deglex))

  @test collect(terms(f, :deglex)) == [ 5z^3, x*y ]
  @test collect(exponent_vectors(f, :deglex)) == [ [ 0, 0, 3 ], [ 1, 1, 0 ] ]
  @test collect(coeffs(f, :deglex)) == [ QQ(5), QQ(1) ]

  Fp = FiniteField(7)
  R, (x, y, z) = PolynomialRing(Fp, 3, ordering = :deglex)
  f = x*y + 5*z^3
  @test collect(monomials(f, :lex)) == [ x*y, z^3 ]
  @test Oscar.leading_monomial(f, :lex) == x*y
  @test Oscar.leading_coeff(f, :lex) == Fp(1)
  @test Oscar.leading_term(f, :lex) == x*y

end

@testset "Polynomial homs" begin
  R, (x, y) = PolynomialRing(QQ, ["x", "y"])
  I1 = x^2 + y^2
  I2 = x^2 * y^2
  I3 = x^3*y - x*y^3
  S, (a,b,c) = PolynomialRing(QQ, ["a", "b", "c"])
  h = hom(S, R, [I1, I2, I3])
  @test kernel(h) == ideal(S, [a^2*b - 4*b^2 - c^2])
  @test h(gen(S, 1)) == I1
  @test image(h, ideal(S, [a,b])) == ideal(R, [I1, I2])
  @test preimage(h, ideal(R, [I2, I3])) == ideal(S, [b, c])
end

@testset "Ideal operations" begin
  R, (x, y) = PolynomialRing(QQ, ["x", "y"])
  f = x^2 + y^2
  g = x^4*y - x*y^3
  I = [f, g]
  @test jacobi_ideal(f) == ideal(R, [2*x, 2*y])
  @test jacobi_matrix(f) == matrix(R, 2, 1, [2*x, 2*y])
  @test jacobi_matrix(I) == matrix(R, 2, 2, [2*x, 4*x^3*y-y^3, 2*y, x^4-3*x*y^2])
end

@testset "Groebner" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  I = ideal([2*x+3*y+4*z-5,3*x+4*y+5*z-2])
  @test groebner_basis(I,:degrevlex) == [y+2*z-11, 3*x+4*y+5*z-2]
  @test groebner_basis(I,:degrevlex, complete_reduction = true) == [y+2*z-11, x-z+14]
  
end
