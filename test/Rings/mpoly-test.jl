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
  S, (a, b, c) = PolynomialRing(QQ, ["a", "b", "c"])
  J = ideal(S, [(c^2+1)*(c^3+2)^2, b-c^2])
  r1 = c^2-b
  r2 = b^2*c+c^3+2*c^2+2
  L = gens(radical(J))

  @test jacobi_ideal(f) == ideal(R, [2*x, 2*y])
  @test jacobi_matrix(f) == matrix(R, 2, 1, [2*x, 2*y])
  @test jacobi_matrix(I) == matrix(R, 2, 2, [2*x, 4*x^3*y-y^3, 2*y, x^4-3*x*y^2])
  @test length(L) == 2
  @test length(findall(x->x==r1, L)) == 1
  @test length(findall(x->x==r2, L)) == 1

  @test issubset(ideal(S, [a]), ideal(S, [a]))
  @test issubset(ideal(S, [a]), ideal(S, [a, b]))
  @test !issubset(ideal(S, [c]), ideal(S, [b]))
  @test !issubset(ideal(S, [a, b, c]), ideal(S, [a*b*c]))
end

@testset "Primary decomposition" begin

  # primary_decomposition
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = ideal(R, [x, y*z^2])
  for method in (:GTZ, :SY)
    j = ideal(R, [R(1)])
    for (q, p) in primary_decomposition(i, alg=method)
      j = intersect(j, q)
      @test isprimary(q)
      @test isprime(p)
      @test p == radical(q)
    end
    @test j == i
  end

  R, (a, b, c, d) = PolynomialRing(ZZ, ["a", "b", "c", "d"])
  i = ideal(R, [9, (a+3)*(b+3)])
  l = primary_decomposition(i)
  @test length(l) == 2

  # minimal_primes
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  i1 = ideal(R, [z^3+2, -z^2+y])
  i2 = ideal(R, [z^2+1, -z^2+y])
  l = minimal_primes(i)
  @test length(l) == 2
  @test l[1] == i1 && l[2] == i2 || l[1] == i2 && l[2] == i1

  l = minimal_primes(i, alg=:charSets)
  @test length(l) == 2
  @test l[1] == i1 && l[2] == i2 || l[1] == i2 && l[2] == i1

  R, (a, b, c, d) = PolynomialRing(ZZ, ["a", "b", "c", "d"])
  i = ideal(R, [R(9), (a+3)*(b+3)])
  i1 = ideal(R, [R(3), a])
  i2 = ideal(R, [R(3), b])
  l = minimal_primes(i)
  @test length(l) == 2
  @test l[1] == i1 && l[2] == i2 || l[1] == i2 && l[2] == i1

  # equidimensional_decomposition_weak
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = intersect(ideal(R, [z]), ideal(R, [x, y]))
  i = intersect(i, ideal(R, [x^2, z^2]))
  i = intersect(i, ideal(R, [x^5, y^5, z^5]))
  l = equidimensional_decomposition_weak(i)
  @test length(l) == 3
  @test l[1] == ideal(R, [z^4, y^5, x^5, x^3*z^3, x^4*y^4])
  @test l[2] == ideal(R, [y*z, x*z, x^2])
  @test l[3] == ideal(R, [z])

 # equidimensional_decomposition_radical
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  l = equidimensional_decomposition_radical(i)
  @test length(l) == 1
  @test l[1] == ideal(R, [z^2-y, y^2*z+z^3+2*z^2+2])
  
  # equidimensional_hull
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = intersect(ideal(R, [z]), ideal(R, [x, y]))
  i = intersect(i, ideal(R, [x^2, z^2]))
  i = intersect(i, ideal(R, [x^5, y^5, z^5]))
  @test equidimensional_hull(i) == ideal(R, [z])

  R, (a, b, c, d) = PolynomialRing(ZZ, ["a", "b", "c", "d"])
  i = ideal(R, [1326*a^2*d^5, 1989*a^2*c^5, 102*b^4*d^5, 153*b^4*c^5,
            663*a^2*c^5*d^5, 51*b^4*c^5*d^5, 78*a^2*d^15, 117*a^2*c^15,
            78*a^15*d^5, 117*a^15*c^5, 6*a^2*b^4*d^15, 9*a^2*b^4*c^15,
            39*a^2*c^5*d^15, 39*a^2*c^15*d^5, 6*a^2*b^15*d^5, 9*a^2*b^15*c^5,
            6*a^15*b^4*d^5, 9*a^15*b^4*c^5, 39*a^15*c^5*d^5, 3*a^2*b^4*c^5*d^15,
            3*a^2*b^4*c^15*d^5, 3*a^2*b^15*c^5*d^5, 3*a^15*b^4*c^5*d^5])
  @test equidimensional_hull(i) == ideal(R, [R(3)])

  # equidimensional_hull_radical
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  @test equidimensional_hull_radical(i) == ideal(R, [z^2-y, y^2*z+z^3+2*z^2+2])

end

@testset "Groebner" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  I = ideal([2*x+3*y+4*z-5,3*x+4*y+5*z-2])
  @test groebner_basis(I,:degrevlex) == [y+2*z-11, 3*x+4*y+5*z-2]
  @test groebner_basis(I,:degrevlex, complete_reduction = true) == [y+2*z-11, x-z+14]

end

@testset "Primary decomposition" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  I = ideal(R, [x, y*z^2])
  J = ideal(R, [x, y^2])
  L = primary_decomposition(I)
  Q = ideal(R, [R(1)])
  @test isprime(I) == false
  @test isprimary(I) == false
  @test isprime(J) == false
  @test isprimary(J) == true
  for (q, p) in L
    Q = intersect(Q, q)
    @test isprimary(q)
    @test isprime(p)
    @test p == radical(q)
  end
  @test Q == I
end
