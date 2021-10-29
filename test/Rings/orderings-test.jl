@testset "Polynomial Orderings construction" begin
   R, (x, y, z) = PolynomialRing(QQ, 3)
   f = x*y + 5*z^3
 
   @test isa(lex([x, y]), MonomialOrdering)
   @test isa(revlex([x, y, z]), MonomialOrdering)
   @test isa(degrevlex([x, z]), MonomialOrdering)
   @test isa(deglex([x, y, z]), MonomialOrdering)
   @test isa(neglex([x, y, z]), MonomialOrdering)
   @test isa(negrevlex([x, y, z]), MonomialOrdering)
   @test isa(negdeglex([x, y, z]), MonomialOrdering)
   @test isa(negdegrevlex([x, y, z]), MonomialOrdering)
 
   @test isa(wdeglex([x, y], [1, 2]), MonomialOrdering)
   @test isa(wdegrevlex([x, y], [1, 2]), MonomialOrdering)
   @test isa(negwdeglex([x, y], [1, 2]), MonomialOrdering)
   @test isa(negwdegrevlex([x, y], [1, 2]), MonomialOrdering)
 
   @test isa(revlex([x, y])*neglex([z]), MonomialOrdering)
 
   @test collect(monomials(f, :lex)) == [ x*y, z^3 ]
   @test collect(monomials(f, :revlex)) == [ z^3, x*y ]
   @test collect(monomials(f, :deglex)) == [ z^3, x*y ]
   @test collect(monomials(f, :degrevlex)) == [ z^3, x*y ]
   @test collect(monomials(f, :neglex)) == [ z^3, x*y ]
   @test collect(monomials(f, :negrevlex)) == [ x*y, z^3 ]
   @test collect(monomials(f, :negdeglex)) == [ x*y, z^3 ]
   @test collect(monomials(f, :negdegrevlex)) == [ x*y, z^3 ]
 
   w = [ 1, 2, 1 ]
   @test collect(monomials(f, :wdeglex, w)) == [ x*y, z^3 ]
   @test collect(monomials(f, :wdegrevlex, w)) == [ x*y, z^3 ]
   @test collect(monomials(f, :negwdeglex, w)) == [ x*y, z^3 ]
   @test collect(monomials(f, :negwdegrevlex, w)) == [ x*y, z^3 ]
 
   M = [ 1 1 1; 1 0 0; 0 1 0 ]
   @test collect(monomials(f, M)) == collect(monomials(f, :deglex))
end

@testset "Polynomial Orderings terms, monomials and coefficients" begin
   R, (x, y, z) = PolynomialRing(QQ, 3)
   f = x*y + 5*z^3
 
   @test collect(terms(f, :deglex)) == [ 5z^3, x*y ]
   @test collect(exponent_vectors(f, :deglex)) == [ [ 0, 0, 3 ], [ 1, 1, 0 ] ]
   @test collect(coefficients(f, :deglex)) == [ QQ(5), QQ(1) ]

   Fp = FiniteField(7)
   R, (x, y, z) = PolynomialRing(Fp, 3, ordering = :deglex)
   f = x*y + 5*z^3
   @test collect(monomials(f, :lex)) == [ x*y, z^3 ]
   @test Oscar.leading_monomial(f, :lex) == x*y
   @test Oscar.leading_coefficient(f, :lex) == Fp(1)
   @test Oscar.leading_term(f, :lex) == x*y
end
 
 @testset "Polynomial Orderings comparison" begin
   R, (x, y, z) = @inferred PolynomialRing(QQ, ["x", "y", "z"])

   @test lex([x])*lex([y,z]) == lex([x, y, z])
   @test lex([z])*lex([y])*lex([x]) == revlex([x, y, z])
   @test degrevlex([x, y, z])*revlex([y]) == degrevlex([x, y, z])
   @test deglex([z])*deglex([x])*deglex([y]) == lex([z])*lex([x, y])
   @test neglex([x, y])*neglex([z]) == neglex([x, y, z])
   @test deglex([x, y, z]) == wdeglex([x, y, z], [1, 1, 1])
   @test negdeglex([x, y, z]) == negwdeglex([x, y, z], [1, 1, 1])
   @test negdegrevlex([x, y, z]) == negwdegrevlex([x, y, z], [1, 1, 1])
   @test neglex([z])*neglex([y])*neglex([x]) == negrevlex([x, y, z])
 end

   
