@testset "Polynomial Orderings construction" begin
   R, (x, y, z) = PolynomialRing(QQ, 3)
   f = x*y + 5*z^3
   g = (1 + x + y + z)^2
 
   @test isa(lex([x, y]), MonomialOrdering)
   @test isa(revlex([x, y, z]), MonomialOrdering)
   @test isa(degrevlex([x, z]), MonomialOrdering)
   @test isa(deglex([x, y, z]), MonomialOrdering)
   @test isa(neglex([x, y, z]), MonomialOrdering)
   @test isa(negrevlex([x, y, z]), MonomialOrdering)
   @test isa(negdeglex([x, y, z]), MonomialOrdering)
   @test isa(negdegrevlex([x, y, z]), MonomialOrdering)
   @test isa(matrix_ordering([x, y, z], matrix(ZZ, [1 1 1; 1 0 0; 0 1 0])), MonomialOrdering)

   @test isa(wdeglex([x, y], [1, 2]), MonomialOrdering)
   @test isa(wdegrevlex([x, y], [1, 2]), MonomialOrdering)
   @test isa(negwdeglex([x, y], [1, 2]), MonomialOrdering)
   @test isa(negwdegrevlex([x, y], [1, 2]), MonomialOrdering)
 
   @test isa(revlex([x, y])*neglex([z]), MonomialOrdering)

   @test collect(monomials(g, :lex)) == [x^2, x*y, x*z, x, y^2, y*z, y, z^2, z, 1] # lp
   @test collect(monomials(g, :revlex)) == [z^2, y*z, x*z, z, y^2, x*y, y, x^2, x, 1] # rp
   @test collect(monomials(g, :deglex)) == [x^2, x*y, x*z, y^2, y*z, z^2, x, y, z, 1] # Dp
   @test collect(monomials(g, :degrevlex)) == [x^2, x*y, y^2, x*z, y*z, z^2, x, y, z, 1] # dp
   @test collect(monomials(g, :neglex)) == [1, z, z^2, y, y*z, y^2, x, x*z, x*y, x^2] # ls
   @test collect(monomials(g, :negrevlex)) == [1, x, x^2, y, x*y, y^2, z, x*z, y*z, z^2] # rs not documented ?
   @test collect(monomials(g, :negdeglex)) == [1, x, y, z, x^2, x*y, x*z, y^2, y*z, z^2] # Ds
   @test collect(monomials(g, :negdegrevlex)) == [1, x, y, z, x^2, x*y, y^2, x*z, y*z, z^2] # ds

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

   @test !is_global(lex([x, y]))
   @test !is_local(lex([x, y]))
   @test !is_mixed(lex([x, y]))
   @test !is_global(neglex([y, z]))
   @test !is_local(neglex([y, z]))
   @test !is_mixed(neglex([y, z]))
   @test !is_global(lex([x, y])*neglex([z]))
   @test !is_local(lex([x, y])*neglex([z]))
   @test is_mixed(lex([x, y])*neglex([z]))

   a = neglex([x, y, z])*degrevlex([x, y, z])
   @test !is_global(a)
   @test is_local(a)
   @test !is_mixed(a)

   a = degrevlex([x, y, z])*neglex([x, y, z])
   @test is_global(a)
   @test !is_local(a)
   @test !is_mixed(a)

   a = neglex([x, y])*degrevlex([y, z])
   @test !is_global(a)
   @test !is_local(a)
   @test is_mixed(a)
end

@testset "Polynomial Orderings terms, monomials and coefficients" begin
   R, (x, y, z) = PolynomialRing(QQ, 3)
   f = x*y + 5*z^3
 
   @test collect(terms(f, :deglex)) == [ 5z^3, x*y ]
   @test collect(exponent_vectors(f, :deglex)) == [ [ 0, 0, 3 ], [ 1, 1, 0 ] ]
   @test collect(coefficients(f, :deglex)) == [ QQ(5), QQ(1) ]

   Fp = GF(7)
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

   m = matrix(ZZ, [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1])
   @test negdeglex(gens(R)) == matrix_ordering(gens(R), m)

   m = matrix(ZZ, [-2 -3 -4; 1 0 0; 0 1 0; 0 0 1])
   @test negwdeglex(gens(R), [2, 3, 4]) == matrix_ordering(gens(R), m)
 end

 @testset "Polynomial Orderings sorting" begin
   R, (x1, x2, x3, x4) = PolynomialRing(QQ, "x".*string.(1:4))
   
   M = [x2^3, x1*x2^2, x1^2*x2, x2^2*x4, x2^2*x3, x2^2, x1^3,
        x1*x2*x4, x1*x2*x3, x1*x2, x1^2*x4, x1^2*x3, x1^2, x2*x4^2, x2*x3*x4,
        x2*x4, x2*x3^2, x2*x3, x2, x1*x4^2, x1*x3*x4, x1*x4, x1*x3^2, x1*x3,
        x1, x4^3, x3*x4^2, x4^2, x3^2*x4, x3*x4, x4, x3^3, x3^2, x3, one(R)]

   f = sum(M)
   o = wdeglex([x1, x2], [1, 2])*revlex([x3, x4])
   @test collect(monomials(f, o)) == M

   M = [x4^3, x3*x4^2, x3^2*x4, x4^2, x3^3, x3*x4, x3^2, x4, x3,
        one(R), x1*x4^2, x1*x3*x4, x1*x3^2, x1*x4, x1*x3, x1, x1^2*x4, x1^2*x3,
        x1^2, x1^3, x2*x4^2, x2*x3*x4, x2*x3^2, x2*x4, x2*x3, x2, x1*x2*x4,
        x1*x2*x3, x1*x2, x1^2*x2, x2^2*x4, x2^2*x3, x2^2, x1*x2^2, x2^3]

   f = sum(M)
   o = negrevlex([x1, x2])*wdegrevlex([x3, x4], [1, 2])
   @test collect(monomials(f, o)) == M
     
   M = [one(R), x3, x3^2, x4, x3^3, x3*x4, x3^2*x4, x4^2, x3*x4^2,
        x4^3, x2, x2*x3, x2*x3^2, x2*x4, x2*x3*x4, x2*x4^2, x2^2, x2^2*x3,
        x2^2*x4, x2^3, x1, x1*x3, x1*x3^2, x1*x4, x1*x3*x4, x1*x4^2, x1*x2,
        x1*x2*x3, x1*x2*x4, x1*x2^2, x1^2, x1^2*x3, x1^2*x4, x1^2*x2, x1^3]
   
   f = sum(M)
   o = neglex([x1, x2])*negwdegrevlex([x3, x4], [1, 2])
   @test collect(monomials(f, o)) == M
        
   M = [x3^3, x3^2*x4, x3^2, x3*x4^2, x3*x4, x3, x4^3, x4^2, x4,
        one(R), x1*x3^2, x1*x3*x4, x1*x3, x1*x4^2, x1*x4, x1, x1^2*x3, x1^2*x4,
        x1^2, x2*x3^2, x2*x3*x4, x2*x3, x2*x4^2, x2*x4, x2, x1^3, x1*x2*x3,
        x1*x2*x4, x1*x2, x1^2*x2, x2^2*x3, x2^2*x4, x2^2, x1*x2^2, x2^3]

   f = sum(M)
   o = negwdeglex([x1, x2], [1, 2])*lex([x3, x4])
   @test collect(monomials(f, o)) == M

   M = [x1^3, x1^2*x2, x1*x2^2, x2^3, x1^2, x1^2*x3, x1^2*x4, x1*x2,
        x1*x2*x3, x1*x2*x4, x2^2, x2^2*x3, x2^2*x4, x1, x1*x3, x1*x4, x1*x3^2,
        x1*x3*x4, x1*x4^2, x2, x2*x3, x2*x4, x2*x3^2, x2*x3*x4, x2*x4^2, one(R),
        x3, x4, x3^2, x3*x4, x4^2, x3^3, x3^2*x4, x3*x4^2, x4^3]

   f = sum(M)
   o = deglex([x1, x2])*negdeglex([x3, x4])
   @test collect(monomials(f, o)) == M
   
   M = [one(R), x1, x2, x3, x4, x1^2, x1*x2, x2^2, x1*x3, x2*x3,
        x3^2, x1*x4, x2*x4, x3*x4, x4^2, x1^3, x1^2*x2, x1*x2^2, x2^3, x1^2*x3,
        x1*x2*x3, x2^2*x3, x1*x3^2, x2*x3^2, x3^3, x1^2*x4, x1*x2*x4, x2^2*x4,
        x1*x3*x4, x2*x3*x4, x3^2*x4, x1*x4^2, x2*x4^2, x3*x4^2, x4^3]

   f = sum(M)
   o = negdegrevlex(gens(R))
   @test collect(monomials(f, o)) == M

   o = matrix_ordering(gens(R), matrix(ZZ, [ 1 1 1 1; 0 0 0 -1; 0 0 -1 0; 0 -1 0 0 ]))
   @test collect(monomials(f, o)) == collect(monomials(f, degrevlex(gens(R))))
 end

@testset "Polynomial Ordering internal conversion to Singular" begin
   R, (x, y, s, t, u) = PolynomialRing(QQ, ["x", "y", "s", "t", "u"])

   O1 = degrevlex(gens(R))
   @test length(string(O1)) > 2
   @test string(singular(O1)) == "ordering_dp(5)"

   O2 = lex([x, y])*deglex([s, t, u])
   @test length(string(O2)) > 2
   @test string(singular(O2)) == "ordering_lp(2) * ordering_Dp(3)"

   O3 = wdeglex(gens(R), [2, 3, 5, 7, 3])
   @test length(string(O3)) > 2
   @test string(singular(O3)) == "ordering_Wp([2, 3, 5, 7, 3])"

   O4 = deglex([x, y, t]) * deglex([y, s, u])
   @test length(string(O4)) > 2
   @test string(singular(O4)) == "ordering_M([1 1 0 1 0; 0 -1 0 -1 0; 0 0 0 -1 0; 0 0 1 0 1; 0 0 0 0 -1])"

   K = FreeModule(R, 3)

   O5 = revlex(gens(K))*degrevlex(gens(R))
   @test length(string(O5)) > 2
   @test string(singular(O5)) == "ordering_c() * ordering_dp(5)"

   O6 = matrix_ordering([x, y], matrix(ZZ, 2, 2, [1 2; 3 4])) * lex(gens(K)) * wdeglex([s, t, u], [1, 2, 3])
   @test length(string(O6)) > 2
   @test string(singular(O6)) == "ordering_M([1 2; 3 4]) * ordering_C() * ordering_Wp([1, 2, 3])"

   O7 = weighted_ordering(gens(R),[-1,2,0,2,0])*degrevlex(gens(R))
   @test length(string(O7)) > 2
   @test string(singular(O7)) == "ordering_a([-1, 2, 0, 2, 0]) * ordering_dp(5)"
end

@testset "Polynomial Ordering misc bugs" begin
  R, (x, y) = QQ["x", "y"]
  @test degrevlex(gens(R)) != degrevlex(Oscar.reverse(gens(R)))

  R, (x, y, z) = QQ["x", "y", "z"]
  @test degrevlex(gens(R)) != degrevlex(Oscar.reverse(gens(R)))

  a = negwdegrevlex([z, x, y], [4, 5, 6])
  @test matrix_ordering([x, y, z], weight_matrix(a)) ==
        matrix_ordering([x, y, z], [-5 -6 -4; 0 -1 0; -1 0 0; 0 0 -1])

  a = weighted_ordering([y, z, x], [4, 6, 8])
  @test canonical_weight_matrix(a) == matrix(ZZ, 1, 3, [4, 2, 3])
end

