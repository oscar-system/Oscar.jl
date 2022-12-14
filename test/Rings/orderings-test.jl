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

   @test collect(monomials(g; ordering = lex(R))) == [x^2, x*y, x*z, x, y^2, y*z, y, z^2, z, 1] # lp
   @test collect(monomials(g; ordering = revlex(R))) == [z^2, y*z, x*z, z, y^2, x*y, y, x^2, x, 1] # rp
   @test collect(monomials(g; ordering = deglex(R))) == [x^2, x*y, x*z, y^2, y*z, z^2, x, y, z, 1] # Dp
   @test collect(monomials(g; ordering = degrevlex(R))) == [x^2, x*y, y^2, x*z, y*z, z^2, x, y, z, 1] # dp
   @test collect(monomials(g; ordering = neglex(R))) == [1, z, z^2, y, y*z, y^2, x, x*z, x*y, x^2] # ls
   @test collect(monomials(g; ordering = negrevlex(R))) == [1, x, x^2, y, x*y, y^2, z, x*z, y*z, z^2] # rs not documented ?
   @test collect(monomials(g; ordering = negdeglex(R))) == [1, x, y, z, x^2, x*y, x*z, y^2, y*z, z^2] # Ds
   @test collect(monomials(g; ordering = negdegrevlex(R))) == [1, x, y, z, x^2, x*y, y^2, x*z, y*z, z^2] # ds

   @test collect(monomials(f; ordering = lex(R))) == [ x*y, z^3 ]
   @test collect(monomials(f; ordering = revlex(R))) == [ z^3, x*y ]
   @test collect(monomials(f; ordering = deglex(R))) == [ z^3, x*y ]
   @test collect(monomials(f; ordering = degrevlex(R))) == [ z^3, x*y ]
   @test collect(monomials(f; ordering = neglex(R))) == [ z^3, x*y ]
   @test collect(monomials(f; ordering = negrevlex(R))) == [ x*y, z^3 ]
   @test collect(monomials(f; ordering = negdeglex(R))) == [ x*y, z^3 ]
   @test collect(monomials(f; ordering = negdegrevlex(R))) == [ x*y, z^3 ]
 
   w = [ 1, 2, 1 ]
   @test collect(monomials(f; ordering = wdeglex(R, w))) == [ x*y, z^3 ]
   @test collect(monomials(f; ordering = wdegrevlex(R, w))) == [ x*y, z^3 ]
   @test collect(monomials(f; ordering = negwdeglex(R, w))) == [ x*y, z^3 ]
   @test collect(monomials(f; ordering = negwdegrevlex(R, w))) == [ x*y, z^3 ]
 
   M = [ 1 1 1; 1 0 0; 0 1 0 ]
   @test collect(monomials(f; ordering = matrix_ordering(R, M))) ==
                                  collect(monomials(f; ordering = deglex(R)))

   a = lex([x, y])
   @test is_global(a)
   @test !is_local(a)
   @test !is_mixed(a)
   @test cmp(a, x, one(R)) > 0
   @test cmp(a, y, one(R)) > 0
   @test_throws ErrorException cmp(a, z, one(R))

   a = neglex([y, z])
   @test !is_global(a)
   @test is_local(a)
   @test !is_mixed(a)
   @test_throws ErrorException cmp(a, x, one(R))
   @test cmp(a, y, one(R)) < 0
   @test cmp(a, z, one(R)) < 0

   a = lex([x, y])*neglex([z])
   @test !is_global(a)
   @test !is_local(a)
   @test is_mixed(a)
   @test cmp(a, x^5*y^6*z^3, x^5*y^6*z^4) > 0

   a = neglex([x, y, z])*degrevlex([x, y, z])
   @test !is_global(a)
   @test is_local(a)
   @test !is_mixed(a)

   a = degrevlex([x, y, z])*neglex([x, y, z])
   @test is_global(a)
   @test !is_local(a)
   @test !is_mixed(a)

   a = neglex([x, y])*degrevlex([y, z])
   @test a == neglex([x, y])*lex([z])
   @test !is_global(a)
   @test !is_local(a)
   @test is_mixed(a)

   @test_throws ArgumentError monomial_ordering(gens(R), :foo)
   @test_throws ArgumentError monomial_ordering(gens(R), :lex, ones(Int, ngens(R)+1))
   @test_throws ArgumentError monomial_ordering(gens(R), :foo, ones(Int, ngens(R)))
   @test_throws ArgumentError matrix_ordering(gens(R), zero_matrix(ZZ, 2, ngens(R) + 1))

   @test lex(R) == lex(gens(R))
   @test deglex(R) == deglex(gens(R))
   @test degrevlex(R) == degrevlex(gens(R))
   @test revlex(R) == revlex(gens(R))
   @test neglex(R) == neglex(gens(R))
   @test negrevlex(R) == negrevlex(gens(R))
   @test negdegrevlex(R) == negdegrevlex(gens(R))
   @test negdeglex(R) == negdeglex(gens(R))
   @test wdeglex(R, [1,2,3]) == wdeglex(gens(R), [1,2,3])
   @test wdegrevlex(R, [1,2,3]) == wdegrevlex(gens(R), [1,2,3])
   @test negwdeglex(R, [1,2,3]) == negwdeglex(gens(R), [1,2,3])
   @test negwdegrevlex(R, [1,2,3]) == negwdegrevlex(gens(R), [1,2,3])

   @test support(lex([x])) == [x]
   @test support(lex([x, y])) == [x,y] || support(lex([x,y])) == [y,x]
   @test 3 == length(support(deglex([x,y])*wdeglex([y,z], [1,2])))
end

@testset "Polynomial Orderings printing" begin
   R, (x, y, z) = PolynomialRing(QQ, 3)
   f = x*y + 5*z^3 + 2
   @test length(string(coefficients(f))) > 2
   @test length(string(coefficients_and_exponents(f))) > 2
   @test length(string(exponents(f))) > 2
   @test length(string(monomials(f))) > 2
   @test length(string(terms(f))) > 2
end

@testset "Polynomial Orderings terms, monomials and coefficients" begin
   R, (x, y, z) = PolynomialRing(QQ, 3)
   f = x*y + 5*z^3
 
   @test collect(terms(f; ordering = deglex(R))) == [ 5z^3, x*y ]
   @test collect(exponents(f; ordering = deglex(R))) == [ [ 0, 0, 3 ], [ 1, 1, 0 ] ]
   @test collect(coefficients(f; ordering = deglex(R))) == [ QQ(5), QQ(1) ]

   @test collect(coefficients_and_exponents(f)) ==
                    map(tuple, collect(coefficients(f)), collect(exponents(f)))

   @test leading_coefficient_and_exponent(f) ==
                   (leading_coefficient(f), leading_exponent(f))

   Fp = GF(7)
   R, (x, y, z) = PolynomialRing(Fp, 3, ordering = :deglex)
   f = x*y + 5*z^3
   @test collect(monomials(f; ordering = lex(R))) == [ x*y, z^3 ]
   @test leading_monomial(f; ordering = lex(R)) == x*y
   @test leading_coefficient(f; ordering = lex(R)) == Fp(1)
   @test leading_term(f; ordering = lex(R)) == x*y
   @test leading_exponent(f) == first(collect(exponents(f)))

   f = f^3 + f
   g = MPolyBuildCtx(R)
   for (c, e) in coefficients_and_exponents(f)
      push_term!(g, c, e)
   end
   @test f == finish(g)

   @test f == leading_term(f) + tail(f)
   @test f == leading_term(f; ordering = lex(R)) + tail(f; ordering = lex(R))
   @test f == leading_term(f; ordering = neglex(R)) + tail(f; ordering = neglex(R))
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

function test_opposite_ordering(a)
  R = base_ring(a)
  b = opposite_ordering(R, a)
  M = matrix(a)
  N = reduce(hcat, [M[:,i] for i in ncols(M):-1:1])
  @test b == matrix_ordering(gens(R), N)
  @test a == opposite_ordering(R, b)
end

@testset "Polynomial Orderings sorting" begin
   R, (x1, x2, x3, x4) = PolynomialRing(QQ, "x".*string.(1:4))
   
   M = [x2^3, x1*x2^2, x1^2*x2, x2^2*x4, x2^2*x3, x2^2, x1^3,
        x1*x2*x4, x1*x2*x3, x1*x2, x1^2*x4, x1^2*x3, x1^2, x2*x4^2, x2*x3*x4,
        x2*x4, x2*x3^2, x2*x3, x2, x1*x4^2, x1*x3*x4, x1*x4, x1*x3^2, x1*x3,
        x1, x4^3, x3*x4^2, x4^2, x3^2*x4, x3*x4, x4, x3^3, x3^2, x3, one(R)]

   f = sum(M)
   o = wdeglex([x1, x2], [1, 2])*revlex([x3, x4])
   test_opposite_ordering(o)
   @test collect(monomials(f; ordering = o)) == M
   for i in 2:length(M)
      @test cmp(o, M[i-1], M[i]) > 0
   end

   M = [x4^3, x3*x4^2, x3^2*x4, x4^2, x3^3, x3*x4, x3^2, x4, x3,
        one(R), x1*x4^2, x1*x3*x4, x1*x3^2, x1*x4, x1*x3, x1, x1^2*x4, x1^2*x3,
        x1^2, x1^3, x2*x4^2, x2*x3*x4, x2*x3^2, x2*x4, x2*x3, x2, x1*x2*x4,
        x1*x2*x3, x1*x2, x1^2*x2, x2^2*x4, x2^2*x3, x2^2, x1*x2^2, x2^3]

   f = sum(M)
   o = negrevlex([x1, x2])*wdegrevlex([x3, x4], [1, 2])
   test_opposite_ordering(o)
   @test collect(monomials(f; ordering = o)) == M
   for i in 2:length(M)
      @test cmp(o, M[i-1], M[i]) > 0
   end
     
   M = [one(R), x3, x3^2, x4, x3^3, x3*x4, x3^2*x4, x4^2, x3*x4^2,
        x4^3, x2, x2*x3, x2*x3^2, x2*x4, x2*x3*x4, x2*x4^2, x2^2, x2^2*x3,
        x2^2*x4, x2^3, x1, x1*x3, x1*x3^2, x1*x4, x1*x3*x4, x1*x4^2, x1*x2,
        x1*x2*x3, x1*x2*x4, x1*x2^2, x1^2, x1^2*x3, x1^2*x4, x1^2*x2, x1^3]
   
   f = sum(M)
   o = neglex([x1, x2])*negwdegrevlex([x3, x4], [1, 2])
   test_opposite_ordering(o)
   @test collect(monomials(f; ordering = o)) == M
   for i in 2:length(M)
      @test cmp(o, M[i-1], M[i]) > 0
   end
        
   M = [x3^3, x3^2*x4, x3^2, x3*x4^2, x3*x4, x3, x4^3, x4^2, x4,
        one(R), x1*x3^2, x1*x3*x4, x1*x3, x1*x4^2, x1*x4, x1, x1^2*x3, x1^2*x4,
        x1^2, x2*x3^2, x2*x3*x4, x2*x3, x2*x4^2, x2*x4, x2, x1^3, x1*x2*x3,
        x1*x2*x4, x1*x2, x1^2*x2, x2^2*x3, x2^2*x4, x2^2, x1*x2^2, x2^3]

   f = sum(M)
   o = negwdeglex([x1, x2], [1, 2])*lex([x3, x4])
   test_opposite_ordering(o)
   @test collect(monomials(f; ordering = o)) == M
   for i in 2:length(M)
      @test cmp(o, M[i-1], M[i]) > 0
   end

   M = [x1^3, x1^2*x2, x1*x2^2, x2^3, x1^2, x1^2*x3, x1^2*x4, x1*x2,
        x1*x2*x3, x1*x2*x4, x2^2, x2^2*x3, x2^2*x4, x1, x1*x3, x1*x4, x1*x3^2,
        x1*x3*x4, x1*x4^2, x2, x2*x3, x2*x4, x2*x3^2, x2*x3*x4, x2*x4^2, one(R),
        x3, x4, x3^2, x3*x4, x4^2, x3^3, x3^2*x4, x3*x4^2, x4^3]

   f = sum(M)
   o = deglex([x1, x2])*negdeglex([x3, x4])
   test_opposite_ordering(o)
   @test collect(monomials(f; ordering = o)) == M
   for i in 2:length(M)
      @test cmp(o, M[i-1], M[i]) > 0
   end
   
   M = [one(R), x1, x2, x3, x4, x1^2, x1*x2, x2^2, x1*x3, x2*x3,
        x3^2, x1*x4, x2*x4, x3*x4, x4^2, x1^3, x1^2*x2, x1*x2^2, x2^3, x1^2*x3,
        x1*x2*x3, x2^2*x3, x1*x3^2, x2*x3^2, x3^3, x1^2*x4, x1*x2*x4, x2^2*x4,
        x1*x3*x4, x2*x3*x4, x3^2*x4, x1*x4^2, x2*x4^2, x3*x4^2, x4^3]

   f = sum(M)
   o = negdegrevlex(gens(R))
   test_opposite_ordering(o)
   @test collect(monomials(f; ordering = o)) == M
   for i in 2:length(M)
      @test cmp(o, M[i-1], M[i]) > 0
   end

   o = matrix_ordering(gens(R), matrix(ZZ, [ 1 1 1 1; 0 0 0 -1; 0 0 -1 0; 0 -1 0 0 ]))
   test_opposite_ordering(o)
   @test collect(monomials(f; ordering = o)) == collect(monomials(f; ordering = degrevlex(gens(R))))

   for a in (matrix_ordering([x1, x2], [1 2; 3 4]),
             negwdegrevlex([x1, x2], [1, 2]),
             wdegrevlex([x1, x2], [1, 2]),
             negdeglex([x1, x2]),
             negdegrevlex([x1, x2]),
             revlex([x1, x2]),
             lex([x1, x2]))
      test_opposite_ordering(a*lex([x3, x4]))
      @test collect(monomials(x3 + x4; ordering = a*lex([x3, x4]))) == [x3, x4]
   end
end

@testset "Polynomial Ordering internal conversion to Singular" begin
   R, (x, y, s, t, u) = PolynomialRing(QQ, ["x", "y", "s", "t", "u"])

   for O in (wdegrevlex([x,y,s],[1,2,3])*revlex([t,u]),
             neglex([x,y,s])*negrevlex([t,u]),
             negdeglex([x,y,s])*negdegrevlex([t,u]),
             negwdeglex([x,y,s],[1,2,3])*negwdegrevlex([t,u],[1,2]))
      @test O == monomial_ordering(R, singular(O))
      @test O == monomial_ordering(R, ordering(singular_poly_ring(R, O)))
   end

   @test_throws ErrorException monomial_ordering(R, Singular.ordering_lp(4))
   @test_throws ErrorException monomial_ordering(R, Singular.ordering_S())
   bad = Singular.ordering_a([1,2,3,4,5,6])*Singular.ordering_rs(5)
   @test_throws ErrorException monomial_ordering(R, bad)

   O1 = degrevlex(gens(R))
   test_opposite_ordering(O1)
   @test monomial_ordering(R, singular(O1)) == O1
   @test length(string(O1)) > 2
   @test string(singular(O1)) == "ordering_dp(5)"

   O2 = lex([x, y])*deglex([s, t, u])
   test_opposite_ordering(O2)
   @test monomial_ordering(R, singular(O2)) == O2
   @test length(string(O2)) > 2
   @test string(singular(O2)) == "ordering_lp(2) * ordering_Dp(3)"

   O3 = wdeglex(gens(R), [2, 3, 5, 7, 3])
   test_opposite_ordering(O3)
   @test monomial_ordering(R, singular(O3)) == O3
   @test length(string(O3)) > 2
   @test string(singular(O3)) == "ordering_Wp([2, 3, 5, 7, 3])"

   O4 = deglex([x, y, t]) * deglex([y, s, u])
   test_opposite_ordering(O4)
   @test monomial_ordering(R, singular(O4)) == O4
   @test length(string(O4)) > 2
   @test string(singular(O4)) == "ordering_M([1 1 0 1 0; 0 -1 0 -1 0; 0 0 0 -1 0; 0 0 1 0 1; 0 0 0 0 -1])"

   K = FreeModule(R, 4)

   O5 = revlex(gens(K))*degrevlex(gens(R))
   @test monomial_ordering(R, singular(O5)) == degrevlex(gens(R))
   @test length(string(O5)) > 2
   @test string(singular(O5)) == "ordering_c() * ordering_dp(5)"

   a = matrix_ordering([x, y], matrix(ZZ, 2, 2, [1 2; 3 4]))
   b = wdeglex([s, t, u], [1, 2, 3])
   O6 = a * lex(gens(K)) * b
   @test monomial_ordering(R, singular(O6)) == a * b
   @test length(string(O6)) > 2
   @test string(singular(O6)) == "ordering_M([1 2; 3 4]) * ordering_C() * ordering_Wp([1, 2, 3])"

   O7 = weight_ordering([-1,2,0,2,0], degrevlex(gens(R)))
   @test monomial_ordering(R, singular(O7)) == O7
   @test length(string(O7)) > 2
   @test string(singular(O7)) == "ordering_a([-1, 2, 0, 2, 0]) * ordering_dp(5)"

   O8 = lex([gen(K,1), gen(K,3), gen(K,4), gen(K,2)]) * degrevlex(gens(R))
   @test_throws ErrorException singular(O8)

   O9 = matrix_ordering([x, y], [1 2; 1 2]; check = false) * lex(gens(R))
   @test monomial_ordering(R, singular(O9)) == O9
   @test singular(O9) isa Singular.sordering
end

@testset "Polynomial Ordering misc bugs" begin
   R, (x, y) = QQ["x", "y"]
   @test degrevlex(gens(R)) != degrevlex(Oscar.reverse(gens(R)))

   R, (x, y, z) = QQ["x", "y", "z"]
   @test degrevlex(gens(R)) != degrevlex(Oscar.reverse(gens(R)))

   a = negwdegrevlex([z, x, y], [4, 5, 6])
   test_opposite_ordering(a)
   @test matrix_ordering([x, y, z], matrix(a)) ==
         matrix_ordering([x, y, z], [-5 -6 -4; 0 -1 0; -1 0 0; 0 0 -1])

   a = matrix_ordering([y, z, x], [4 6 8; 1 0 0; 0 1 0])
   test_opposite_ordering(a)
   @test canonical_matrix(a) == matrix(ZZ, 3, 3, [4 2 3; 0 1 0; 0 0 1])
   @test simplify(a) isa MonomialOrdering
end

@testset "Polynomial Ordering elimination" begin
   R, (x, y, z, w) = QQ["x", "y", "z", "w"]
   @test is_elimination_ordering(lex(R), [x])
   @test is_elimination_ordering(lex(R), [x,y])
   @test is_elimination_ordering(lex(R), [x,y,z])
   @test !is_elimination_ordering(lex(R), [y,z])
   @test is_elimination_ordering(deglex([x,y])*deglex([y,z,w]), [x,y])
end
