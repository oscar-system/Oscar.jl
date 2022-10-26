@testset "groebner" begin
    R, (x, y) = PolynomialRing(QQ, ["x", "y"])
    I = ideal(R,[x*y^2 - x, x^3 - 2*y^5])
    @test leading_ideal(I, ordering=degrevlex(gens(R))) == ideal(R,[x*y^2, x^4, y^5])
    @test leading_ideal(I, ordering=lex(gens(R))) == ideal(R,[y^7, x*y^2, x^3])
    R, (x, y) = PolynomialRing(GF(5), ["x", "y"])
    I = ideal(R, [x])
    gb = f4(I)
    @test normal_form(y, I) == y
end

@testset "groebner leading ideal" begin
   R, (t, x, y, z) = PolynomialRing(QQ, ["t", "x", "y", "z"])

   I = ideal(R, [x + y + t + z, x^2 + y^2 + t^3])

   @test leading_ideal(groebner_basis(I, ordering=lex([t, x, y, z]))) == ideal(R, [x^3, t])
   @test leading_ideal(groebner_basis(I, ordering=degrevlex([t, x, y, z]))) == ideal(R, [t, x^3])
   @test leading_ideal(groebner_basis(I, ordering=lex([t, x, y, z]))) == ideal(R, [x^3, t])
   @test leading_ideal(groebner_basis(I, ordering=revlex([t, x, y, z]))) == ideal(R, [y^2, z])
   @test leading_ideal(groebner_basis(I, ordering=wdeglex([t, x, y, z], [2, 3, 1, 4]))) == ideal(R, [z, t^3])
   @test leading_ideal(groebner_basis(I, ordering=wdegrevlex([t, x, y, z], [2, 1, 1, 1]))) == ideal(R, [t, x^3])
end

@testset "groebner orderings" begin
   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

   I = ideal(R, [x + y + z, x^2 + y^2 + z^3])

   @test groebner_basis(I, ordering=lex([x])*lex([y,z])) == groebner_basis(I, ordering=lex([x, y, z]))
   @test groebner_basis(I, ordering=lex([z])*lex([y])*lex([x])) == groebner_basis(I, ordering=revlex([x, y, z]))
   @test groebner_basis(I, ordering=degrevlex([x, y, z])*revlex([y])) == groebner_basis(I, ordering=degrevlex([x, y, z]))
   @test groebner_basis(I, ordering=deglex([z])*deglex([x])*deglex([y])) == groebner_basis(I, ordering=lex([z])*lex([x, y]))
   @test groebner_basis(I, ordering=deglex([x, y, z])) == groebner_basis(I, ordering=wdeglex([x, y, z], [1, 1, 1]))
   M = matrix_ordering([x, y, z], [1 1 1; 0 1 0; 1 0 0])
   @test gens(groebner_basis(I, ordering = M)) == [ x + y + z, 2*x^2 + 2*x*z + z^3 + z^2 ]
   @test_throws ErrorException groebner_basis(I, ordering = negdeglex([x, y, z]))
   @test gens(standard_basis(I, negdeglex([x, y, z]))) == [ x + y + z, 2*y^2 + 2*y*z + z^3 + z^2 ]

   @test groebner_basis_with_transformation_matrix(I, ordering=lex([x])*lex([y,z])) == groebner_basis_with_transformation_matrix(I, ordering=lex([x, y, z]))
   @test groebner_basis_with_transformation_matrix(I, ordering=lex([z])*lex([y])*lex([x])) == groebner_basis_with_transformation_matrix(I, ordering=revlex([x, y, z]))
   @test groebner_basis_with_transformation_matrix(I, ordering=degrevlex([x, y, z])*revlex([y])) == groebner_basis_with_transformation_matrix(I, ordering=degrevlex([x, y, z]))
   @test groebner_basis_with_transformation_matrix(I, ordering=deglex([z])*deglex([x])*deglex([y])) == groebner_basis_with_transformation_matrix(I, ordering=lex([z])*lex([x, y]))
   @test groebner_basis_with_transformation_matrix(I, ordering=deglex([x, y, z])) == groebner_basis_with_transformation_matrix(I, ordering=wdeglex([x, y, z], [1, 1, 1]))
end

@testset "non-global orderings" begin
  R, (x, y) = QQ["x", "y"]
  I = ideal(R, [x^2*(x-1)-y^2, y^3*(x+y-6)])
  o = negdegrevlex(gens(R))
  G = standard_basis(I, o)
  @test normal_form(x^5-5, I, o) == -5
  u = negdegrevlex([x])*negdegrevlex([y])
  @test ideal_membership(x^4, I, ordering=u)
end

@testset "f4" begin
    R, (x1,x2,x3,x4) = PolynomialRing(GF(next_prime(2^28)), ["x1", "x2", "x3", "x4"], ordering=:degrevlex)
    I = ideal(R,[x1+2*x2+2*x3+2*x4-1,
            x1^2+2*x2^2+2*x3^2+2*x4^2-x1,
            2*x1*x2+2*x2*x3+2*x3*x4-x2,
            x2^2+2*x1*x3+2*x2*x4-x3])
    H = f4(I);
    G = gfp_mpoly[x1 + 2*x2 + 2*x3 + 2*x4 + 268435458
                  x3^2 + 2*x2*x4 + 76695850*x3*x4 + 115043772*x4^2 + 115043768*x2 + 191739613*x3 + 230087535*x4
                  x2*x3 + 268435457*x2*x4 + 230087533*x3*x4 + 76695842*x4^2 + 210913575*x2 + 38347923*x3 + 153391692*x4
                  x2^2 + 2*x2*x4 + 153391692*x3*x4 + 230087538*x4^2 + 230087536*x2 + 115043768*x3 + 191739613*x4
                  x3*x4^2 + 238609298*x4^3 + 14913081*x2*x4 + 56338306*x3*x4 + 129246702*x4^2 + 263464432*x2 + 260150414*x3 + 258493405*x4
                  x2*x4^2 + 89478486*x4^3 + 29826162*x2*x4 + 263464432*x3*x4 + 238609297*x4^2 + 141674270*x2 + 9942054*x3
                  x4^4 + 35851649*x4^3 + 232885084*x2*x4 + 251179133*x3*x4 + 231479137*x4^2 + 191484965*x2 + 85955250*x3 + 167408122*x4]
    @test elements(H) == G
    @test isdefined(I, :gb)
    @test I.gb[degrevlex(gens(base_ring(I)))].O == G
    H = f4(I, eliminate=2);
    G = gfp_mpoly[x3^2*x4 + 73209671*x3*x4^2 + 260301051*x4^3 + 188447115*x3^2 + 167207272*x3*x4 + 120660383*x4^2 + 210590781*x3 + 109814506*x4
                  x3^3 + 156877866*x3*x4^2 + 59264971*x4^3 + 224858274*x3^2 + 183605206*x3*x4 + 130731555*x4^2 + 110395535*x3 + 158620953*x4
                  x4^4 + 167618101*x3*x4^2 + 102789335*x4^3 + 193931678*x3^2 + 156155981*x3*x4 + 60823186*x4^2 + 239040667*x3 + 127377432*x4
                  x3*x4^3 + 99215126*x3*x4^2 + 261328123*x4^3 + 132228634*x3^2 + 93598185*x3*x4 + 85654356*x4^2 + 3613010*x3 + 240673711*x4]
    @test elements(H) == G
    @test I.gb[degrevlex(gens(base_ring(I))[3:end])].O == G
end
