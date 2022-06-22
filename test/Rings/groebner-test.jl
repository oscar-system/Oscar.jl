@testset "groebner" begin
    R, (x, y) = PolynomialRing(QQ, ["x", "y"])
    I = ideal(R,[x*y^2 - x, x^3 - 2*y^5])
    @test leading_ideal(I, ordering=degrevlex(gens(R))) == ideal(R,[x*y^2, x^4, y^5])
    @test leading_ideal(I, ordering=lex(gens(R))) == ideal(R,[y^7, x*y^2, x^3])
end

@testset "groebner leading ideal" begin
   R, (t, x, y, z) = PolynomialRing(QQ, ["t", "x", "y", "z"])

   I = ideal(R, [x + y + t + z, x^2 + y^2 + t^3])

   @test leading_ideal(groebner_basis(I, ordering=lex([t, x, y, z])), ordering=lex([t, x, y, z])) == ideal(R, [x^3, t])
   @test leading_ideal(groebner_basis(I, ordering=degrevlex([t, x, y, z])), ordering=degrevlex([t, x, y, z])) == ideal(R, [t, x^3])
   @test leading_ideal(groebner_basis(I, ordering=lex([t, x, y, z])), ordering=lex([t, x, y, z])) == ideal(R, [x^3, t])
   @test leading_ideal(groebner_basis(I, ordering=revlex([t, x, y, z])), ordering=revlex([t, x, y, z])) == ideal(R, [y^2, z])
   @test leading_ideal(groebner_basis(I, ordering=wdeglex([t, x, y, z], [2, 3, 1, 4])), ordering=wdeglex([t, x, y, z], [2, 3, 1, 4])) == ideal(R, [z, t^3])
   @test leading_ideal(groebner_basis(I, ordering=wdegrevlex([t, x, y, z], [2, 1, 1, 1])), ordering=wdegrevlex([t, x, y, z], [2, 1, 1, 1])) == ideal(R, [t, x^3])
end

@testset "groebner orderings" begin
   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

   I = ideal(R, [x + y + z, x^2 + y^2 + z^3])

   @test groebner_basis(I, ordering=lex([x])*lex([y,z])) == groebner_basis(I, ordering=lex([x, y, z]))
   @test groebner_basis(I, ordering=lex([z])*lex([y])*lex([x])) == groebner_basis(I, ordering=revlex([x, y, z]))
   @test groebner_basis(I, ordering=degrevlex([x, y, z])*revlex([y])) == groebner_basis(I, ordering=degrevlex([x, y, z]))
   @test groebner_basis(I, ordering=deglex([z])*deglex([x])*deglex([y])) == groebner_basis(I, ordering=lex([z])*lex([x, y]))
   @test groebner_basis(I, ordering=deglex([x, y, z])) == groebner_basis(I, ordering=wdeglex([x, y, z], [1, 1, 1]))
   M = Oscar.Orderings.MonomialOrdering(R, Oscar.Orderings.ordering([ x, y, z ], matrix(ZZ, [ 1 1 1 ; 0 1 0 ; 1 0 0 ])))
   @test groebner_basis(I, ordering = M) == [ x + y + z, 2*x^2 + 2*x*z + z^3 + z^2 ]
   @test_throws ErrorException groebner_basis(I, ordering = negdeglex([x, y, z]))
   @test groebner_basis(I, ordering = negdeglex([x, y, z]), enforce_global_ordering = false) == [ x + y + z, 2*y^2 + 2*y*z + z^3 + z^2 ]

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
  G = std_basis(I, o)
  @test normal_form(x^5-5, I, o) == -5
  u = negdegrevlex([x])*negdegrevlex([y])
  @test ideal_membership(x^4, I, ordering=u)
end
