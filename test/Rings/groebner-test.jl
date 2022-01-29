@testset "groebner" begin
    R, (x, y) = PolynomialRing(QQ, ["x", "y"])
    I = ideal(R,[x*y^2 - x, x^3 - 2*y^5])
    @test leading_ideal(I) == ideal(R,[x*y^2, x^4, y^5])
    @test leading_ideal(I, :lex) == ideal(R,[y^7, x*y^2, x^3])
end

@testset "groebner leading ideal" begin
   R, (t, x, y, z) = PolynomialRing(QQ, ["t", "x", "y", "z"])

   I = ideal(R, [x + y + t + z, x^2 + y^2 + t^3])

   @test leading_ideal(groebner_basis(I, lex([t, x, y, z])), lex([t, x, y, z])) == ideal(R, [x^3, t])
   @test leading_ideal(groebner_basis(I, degrevlex([t, x, y, z])), degrevlex([t, x, y, z])) == ideal(R, [t, x^3])
   @test leading_ideal(groebner_basis(I, lex([t, x, y, z])), lex([t, x, y, z])) == ideal(R, [x^3, t])
   @test leading_ideal(groebner_basis(I, revlex([t, x, y, z])), revlex([t, x, y, z])) == ideal(R, [y^2, z])
   @test leading_ideal(groebner_basis(I, wdeglex([t, x, y, z], [2, 3, 1, 4])), wdeglex([t, x, y, z], [2, 3, 1, 4])) == ideal(R, [z, t^3])
   @test leading_ideal(groebner_basis(I, wdegrevlex([t, x, y, z], [2, 1, 1, 1])), wdegrevlex([t, x, y, z], [2, 1, 1, 1])) == ideal(R, [t, x^3])
end

@testset "groebner orderings" begin
   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

   I = ideal(R, [x + y + z, x^2 + y^2 + z^3])

   @test groebner_basis(I, lex([x])*lex([y,z])) == groebner_basis(I, lex([x, y, z]))
   @test groebner_basis(I, lex([z])*lex([y])*lex([x])) == groebner_basis(I, revlex([x, y, z]))
   @test groebner_basis(I, degrevlex([x, y, z])*revlex([y])) == groebner_basis(I, degrevlex([x, y, z]))
   @test groebner_basis(I, deglex([z])*deglex([x])*deglex([y])) == groebner_basis(I, lex([z])*lex([x, y]))
   @test groebner_basis(I, deglex([x, y, z])) == groebner_basis(I, wdeglex([x, y, z], [1, 1, 1]))

   @test groebner_basis_with_transformation_matrix(I, lex([x])*lex([y,z])) == groebner_basis_with_transformation_matrix(I, lex([x, y, z]))
   @test groebner_basis_with_transformation_matrix(I, lex([z])*lex([y])*lex([x])) == groebner_basis_with_transformation_matrix(I, revlex([x, y, z]))
   @test groebner_basis_with_transformation_matrix(I, degrevlex([x, y, z])*revlex([y])) == groebner_basis_with_transformation_matrix(I, degrevlex([x, y, z]))
   @test groebner_basis_with_transformation_matrix(I, deglex([z])*deglex([x])*deglex([y])) == groebner_basis_with_transformation_matrix(I, lex([z])*lex([x, y]))
   @test groebner_basis_with_transformation_matrix(I, deglex([x, y, z])) == groebner_basis_with_transformation_matrix(I, wdeglex([x, y, z], [1, 1, 1]))
end