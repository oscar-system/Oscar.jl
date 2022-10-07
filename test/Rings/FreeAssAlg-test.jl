@testset "FreeAssAlg.orderings" begin
  R, (x, y, z) = FreeAssociativeAlgebra(QQ, ["x", "y", "z"])

  a = (1 + x + y + z)^3

  @test collect(monomials(a, degree_left_lex(R))) ==
    [x^3, x^2*y, x^2*z, x*y*x, x*y^2, x*y*z, x*z*x, x*z*y, x*z^2, y*x^2, y*x*y,
    y*x*z, y^2*x, y^3, y^2*z, y*z*x, y*z*y, y*z^2, z*x^2, z*x*y, z*x*z, z*y*x,
    z*y^2, z*y*z, z^2*x, z^2*y, z^3, x^2, x*y, x*z, y*x, y^2, y*z, z*x, z*y, z^2,
    x, y, z, 1]

  @test collect(monomials(a, degree_right_lex(R))) ==
    [x^3, y*x^2, z*x^2, x*y*x, y^2*x, z*y*x, x*z*x, y*z*x, z^2*x, x^2*y, y*x*y,
    z*x*y, x*y^2, y^3, z*y^2, x*z*y, y*z*y, z^2*y, x^2*z, y*x*z, z*x*z, x*y*z,
    y^2*z, z*y*z, x*z^2, y*z^2, z^3, x^2, y*x, z*x, x*y, y^2, z*y, x*z, y*z, z^2,
    x, y, z, 1]

  @test collect(monomials(a, left_elimination_ordering(R))) ==
    [x^3, x^2*y, x*y*x, y*x^2, x^2*z, x*z*x, z*x^2, x^2, x*y^2, y*x*y, y^2*x,
    x*y*z, x*z*y, y*x*z, y*z*x, z*x*y, z*y*x, x*y, y*x, x*z^2, z*x*z, z^2*x, x*z,
    z*x, x, y^3, y^2*z, y*z*y, z*y^2, y^2, y*z^2, z*y*z, z^2*y, y*z, z*y, y, z^3,
    z^2, z, 1]

  @test collect(monomials(a, right_elimination_ordering(R))) ==
    [z^3, y*z^2, z*y*z, z^2*y, x*z^2, z*x*z, z^2*x, z^2, y^2*z, y*z*y, z*y^2,
    x*y*z, y*x*z, x*z*y, z*x*y, y*z*x, z*y*x, y*z, z*y, x^2*z, x*z*x, z*x^2, x*z,
    z*x, z, y^3, x*y^2, y*x*y, y^2*x, y^2, x^2*y, x*y*x, y*x^2, x*y, y*x, y, x^3,
    x^2, x, 1]

  @test collect(monomials(a, weighted_degree_left_lex(R, [1,2,1]))) ==
    [y^3, x*y^2, y*x*y, y^2*x, y^2*z, y*z*y, z*y^2, x^2*y, x*y*x, x*y*z, x*z*y,
    y*x^2, y*x*z, y^2, y*z*x, y*z^2, z*x*y, z*y*x, z*y*z, z^2*y, x^3, x^2*z, x*y,
    x*z*x, x*z^2, y*x, y*z, z*x^2, z*x*z, z*y, z^2*x, z^3, x^2, x*z, y, z*x, z^2,
    x, z, 1]
end
