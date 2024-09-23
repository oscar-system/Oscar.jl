@testset "Monomials of a given degree" begin
  for R in [QQ[:x, :y, :z][1], graded_polynomial_ring(QQ, [:x, :y, :z])[1]]
    x, y, z = gens(R)
    @test collect(monomials_of_degree(R, 2)) == [x^2, x * y, x * z, y^2, y * z, z^2]
    @test collect(monomials_of_degree(R, 2, 1:3)) == [x^2, x * y, x * z, y^2, y * z, z^2]
    @test collect(monomials_of_degree(R, 3, [1, 3])) == [x^3, x^2 * z, x * z^2, z^3]
    @test collect(monomials_of_degree(R, 3, [3, 1, 1])) == [x^3, x^2 * z, x * z^2, z^3]
    @test collect(monomials_of_degree(R, 3, 1:2)) == [x^3, x^2 * y, x * y^2, y^3]
    @test collect(monomials_of_degree(R, 3, 1:1)) == [x^3]
    @test isempty(monomials_of_degree(R, 2, Int[]))
  end
end
