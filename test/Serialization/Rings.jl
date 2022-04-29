@testset "Rings" begin

    mktempdir() do path
        @testset "Polynomials" begin
            R, (x, y) = PolynomialRing(QQ, ["x", "y"])
            p = x^2 + 3//4*x*y + y + 1//2
            filename = joinpath(path, "polynomial.mv")
            save(p)
            loaded = load(filename)
            @test p == loaded
        end
end
