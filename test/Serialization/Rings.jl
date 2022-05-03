@testset "Rings" begin
    mktempdir() do path
        @testset "Polynomials Over QQ" begin
            R, (x, y) = PolynomialRing(QQ, ["x", "y"])
            p = x^2 + 3//4*x*y + y + 1//2
            filename = joinpath(path, "polynomial_q.mv")
            save(p, filename)
            loaded = load(filename)
            @test p == loaded
        end

        @testset "Polynomials Over ZZ" begin
            R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
            p = x^2 + 3*x*y + y + 12
            filename = joinpath(path, "polynomial_z.mv")
            save(p, filename)
            loaded = load(filename)
            @test p == loaded
        end

        @testset "Polynomials Over ZZ Residue Ring" begin
            R = ResidueRing(ZZ, 6)
            MR, (x, y) = PolynomialRing(R, ["x", "y"])
            p = 3 * x^2 - 2 * y + 5
            filename = joinpath(path, "polynomial_rr.mv")
            save(p, filename)
            loaded = load(filename)
            @test p == loaded
        end
    end
end
