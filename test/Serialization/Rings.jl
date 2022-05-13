@testset "Rings" begin
    mktempdir() do path
        @testset "MV Polynomial Over QQ" begin
            R, (x, y) = PolynomialRing(QQ, ["x", "y"])
            p = x^2 + 3//4*x*y + y + 1//2
            filename = joinpath(path, "polynomial_q.mv")
            save(p, filename)
            loaded = load(filename)
            @test p == loaded

        end

        @testset "MV One Variable Polynomial Over ZZ" begin
            R, x = PolynomialRing(ZZ, ["x"])
            p = x[1]^2 + 3*x[1]+ 12
            filename = joinpath(path, "polynomial_z.mv")
            save(p, filename)
            loaded = load(filename)
            @test p == loaded
            @test typeof(loaded) == fmpz_mpoly

        end

        @testset "UV Polynomial Over ZZ Residue Ring" begin
            R = ResidueRing(ZZ, 6)
            MR, y = PolynomialRing(R, "y")
            p = 3 * y^2 - 2 * y + 5
            filename = joinpath(path, "polynomial_rr.uv")
            save(p, filename)
            loaded = load(filename)
            @test p == loaded
        end

        @testset "MV Polynomial Ideal" begin
            R = ResidueRing(ZZ, 6)
            MR, (x, y) = PolynomialRing(R, ["x", "y"])
            p = 3 * x^2 - 2 * y + 5
            filename = joinpath(path, "polynomial_rr.mv")
            save(p, filename)
            loaded = load(filename)
            @test p == loaded

            q = y^2 - x
            i = ideal(MR, [p, q])
            filename = joinpath(path, "ideal_q.i")
            save(i, filename)
            loaded_i = load(filename)
            @test i == loaded_i
        end
    end
end
