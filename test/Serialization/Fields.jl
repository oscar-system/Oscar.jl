@testset "Fields" begin
    mktempdir() do path
        @testset "UV Polynomial Over Simple Extension" begin
            R, x = PolynomialRing(QQ, "x")
            p = x^2 + 3//4
            K, a = NumberField(p)
            Ky, y = K["y"]
            q = a * y^2 + y
            filename = joinpath(path, "polynomial_k.uv")
            save(q, filename)
            loaded = load(filename)
            @test q == loaded
        end

        @testset "MV Polynomial Over Non Simple Extension" begin
            R, t = PolynomialRing(QQ, "t")
            K, a = NumberField([t^2 + 5, t^2 + 7])
            Kxy, (x, y) = K["x", "y"]
            q = a[1] * x^2 + a[2] * y
            filename = joinpath(path, "polynomial_k.mv")
            save(q, filename)
            loaded = load(filename)
            @test q == loaded
        end
    end
end
