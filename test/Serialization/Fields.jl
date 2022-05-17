@testset "Fields" begin
    mktempdir() do path
        @testset "UV Polynomial Over Simple Extension" begin
            R, x = PolynomialRing(QQ, "x")
            p = x^2 + 3//4
            K, a = NumberField(p)
            Ky, y = K["y"]
            q = a * y^2 + y
            filename = joinpath(path, "polynomial_se.uv")
            save(q, filename)
            loaded = load(filename)
            m = is_isomorphic_with_map(K, prime_field(parent(loaded)))
        end

        @testset "MV Polynomial Over Non Simple Extension" begin
            R, t = PolynomialRing(QQ, "t")
            K, a = NumberField([t^2 + 5, t^2 + 7], "a")
            Kxy, (x, y) = K["x", "y"]
            q = a[1] * x^2 + a[2] * y
            filename = joinpath(path, "polynomial_nse.mv")
            save(q, filename)
            loaded = load(filename)
        end

        @testset "UV Polynomial Over Finite Field Extension" begin
            T, t = PolynomialRing(ResidueRing(ZZ, 2), "t")
            F, a = FiniteField(t^2 + t + 1)
            Fx, x = PolynomialRing(F, "x")
            p = a * x^2 + a^2
            filename = joinpath(path, "polynomial_ffe.uv")
            save(p, filename)
            loaded = load(filename)
        end

        @testset "UV Polynomial Over Field Tower Extension" begin
            Qx, x = QQ["x"]
            K, a = NumberField(x^3 - 2, "a", cached=true)
            Ky, y = K["y"]
            L, b = NumberField(y^2 + 1, "b", cached=true)
            Lz, z = L["z"]
            p = a^2 * b * z^2 + a * z + b
            filename = joinpath(path, "polynomial_fte.uv")
            save(p, filename)
            loaded = load(filename)
        end
    end
end
