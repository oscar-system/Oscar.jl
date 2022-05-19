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
            L = base_ring(parent(loaded))
            h = hom(K, L, gen(L))
            @test [h(c) for c in coefficients(q)] == collect(coefficients(loaded))
        end

        @testset "MV Polynomial Over Non Simple Extension" begin
            R, t = PolynomialRing(QQ, "t")
            K, a = NumberField([t^2 + 5, t^2 + 7], "a")
            Kxy, (x, y) = PolynomialRing(K, ["x", "y"])
            p = a[1] * x^2 + a[2] * y
            filename = joinpath(path, "polynomial_nse.mv")
            save(p, filename)
            loaded = load(filename)
            loaded_parent = parent(loaded)
            loaded_coeff_ring = coefficient_ring(loaded_parent)
            h = hom(K, loaded_coeff_ring, gens(loaded_coeff_ring))
            @test [h(c) for c in coefficients(p)] == collect(coefficients(loaded))
        end

        @testset "UV Polynomial Over Finite Field Extension" begin
            T, t = PolynomialRing(ResidueRing(ZZ, 2), "t")
            F, a = FiniteField(t^2 + t + 1)
            Fx, x = PolynomialRing(F, "x")
            p = a * x^2 + a^2
            filename = joinpath(path, "polynomial_ffe.uv")
            save(p, filename)
            loaded = load(filename)
            L = base_ring(parent(loaded))
            h = hom(F, L, gen(L))
            @test [h(c) for c in coefficients(p)] == collect(coefficients(loaded))
        end

        @testset "UV Polynomial Over Field Tower Extension" begin
            Qx, x = QQ["x"]
            K, a = NumberField(x^3 - 2, "a")
            Ky, y = K["y"]
            L, b = NumberField(y^2 + 1, "b")
            Lz, z = L["z"]
            p = a^2 * b * z^2 + a * z + b
            filename = joinpath(path, "polynomial_fte.uv")
            save(p, filename)
            loaded = load(filename)
            loaded_base = base_ring(parent(loaded))
            F = base_field(loaded_base)
            h_1 = hom(K, F, gen(F))
            h_2 = hom(L, loaded_base, h_1, gen(loaded_base))
            @test [h_2(c) for c in coefficients(p)] == collect(coefficients(loaded))
        end

        @testset "Non Simple Field Extension over Rel Extension" begin
            Qx, x = QQ["x"]
            K, a = NumberField(x^3 - 2, "a")
            Ky, y = K["y"]
            L, b = NumberField([y^2 - 5 * a, y^2 - 7 * a])
            filename = joinpath(path, "nse_nfrel.elem")
            elem = b[1] + b[2] * a^2
            save(elem, filename)
            elem_loaded = load(filename)
            PF = parent(elem_loaded)
            F = base_field(PF)
            h_1 = hom(K, F, gen(F))
            h_2 = hom(L, PF, h_1, gens(PF))
            @test h_2(elem) == elem_loaded
        end
    end
end
