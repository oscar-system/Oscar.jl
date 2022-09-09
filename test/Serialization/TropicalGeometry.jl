@testset "Tropical Geometry" begin
    mktempdir() do path
        @testset "Tropical Curves" begin
            IM = IncidenceMatrix([[1,2],[1,3],[1,4]])
            @testset "Abstract" begin
                abs_TC = TropicalCurve(IM)
                test_save_load_roundtrip(path, abs_TC) do loaded
                    @test graph(loaded) == graph(abs_TC)
                end
            end

            @testset "Embedded" begin
                VR = [0 0; 1 0; -1 0; 0 1]
                PC = PolyhedralComplex{fmpq}(IM, VR)
                TC = TropicalCurve(PC)
                test_save_load_roundtrip(path, TC) do loaded
                    loaded_PC = underlying_polyhedral_complex(loaded)
                    @test nrays(PC) == nrays(loaded_PC)
                    @test n_maximal_polyhedra(PC) == n_maximal_polyhedra(loaded_PC)
                    @test dim(PC) == dim(loaded_PC)
                end
            end
        end

        @testset "Tropical Hypersurfaces" begin
            T = TropicalSemiring(min)
            Txy,(x,y) = T["x","y"]
            f = x + y^2
            Tf = TropicalHypersurface(f)
            test_save_load_roundtrip(path, Tf) do loaded
                @test test_equality(f, polynomial(loaded))
            end
        end
    end
end
