@testset "Tropical Geometry" begin
    mktempdir() do path
        @testset "Tropical Curves" begin
            @testset "Abstract" begin
                Sigma = graph_from_adjacency_matrix(Undirected,[0 1 1; 1 0 1; 1 1 0]);
                abs_TC = tropical_curve(Sigma, ones(ZZRingElem, n_edges(Sigma)))
                test_save_load_roundtrip(path, abs_TC) do loaded
                    @test graph(loaded) == graph(abs_TC)
                end
            end

            @testset "Embedded" begin
                IM = incidence_matrix([[1,2],[1,3],[1,4]])
                VR = [0 0; 1 0; -1 0; 0 1]
                PC = polyhedral_complex(QQFieldElem, IM, VR)
                TC = tropical_curve(PC)
                test_save_load_roundtrip(path, TC) do loaded
                    loaded_PC = polyhedral_complex(loaded)
                    @test n_rays(PC) == n_rays(loaded_PC)
                    @test number_of_maximal_polyhedra(PC) == number_of_maximal_polyhedra(loaded_PC)
                    @test dim(PC) == dim(loaded_PC)
                end
            end
        end

        @testset "Tropical Hypersurfaces" begin
            T = tropical_semiring(min)
            Txy,(x,y) = T[:x, :y]
            f = x + y^2
            Tf = tropical_hypersurface(f)

            test_save_load_roundtrip(path, inf(T)) do loaded
              @test inf(T) == loaded
            end

            test_save_load_roundtrip(path, Tf) do loaded
                @test f == tropical_polynomial(loaded)
            end
        end
    end
end
