@testset "QuadForm" begin
    
    mktempdir() do path


        @testset "ZLat" begin
            L = Zlattice(ZZ[1 2;],gram=ZZ[0 1;1 0])
            V = ambient_space(L)

            test_save_load_roundtrip(path, L) do loaded
              @test L == loaded
            end
            test_save_load_roundtrip(path, V) do loaded
              @test V == loaded
            end
        end

        @testset "QuadSpace" begin
            F,a = CyclotomicField(5)
            V = quadratic_space(F, F[a;])

            test_save_load_roundtrip(path, V) do loaded
              F1 = base_ring(loaded)
              iso = hom(F1,F,gen(F))
              @test gram_matrix(V) == map_entries(iso,gram_matrix(loaded))
            end
        end
    end
end
