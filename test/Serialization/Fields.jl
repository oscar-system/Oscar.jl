@testset "Fields" begin
    mktempdir() do path
        @testset "Field Embeddings" begin
            Qx, x = QQ["x"]
            K, _ = NumberField(x^2 + 5)
            E = complex_embeddings(K)[1]

            test_save_load_roundtrip(path, E) do loaded
                loaded_K = number_field(loaded)
                g = gen(loaded_K)
                @test overlaps(loaded(g), E(gen(K)))
                @test overlaps(E(gen(K)), loaded(g))
            end
        end

        @testset "Non Simple Field Embeddings" begin
            Qx, x = QQ["x"]
            K, _ = NumberField([x^2 + 5, x^2 + 7])
            E = complex_embeddings(K)[1]

            test_save_load_roundtrip(path, E) do loaded
                loaded_K = number_field(loaded)
                loaded_gens = map(loaded, gens(loaded_K))
                E_gens = map(E, gens(K))
                compare_vector_1 = collect(zip(E_gens, loaded_gens))
                compare_vector_2 = collect(zip(loaded_gens, E_gens))
                check_overlaps = y -> isempty(filter(x -> !overlaps(x...), y))
                @test check_overlaps(compare_vector_1)
                @test check_overlaps(compare_vector_2)
            end
        end
    end
end
