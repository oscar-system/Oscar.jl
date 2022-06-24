@testset "Fields" begin
    mktempdir() do path
        @testset "Field Embeddings" begin
            Qx, x = QQ["x"]
            K, _ = NumberField(x^2 + 5)
            E = complex_embeddings(K)[1]

            test_save_load_roundtrip(path, E) do loaded
                loaded_K = number_field(loaded)
                g = gen(loaded_K)
                @test contains(loaded(g), E(gen(K)))
                @test contains(E(gen(K)), loaded(g))
            end
        end
    end
end
