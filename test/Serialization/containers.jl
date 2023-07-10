@testset "Serialization.Containers" begin
    mktempdir() do path
        @testset "Vector{LinearProgram}" begin
            c = cube(3)
            LP0 = linear_program(c, [2,2,-3])
            LP1 = linear_program(c, [2,2,4])
            v = [LP0, LP1]
            test_save_load_roundtrip(path, v) do loaded
              @test length(v) == length(loaded)
              @test feasible_region(loaded[1]) == feasible_region(loaded[2])
              @test feasible_region(loaded[1]) == feasible_region(LP0)
              @test objective_function(loaded[1]) == objective_function(v[1])
              @test objective_function(loaded[2]) == objective_function(v[2])
              @test optimal_value(loaded[1]) == optimal_value(v[1])
              @test optimal_value(loaded[2]) == optimal_value(v[2])
            end
        end

        @testset "Vector{FpFieldElem}" begin
            F = GF(ZZRingElem(77777732222322222232222222223))
            one = F(1)
            minusone = F(-1)
            v = [one, minusone]
            test_save_load_roundtrip(path, v) do loaded
              @test v == loaded
            end
            test_save_load_roundtrip(path, v; parent=F) do loaded
              @test v == loaded
            end
        end
        
        @testset "Vector{fpFieldElem}" begin
            F = GF(7)
            one = F(1)
            minusone = F(-1)
            v = [one, minusone]
            test_save_load_roundtrip(path, v) do loaded
              @test v == loaded
            end
            test_save_load_roundtrip(path, v; parent=F) do loaded
              @test v == loaded
            end
        end

        @testset "Vector{Any}" begin
            c = cube(3)
            LP0 = linear_program(c, [2,2,-3])
            v = [c, LP0]
            test_save_load_roundtrip(path, v) do loaded
              @test length(v) == length(loaded)
              @test loaded[1] isa Polyhedron
              @test loaded[2] isa LinearProgram
              @test loaded isa Vector{Any}
            end
        end
        
        @testset "Vector{Union{Polyhedron, LinearProgram}}" begin
            c = cube(3)
            LP0 = linear_program(c, [2,2,-3])
            v = Vector{Union{Polyhedron, LinearProgram}}([c, LP0])
            test_save_load_roundtrip(path, v) do loaded
              @test length(v) == length(loaded)
              @test loaded[1] isa Polyhedron
              @test loaded[2] isa LinearProgram
              @test loaded isa Vector
            end
        end
        
        @testset "Testing (de)serialization of Vector{$(T)}" for T in 
            (
                UInt, UInt128, UInt16, UInt32, UInt64, UInt8,
                Int, Int128, Int16, Int32, Int64, Int8,
                Float16, Float32, Float64
            )
            original = [T(1), T(2)]
            test_save_load_roundtrip(path, original) do loaded
              @test loaded isa Vector{T}
              @test original == loaded
            end
        end

        @testset "(de)serialization NamedTuple{$(S), $(T)}" for (S, T) in
            (
                (UInt, UInt128), (UInt16, UInt32), (UInt64, UInt8),
                (Int, Int128), (Int16, Int32), (Int64, Int8),
                (Float16, Float32)
            )
            original = (first = S(30), second = T(3.0))
            test_save_load_roundtrip(path, original) do loaded
                @test loaded isa NamedTuple{(:first, :second), Tuple{S, T}}
                @test original == loaded
            end
        end

        @testset "Test for backwards compatibility" begin
            loaded_container = load(joinpath(@__DIR__, "old-containers.json"))
            @test loaded_container == (r = QQFieldElem(1, 2), m = QQFieldElem[1//2 1; 0 1], t = (1, 2, 3))             
        end
    end
end
