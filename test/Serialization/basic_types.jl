@testset "basic_types" begin
    
    mktempdir() do path
        @testset "Testing (de)serialization of $T" for T in 
            [
                UInt, UInt128, UInt16, UInt32, UInt64, UInt8,
                Int, Int128, Int16, Int32, Int64, Int8,
                Float16, Float32, Float64
            ]
            original = T(1)
            filename = joinpath(path, string(T)*".json")
            save(original, filename)
            loaded = load(filename)
            @test loaded isa T
            @test original == loaded
        end
        
        @testset "Testing (de)serialization of Vector{$T}" for T in 
            [
                UInt, UInt128, UInt16, UInt32, UInt64, UInt8,
                Int, Int128, Int16, Int32, Int64, Int8,
                Float16, Float32, Float64
            ]
            original = [T(1), T(2)]
            filename = joinpath(path, string(T)*"_vec.json")
            save(original, filename)
            loaded = load(filename)
            @test loaded isa Vector{T}
            @test original == loaded
        end
    end
end
