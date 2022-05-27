@testset "basic_types" begin
    
    mktempdir() do path
        @testset "Testing (de)serialization of $(T)" for T in 
            (
                UInt, UInt8, UInt16, UInt32, UInt64, UInt128,
                Int, Int8, Int16, Int32, Int64, Int128,
                Float16, Float32, Float64,
                BigInt,
                fmpz,
            )
            original = T(1)
            filename = joinpath(path, string(T)*".json")
            save(original, filename)
            loaded = load(filename)
            @test loaded isa T
            @test original == loaded
        end

        @testset "String" begin
            original = "original"
            filename = joinpath(path, "original.json")
            save(original, filename)
            loaded = load(filename)
            @test loaded isa String
            @test loaded == original
        end

        @testset "Symbol" begin
            original = :original
            filename = joinpath(path, "original.json")
            save(original, filename)
            loaded = load(filename)
            @test loaded isa Symbol
            @test loaded == original
        end

        @testset "Singleton types" begin
            original = [ZZ, QQ]
            filename = joinpath(path, "original.json")
            save(original, filename)
            loaded = load(filename)
            @test loaded[1] isa FlintIntegerRing
            @test loaded[1] === ZZ
            @test loaded[2] isa FlintRationalField
            @test loaded[2] === QQ
        end
    end
end
