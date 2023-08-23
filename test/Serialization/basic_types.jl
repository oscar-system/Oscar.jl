function test_save_load_roundtrip(func, path, original::T; params=nothing) where T
  # save and load from a file
  filename = joinpath(path, "original.json")
  save(filename, original)
  loaded = load(filename; params=params)
  if T <: Vector
    @test loaded isa Vector
  else
    @test loaded isa T
  end
  func(loaded)

  # save and load from an IO buffer
  io = IOBuffer()
  save(io, original)
  seekstart(io)
  loaded = load(io; params=params)

  if T <: Vector
    @test loaded isa Vector
  else
    @test loaded isa T
  end
  func(loaded)

  # save and load from an IO buffer, with prescribed type
  io = IOBuffer()
  save(io, original)
  seekstart(io)
  loaded = load(io, type=T, params=params)
  if T <: Vector
    @test loaded isa Vector
  else
    @test loaded isa T
  end
  func(loaded)

  # test loading on a empty state
  save(filename, original)
  reset_global_serializer_state()
  loaded = load(filename; params=params)
end

@testset "basic_types" begin
    
    mktempdir() do path
        @testset "Testing (de)serialization of $(T)" for T in 
            (
                UInt, UInt8, UInt16, UInt32, UInt64, UInt128,
                Int, Int8, Int16, Int32, Int64, Int128,
                Float16, Float32, Float64,
                BigInt,
                ZZRingElem,
                QQFieldElem,
                residue_ring(ZZ, 6),
                Nemo.fpField(UInt(7)),
                Nemo.FpField(ZZRingElem(7)),
                #PadicField(7, 30),
                #TropicalSemiring()
            )
            original = T(1)
            test_save_load_roundtrip(path, original) do loaded
              @test loaded == original
            end
        end
        
        @testset "String" begin
            original = "original"
            test_save_load_roundtrip(path, original) do loaded
              @test loaded == original
            end
        end

        @testset "Symbol" begin
            original = :original
            test_save_load_roundtrip(path, original) do loaded
              @test loaded == original
            end
        end

        @testset "Singleton types" begin
            test_save_load_roundtrip(path, ZZ) do loaded
              @test loaded === ZZ
            end
            test_save_load_roundtrip(path, QQ) do loaded
              @test loaded === QQ
            end
        end
    end
end
