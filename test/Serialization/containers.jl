@testset "Serialization.Containers" begin
  mktempdir() do path
    @testset "Empty Containers" begin
      v = Int[]
      test_save_load_roundtrip(path, v) do loaded
        v == loaded
      end

      t = Tuple{Vector{Int}, Vector{Int}}([v, [0]])
      test_save_load_roundtrip(path, t) do loaded
        @test t == loaded
      end

      nt = (a = v, b = t)
      test_save_load_roundtrip(path, nt) do loaded
        @test nt == loaded
      end
    end
  
    @testset "ids in containers" begin
      R, x = QQ[:x]
      test_save_load_roundtrip(path, (x^2, x + 1, R)) do loaded
        @test loaded[3] == R
        @test parent(loaded[1]) == parent(loaded[2]) == loaded[3]
      end
    end
    
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
      test_save_load_roundtrip(path, v; params=F) do loaded
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
      test_save_load_roundtrip(path, v; params=F) do loaded
        @test v == loaded
      end
    end

    @testset "Tuple" begin
      c = cube(3)
      LP0 = linear_program(c, [2,2,-3])
      v = (c, LP0)
      test_save_load_roundtrip(path, v) do loaded
        @test length(v) == length(loaded)
        @test loaded[1] isa Polyhedron
        @test loaded[2] isa LinearProgram
        @test loaded isa Tuple
      end
    end

    # Does it make sense to save such types?
    # I feel that these types should be discouraged?
    # @testset "Vector{Union{Polyhedron, LinearProgram}}" begin
    #   c = cube(3)
    #   LP0 = linear_program(c, [2,2,-3])
    #   v = Vector{Union{Polyhedron, LinearProgram}}([c, LP0])
    #   test_save_load_roundtrip(path, v) do loaded
    #     @test length(v) == length(loaded)
    #     @test loaded[1] isa Polyhedron
    #     @test loaded[2] isa LinearProgram
    #     @test loaded isa Vector
    #   end
    # end
    
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
