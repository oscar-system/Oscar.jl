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
      v = LinearProgram{QQFieldElem}[LP0, LP1]
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
      F = FpField(ZZRingElem(77777732222322222232222222223))
      one = F(1)
      minusone = F(-1)
      v = [one, minusone]
      test_save_load_roundtrip(path, v) do loaded
        @test v == loaded
      end
      test_save_load_roundtrip(path, v; params=F) do loaded
        @test v == loaded
      end

      # tests loading into other types
      filename = joinpath(path, "original.json")
      save(filename, v;)
      loaded = load(filename; params=GF(ZZRingElem(77777732222322222232222222223)),
                    type=Vector{FqFieldElem})
      @test loaded isa Vector{FqFieldElem}
    end
    
    @testset "Vector{fpFieldElem}" begin
      F = fpField(UInt(7))
      one = F(1)
      minusone = F(-1)
      v = [one, minusone]
      test_save_load_roundtrip(path, v) do loaded
        @test v == loaded
      end
      test_save_load_roundtrip(path, v; params=F) do loaded
        @test v == loaded
      end
      # tests loading into other types
      filename = joinpath(path, "original.json")
      save(filename, v;)
      loaded = load(filename; params=GF(ZZRingElem(77777732222322222232222222223)),
                    type=Vector{FqFieldElem})
      @test loaded isa Vector{FqFieldElem}
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

    @testset "Testing (de)serialization of Multidimensional Arrays" begin
      original = [1; 2;; 3; 4;; 5; 6;;;
                  7; 8;; 9; 10;; 11; 12]
      test_save_load_roundtrip(path, original) do loaded
        @test original == loaded
      end

      Qxy,(x, y) = QQ[:x, :y]
      original = Array{MPolyRingElem, 4}(
        reshape(reduce(hcat, [x .* Qxy.(arr), Qxy.(arr)]), 2, 3, 2, 2)
      )
      test_save_load_roundtrip(path, original) do loaded
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

    @testset "(de)serialization Dict{$S, Any}" for (S, keys) in
      (
        (String, ["a", "b"]), (Int, [1, 2]), (Symbol, [:a, :b])
      )
      Qx, x = QQ[:x]
      p = x^2 + 1
      original = Dict{S, Any}(keys[1] => cube(2), keys[2] => p)
      test_save_load_roundtrip(path, original) do loaded
        @test original == loaded
      end
    end

    @testset "(de)serialization Dict{Symbol, T}" begin
      Qx, x = QQ[:x]
      for (T, values) in ((Int, [1, 2]), (QQPolyRingElem, [x^2, x - 1]))
        original = Dict{Symbol, T}(:a => values[1], :b => values[2])
        test_save_load_roundtrip(path, original) do loaded
          @test original == loaded
        end
      end

      original = Dict{Symbol, Int}()
      test_save_load_roundtrip(path, original) do loaded
        @test original == loaded
      end
    end

    @testset "Testing (de)serialization of Set" begin
      original = Set{Set{Int}}([Set{Int}([1, 2])])
      test_save_load_roundtrip(path, original) do loaded
        @test original == loaded
      end

      Qx, x = QQ[:x]
      p = x^2 + 1
      q = x
      original = Set{PolyRingElem}([p, q])
      test_save_load_roundtrip(path, original) do loaded
        @test original == loaded
      end

      test_save_load_roundtrip(path, original; params=Qx) do loaded
        @test original == loaded
      end
    end
  end
end
