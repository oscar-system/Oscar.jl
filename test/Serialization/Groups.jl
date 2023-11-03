@testset "Groups" begin
  mktempdir() do path
    @testset "Permutation groups" begin
      G = symmetric_group(5)

      # single element
      x = gen(G, 1)
      test_save_load_roundtrip(path, x) do loaded
        @test x == loaded
      end

      # full symmetric group
      test_save_load_roundtrip(path, G) do loaded
        @test G == loaded
      end

      # subgroup of a symmetric group
      U = sylow_subgroup(G, 2)[1]
      test_save_load_roundtrip(path, U) do loaded
        @test U == loaded
      end
      y = gen(U, 1)
      test_save_load_roundtrip(path, y) do loaded
        @test y == loaded
      end

      # elements and groups together (`Vector{Any}` is not supported)
      v = (x, y, G, U)
      filenamev = joinpath(path, "v")
      save(filenamev, v)
      loadedv = load(filenamev)
      @test v == loadedv

      # three elements from two different groups
      w = (x, x^2, y)
      filenamew = joinpath(path, "w")
      save(filenamew, w)
      loadedw = load(filenamew)
      @test w == loadedw
      @test parent(loadedv[1]) === parent(loadedw[1])

      # simulate loading into a fresh Julia session
      Oscar.reset_global_serializer_state()
      loadedv = load(filenamev)
      @test parent(loadedv[1]) === loadedv[3]

      loadedw = load(filenamew)
      @test parent(loadedw[1]) === parent(loadedw[2])
      @test parent(loadedw[1]) !== parent(loadedw[3])
      @test parent(loadedw[1]) !== G
      @test parent(loadedw[1]) == G
      @test parent(loadedw[3]) !== U
      @test parent(loadedw[3]) == U
      @test loadedw[1] == loadedv[1]
      @test parent(loadedw[1]) === parent(loadedv[1])
    end

    @testset "Free groups" begin
      F = free_group(2)

      # single element
      x = gen(F, 1)
      test_save_load_roundtrip(path, x) do loaded
        @test x == loaded
      end

      # full free group
      test_save_load_roundtrip(path, F) do loaded
        @test F === loaded
      end

      # subgroup of free group
      U = sub(F, [gen(F, 1)])[1]
      test_save_load_roundtrip(path, U) do loaded
        @test U === loaded
      end
      y = gen(U, 1)
      test_save_load_roundtrip(path, y) do loaded
        @test y == loaded
      end

#T missing: serialize a subgroup, deser. must be compatible with full group!

#T missing: two subgroups of a free group (or full group plus subgroup);
#T must be compatible after deser. into new session!

      # elements and groups together (`Vector{Any}` is not supported)
      v = (x, y, F, U)
      filenamev = joinpath(path, "v")
      save(filenamev, v)
      loadedv = load(filenamev)
      @test v == loadedv

      # three elements from two different groups
      w = (x, x^2, y)
      filenamew = joinpath(path, "w")
      save(filenamew, w)
      loadedw = load(filenamew)
      @test w == loadedw
      @test parent(loadedv[1]) === parent(loadedw[1])

      # simulate loading into a fresh Julia session
       Oscar.reset_global_serializer_state()
       loadedv = load(filenamev)
       @test parent(loadedv[1]) === loadedv[3]

       loadedw = load(filenamew)
       @test parent(loadedw[1]) === parent(loadedw[2])
       @test parent(loadedw[1]) !== parent(loadedw[3])
       @test parent(loadedw[1]) !== F
       @test parent(loadedw[3]) !== U
       @test loadedw[1] == loadedv[1]
       @test parent(loadedw[1]) === parent(loadedv[1])
    end

    @testset "Finitely presented groups" begin
      F = free_group(2)
      x1 = gen(F, 1)
      x2 = gen(F, 2)
      G = quo(F, [x1^2, x2^2, comm(x1, x2)])[1]

      # single element
      x = gen(G, 1)
      test_save_load_roundtrip(path, x) do loaded
        @test x == loaded
      end

      # full f.p. group
      test_save_load_roundtrip(path, G) do loaded
        @test G == loaded
      end

      # subgroup of f.p. group
      U = sub(G, [gen(G, 1)])[1]
      test_save_load_roundtrip(path, U) do loaded
        @test gens(U) == gens(loaded)
      end
      y = gen(U, 1)
      test_save_load_roundtrip(path, y) do loaded
        @test y == loaded
      end

       # elements and groups together (`Vector{Any}` is not supported)
       v = (x, y, G, U)
       filenamev = joinpath(path, "v")
       save(filenamev, v)
       loadedv = load(filenamev)
       @test v == loadedv

       # three elements from two different groups
       w = (x, x^2, y)
       filenamew = joinpath(path, "w")
       save(filenamew, w)
       loadedw = load(filenamew)
       @test w == loadedw
       @test parent(loadedv[1]) === parent(loadedw[1])

       # simulate loading into a fresh Julia session
       Oscar.reset_global_serializer_state()
       loadedv = load(filenamev)
       @test parent(loadedv[1]) === loadedv[3]

       loadedw = load(filenamew)
       @test parent(loadedw[1]) === parent(loadedw[2])
       @test parent(loadedw[1]) !== parent(loadedw[3])
       @test parent(loadedw[1]) !== G
       @test parent(loadedw[3]) !== U
       @test loadedw[1] == loadedv[1]
       @test parent(loadedw[1]) === parent(loadedv[1])
    end

    @testset "Pc groups" begin
      paras = [(1, 1), (5, 1), (24, 12)]
      for (n, i) in paras
        G = small_group(n, i)

        # single element
        x = rand(G)
        test_save_load_roundtrip(path, x) do loaded
          @test x == loaded
        end

        # full pc group
        test_save_load_roundtrip(path, G) do loaded
          @test G == loaded
        end

        # subgroup of pc group
        U = sylow_subgroup(G, 2)[1]
        test_save_load_roundtrip(path, U) do loaded
          @test gens(U) == gens(loaded)
        end
        y = rand(U)
        test_save_load_roundtrip(path, y) do loaded
          @test y == loaded
        end

        # elements and groups together (`Vector{Any}` is not supported)
        v = (x, y, G, U)
        filenamev = joinpath(path, "v")
        save(filenamev, v)
        loadedv = load(filenamev)
        @test v == loadedv

        # three elements from two different groups
        w = (x, x^2, y)
        filenamew = joinpath(path, "w")
        save(filenamew, w)
        loadedw = load(filenamew)
        @test w == loadedw
        @test parent(loadedv[1]) === parent(loadedw[1])

        # simulate loading into a fresh Julia session
        Oscar.reset_global_serializer_state()
        loadedv = load(filenamev)
        @test parent(loadedv[1]) === loadedv[3]
        loadedw = load(filenamew)
        @test parent(loadedv[1]) === parent(loadedw[2])
        @test parent(loadedw[1]) === parent(loadedw[2])
        @test parent(loadedw[1]) !== parent(loadedw[3])
        @test parent(loadedw[1]) !== G
        @test parent(loadedw[3]) !== U
        @test loadedw[1] == loadedv[1]
        @test parent(loadedw[1]) === parent(loadedv[1])
      end
    end
  end
end
