@testset "GAP objects" begin
  mktempdir() do path
    @testset "IsObject" begin
      x = GAP.evalstr("Objectify(NewType(NewFamily(\"TestFamily\"),IsObject),rec())")
      filenamex = joinpath(path, "x")
      @test_throws ErrorException Oscar.save(filenamex, x)
    end

    @testset "IsFamily" begin
      x = GAP.Globals.FamilyObj(GapObj(symmetric_group(5)))
      filenamex = joinpath(path, "x")
      @test_throws ErrorException Oscar.save(filenamex, x)
    end

    @testset "IsFreeGroup" begin
      wfilts = [:IsSyllableWordsFamily, :IsLetterWordsFamily]
      for wfilt in [getproperty(GAP.Globals, x) for x in wfilts]
        # full free group of infinite rank
        F = GAP.Globals.FreeGroup(wfilt, GAP.Globals.infinity)
        test_save_load_roundtrip(path, F) do loaded
          @test F == loaded
          @test wfilt(GAP.Globals.ElementsFamily(GAP.Globals.FamilyObj(loaded)))
        end

        # full free group of finite rank
        F = GAP.Globals.FreeGroup(wfilt, 2)
        test_save_load_roundtrip(path, F) do loaded
          @test F == loaded
          @test wfilt(GAP.Globals.ElementsFamily(GAP.Globals.FamilyObj(loaded)))
        end

        # subgroup of full free group
        Fgens = GAP.Globals.GeneratorsOfGroup(F)
        U = GAP.Globals.Subgroup(F, GapObj([Fgens[1]]))
        test_save_load_roundtrip(path, U) do loaded
          @test U == loaded
        end

        # full group and subgroup in two different data files,
        # in the same Julia session
        filenameF = joinpath(path, "F")
        Oscar.save(filenameF, F)
        filenameU = joinpath(path, "U")
        Oscar.save(filenameU, U)
        loadedF = Oscar.load(filenameF)
        loadedU = Oscar.load(filenameU)
        @test GAP.Globals.GeneratorsOfGroup(loadedU)[1] in loadedF

        # full group and subgroup in two different data files,
        # in a new Julia session
        Oscar.reset_global_serializer_state()
        loadedF = Oscar.load(filenameF)
        loadedU = Oscar.load(filenameU)
        @test GAP.Globals.GeneratorsOfGroup(loadedU)[1] in loadedF

        # full group and subgroup together in the same Julia session
        V = GAP.Globals.Subgroup(F, GapObj([Fgens[2]]))
        v = (F, U, V)
        test_save_load_roundtrip(path, v) do loaded
          @test v == loaded
        end

        # full group and subgroup together in a new Julia session
        filenamev = joinpath(path, "v")
        Oscar.save(filenamev, v)
        Oscar.reset_global_serializer_state()
        loaded = Oscar.load(filenamev)
        @test GAP.Globals.GeneratorsOfGroup(loaded[2])[1] in loaded[1]
        @test GAP.Globals.GeneratorsOfGroup(loaded[3])[1] in loaded[1]
      end
    end

    @testset "IsSubgroupFpGroup" begin
      # full f.p. group
      F = GAP.Globals.FreeGroup(2)
      Fgens = GAP.Globals.GeneratorsOfGroup(F)
      G = F/GapObj([x^2 for x in Fgens])
      test_save_load_roundtrip(path, G) do loaded
        @test G == loaded
      end

      # subgroup of full f.p. group
      Ggens = GAP.Globals.GeneratorsOfGroup(G)
      U = GAP.Globals.Subgroup(G, GapObj([Ggens[1]]))
      test_save_load_roundtrip(path, U) do loaded
        @test U == loaded
      end

      # full group and subgroup in two different data files,
      # in the same Julia session
      filenameG = joinpath(path, "G")
      Oscar.save(filenameG, G)
      filenameU = joinpath(path, "U")
      Oscar.save(filenameU, U)
      loadedG = Oscar.load(filenameG)
      loadedU = Oscar.load(filenameU)
      @test GAP.Globals.GeneratorsOfGroup(loadedU)[1] in loadedG

      # full group and subgroup in two different data files,
      # in a new Julia session
      Oscar.reset_global_serializer_state()
      loadedG = Oscar.load(filenameG)
      loadedU = Oscar.load(filenameU)
      @test GAP.Globals.GeneratorsOfGroup(loadedU)[1] in loadedG

      # full group and subgroup together in the same Julia session
      V = GAP.Globals.Subgroup(G, GapObj([Ggens[2]]))
      v = (G, U, V)
      test_save_load_roundtrip(path, v) do loaded
        @test v == loaded
      end

      # full group and subgroup together in a new Julia session
      filenamev = joinpath(path, "v")
      Oscar.save(filenamev, v)
      Oscar.reset_global_serializer_state()
      loaded = Oscar.load(filenamev)
      @test GAP.Globals.GeneratorsOfGroup(loaded[2])[1] in loaded[1]
      @test GAP.Globals.GeneratorsOfGroup(loaded[3])[1] in loaded[1]
    end

    @testset "IsPcGroup" begin
      paras = [(1, 1), (5, 1), (24, 12)]
      for (n, i) in paras
        # full pc group
        G = GAP.Globals.SmallGroup(n, i)
        test_save_load_roundtrip(path, G) do loaded
          @test G == loaded
        end

        # subgroup of full pc group
        U = GAP.Globals.SylowSubgroup(G, 2)
        test_save_load_roundtrip(path, U) do loaded
          @test U == loaded
        end

        # full group and subgroup in two different data files,
        # in the same Julia session
        filenameG = joinpath(path, "G")
        Oscar.save(filenameG, G)
        filenameU = joinpath(path, "U")
        Oscar.save(filenameU, U)
        loadedG = Oscar.load(filenameG)
        loadedU = Oscar.load(filenameU)
        @test loadedG === G
        @test loadedU === U

        # full group and subgroup in two different data files,
        # in a new Julia session
        Oscar.reset_global_serializer_state()
        loadedG = Oscar.load(filenameG)
        loadedU = Oscar.load(filenameU)
        @test loadedG !== G
        @test loadedU !== U
        @test GAP.Globals.IsSubset(loadedG, loadedU)
        @test GAP.Globals.IsomorphismGroups(loadedG, G) != GAP.Globals.fail
        @test GAP.Globals.IsomorphismGroups(loadedU, U) != GAP.Globals.fail

        # full group and subgroup together in the same Julia session
        V = GAP.Globals.SylowSubgroup(G, 3)
        v = (G, U, V)
        test_save_load_roundtrip(path, v) do loaded
          @test v == loaded
        end

        # full group and subgroup together in a new Julia session
        filenamev = joinpath(path, "v")
        Oscar.save(filenamev, v)
        Oscar.reset_global_serializer_state()
        loaded = Oscar.load(filenamev)
        @test GAP.Globals.One(loaded[2]) in loaded[1]
        @test GAP.Globals.One(loaded[3]) in loaded[1]
      end
    end
  end
end
