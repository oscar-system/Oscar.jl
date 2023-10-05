@testset "GAP objects" begin
  mktempdir() do path
    @testset "IsObject" begin
      x = GAP.evalstr("Objectify(NewType(NewFamily(\"TestFamily\"),IsObject),rec())")
      filenamex = joinpath(path, "x")
      @test_throws ErrorException Oscar.save(filenamex, x)
    end

    @testset "IsFamily" begin
      x = GAP.Globals.FamilyObj(symmetric_group(5).X)
      filenamex = joinpath(path, "x")
      @test_throws ErrorException Oscar.save(filenamex, x)
    end

    @testset "IsFreeGroup" begin
      F = GAP.Globals.FreeGroup(2)
      test_save_load_roundtrip(path, F) do loaded
        @test F == loaded
      end
      Fgens = GAP.Globals.GeneratorsOfGroup(F)
      U = GAP.Globals.Subgroup(F, GAP.GapObj([Fgens[1]]))
      test_save_load_roundtrip(path, U) do loaded
        @test U == loaded
      end
#TODO: different internal representations
    end

    @testset "IsSubgroupFpGroup" begin
      F = GAP.Globals.FreeGroup(2)
      Fgens = GAP.Globals.GeneratorsOfGroup(F)
      G = F/GAP.GapObj([x^2 for x in Fgens])
      test_save_load_roundtrip(path, G) do loaded
        @test G == loaded
      end
      Ggens = GAP.Globals.GeneratorsOfGroup(G)
      U = GAP.Globals.Subgroup(G, GAP.GapObj([Ggens[1]]))
      test_save_load_roundtrip(path, U) do loaded
        @test U == loaded
      end
    end
  end
end
