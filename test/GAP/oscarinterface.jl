@testset "OscarInterface" begin

  # Test the OscarInterface package
  GAP_assertion_level = GAP.Globals.AssertionLevel()
  @test GAP.Globals.TestDirectory(GAP.Globals.DirectoriesPackageLibrary(GAP.Obj("OscarInterface"), GAP.Obj("tst")))
  GAP.Globals.SetAssertionLevel(GAP_assertion_level)

end
