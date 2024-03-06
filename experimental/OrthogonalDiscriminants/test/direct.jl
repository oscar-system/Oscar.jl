@testset "compute via the Atlas of Group Representations" begin
#TODO: How can we test functions that need web access?
#  # Check equality of positive results in small cases.
#  l = all_od_infos(comment_matches => "AGR", dim => 1:8);
#  for entry in l
#    chi = Oscar.OrthogonalDiscriminants.character_of_entry(entry)
#    comp = Oscar.OrthogonalDiscriminants.od_from_atlas_group(chi)
#    @test comp == (true, entry[:valuestring])
#  end
#
#  # Check negative results.
#  t = character_table("L2(17)")
#  chi = t[10]
#  @test Oscar.OrthogonalDiscriminants.od_from_atlas_group(chi) == (false, "")
#  t = character_table("L2(13)", 3)
#  chi = t[4]
#  @test Oscar.OrthogonalDiscriminants.od_from_atlas_group(chi) == (false, "")
end
