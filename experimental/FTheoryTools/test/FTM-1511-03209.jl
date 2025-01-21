@testset "Test Downloading Artifact and elementary properties" begin
  t = literature_model(arxiv_id = "1511.03209")
  t_resolved = resolve(h, 1)

  @test n_rays(ambient_space(t)) == 104
  @test n_rays(ambient_space(t_resolved)) == 310
  @test typeof(get_attribute(t_resolved, :inter_dict)) == Dict{NTuple{4, Int64}, ZZRingElem}
end
