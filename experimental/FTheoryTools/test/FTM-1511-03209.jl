@testset "Test Downloading Artifact and elementary properties" begin
  h = literature_model(arxiv_id = "1511.03209")
  h_resolved = resolve(h, 1)

  @test n_rays(ambient_space(h)) == 104
  @test n_rays(ambient_space(h_resolved)) == 310
end
