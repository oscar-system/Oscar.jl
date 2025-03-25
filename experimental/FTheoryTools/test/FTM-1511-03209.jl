@testset "Test Downloading Artifact and elementary properties" begin
  t = literature_model(arxiv_id = "1511.03209")
  fully_resolved_big_model = resolve(t, 1)
  f1 = special_flux_family(fully_resolved_big_model, check = false)
  g1 = random_flux_instance(f1, check = false)
  f2 = special_flux_family(fully_resolved_big_model, not_breaking = true, check = false)
  g2 = random_flux_instance(f2, check = false)
  @test n_rays(ambient_space(t)) == 104
  @test n_rays(ambient_space(fully_resolved_big_model)) == 313
  @test typeof(get_attribute(fully_resolved_big_model, :inter_dict)) == Dict{NTuple{4, Int64}, ZZRingElem}
  @test length(chosen_g4_flux_basis(fully_resolved_big_model)) == 629
  @test is_well_quantized(g1) == true
  @test breaks_non_abelian_gauge_group(g2) == false
  @test size(matrix_integral(f1)) == (629, 224)
  @test size(matrix_rational(f1)) == (629, 127)
  @test size(matrix_integral(f2)) == (629, 1)
  @test size(matrix_rational(f2)) == (629, 127)
end
