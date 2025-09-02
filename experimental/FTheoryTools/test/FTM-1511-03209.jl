#############################################################
# Tests for the model from 1511.03209
#############################################################

@testset "Test Downloading Artifact and elementary properties" begin
  t = literature_model(arxiv_id = "1511.03209")
  fully_resolved_big_model = resolve(t, 1)
  f1 = special_flux_family(fully_resolved_big_model, completeness_check = false)
  g1 = random_flux_instance(f1, completeness_check = false, consistency_check = false)
  f2 = special_flux_family(fully_resolved_big_model, not_breaking = true, completeness_check = false)
  g2 = random_flux_instance(f2, completeness_check = false, consistency_check = false)
  @test n_rays(ambient_space(t)) == 104
  @test n_rays(ambient_space(fully_resolved_big_model)) == 313
  @test typeof(get_attribute(fully_resolved_big_model, :inter_dict)) == Dict{NTuple{4, Int64}, ZZRingElem}
  @test length(chosen_g4_flux_gens(fully_resolved_big_model)) == 629
  @test is_well_quantized(g1) == true
  @test breaks_non_abelian_gauge_group(g2) == false
  @test size(matrix_integral(f1)) == (629, 224)
  @test size(matrix_rational(f1)) == (629, 127)
  @test size(matrix_integral(f2)) == (629, 1)
  @test size(matrix_rational(f2)) == (629, 127)
  @test has_attribute(fully_resolved_big_model, :exceptional_classes)
  @test has_attribute(fully_resolved_big_model, :exceptional_divisor_indices)
  @test (104 in exceptional_divisor_indices(fully_resolved_big_model)) == false
  @test 105 in exceptional_divisor_indices(fully_resolved_big_model)
  @test length(exceptional_divisor_indices(fully_resolved_big_model)) == 206
  @test length(fully_resolved_big_model.__attrs) == 48
  @test length(fully_resolved_big_model.__attrs[:inter_dict]) == 14154797
  @test maximum(values(fully_resolved_big_model.__attrs[:inter_dict])) == 407568
  @test fully_resolved_big_model.__attrs[:inter_dict][(103,103,103,103)] == 407568
  @test fully_resolved_big_model.__attrs[:inter_dict][(104,104,104,104)] == -6654
  @test length(fully_resolved_big_model.__attrs[:s_inter_dict]) == 66
  @test paper_buzzwords(t) == [ "Tate", "Most flux vacua"]
  @test paper_buzzwords(fully_resolved_big_model) == [ "Tate", "Most flux vacua"]
  @test paper_authors(fully_resolved_big_model) == ["Washington Taylor", "Yi-Nan Wang"]
end
