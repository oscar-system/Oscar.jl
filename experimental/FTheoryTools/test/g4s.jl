##############################################################
# 1: Compute well-quantized fluxes and make consistency checks
##############################################################

qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
g4_base = ambient_space_models_of_g4_fluxes(qsm_model, check = false)
res = well_quantized_ambient_space_models_of_g4_fluxes(qsm_model, check = false)
M=res[1]
g4_classes = [sum(M[i,j]*g4_base[i] for i in 1:length(g4_base)) for j in 1:size(M,2)]
g4_list = [g4_flux(qsm_model, cl, check = false) for cl in g4_classes]
g4_dummy = g4_flux(qsm_model, sum(rand(-100:100)*g for g in g4_classes), check = false)

@testset "Execute well-quantized tests" begin
  @test length(g4_base) == 37
  @test nrows(res[1]) == 37
  @test ncols(res[1]) == 37
  @test nrows(res[2]) == 37
  @test ncols(res[2]) == 0
  @test all(k -> passes_elementary_quantization_checks(k), g4_list) == true
  @test passes_elementary_quantization_checks(g4_dummy)
end



##############################################################
# 2: Compute well-quantized and vertical fluxes
##############################################################

res2 = well_quantized_and_vertical_ambient_space_models_of_g4_fluxes(qsm_model, check = false)
M2=res2[1]
g4_classes = [sum(M2[i,j]*g4_base[i] for i in 1:length(g4_base)) for j in 1:size(M2,2)]
g4_list = [g4_flux(qsm_model, cl, check = false) for cl in g4_classes]
g4_dummy = g4_flux(qsm_model, sum(rand(-100:100)*g for g in g4_classes), check = false)

@testset "Execute well-quantized tests" begin
  @test length(g4_base) == 37
  @test nrows(res2[1]) == 37
  @test ncols(res2[1]) == 25
  @test nrows(res2[2]) == 37
  @test ncols(res2[2]) == 0
  @test all(k -> passes_elementary_quantization_checks(k), g4_list) == true
  @test all(k -> passes_verticality_checks(k), g4_list) == true
  @test passes_elementary_quantization_checks(g4_dummy)
  @test passes_verticality_checks(g4_dummy)
end
