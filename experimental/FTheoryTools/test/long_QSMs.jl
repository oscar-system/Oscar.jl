#############################################################
# Extra long QSM model tests
#############################################################

using Random
our_rng = Random.Xoshiro(1234)

@testset "Advanced intersection theory and QSM-fluxes" begin
  for k in 1:5000
    qsm_model = try literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => k), rng = our_rng) catch e continue end
    h22_converter_dict = converter_dict_h22_ambient(qsm_model, completeness_check = false)
    coh_ring = cohomology_ring(ambient_space(qsm_model), completeness_check = false)
    coh_ring_gens = gens(coh_ring)
    for (key, value) in h22_converter_dict
      obj1 = coh_ring_gens[key[1]] * coh_ring_gens[key[2]]
      obj2 = sum(value[k][1] * coh_ring_gens[value[k][2][1]] * coh_ring_gens[value[k][2][2]] for k in 1:length(value))
      @test obj1 == obj2
    end
    qsm_g4_flux = qsm_flux(qsm_model)
    h22_basis = gens_of_h22_hypersurface_indices(qsm_model, completeness_check = false)
    flux_poly_str = string(polynomial(cohomology_class(qsm_g4_flux)))
    ring = base_ring(parent(polynomial(cohomology_class(qsm_g4_flux))))
    flux_poly = Oscar.eval_poly(flux_poly_str, ring)
    coeffs = collect(coefficients(flux_poly))
    raw_exps = collect(exponents(flux_poly))
    massaged_exps = [length(pos) == 1 ? (pos[1], pos[1]) : Tuple(sort(pos)) for pos in (findall(!=(0), mon) for mon in raw_exps)]
    flux_vector = fill(QQ(0), length(h22_basis))
    for (i, exp_pair) in enumerate(massaged_exps)
      idx = findfirst(==(exp_pair), h22_basis)
      flux_vector[idx] = coeffs[i]
    end
    flux_vector = transpose(matrix(QQ, [flux_vector]))
    fg = special_flux_family(qsm_model; not_breaking = true, completeness_check = false, algorithm = "special")
    @test ncols(matrix_integral(fg)) == 1
    @test nrows(matrix_integral(fg)) == nrows(matrix_rational(fg))
    @test unique(offset(fg)) == [0]
    M1 = matrix_integral(fg)
    M2 = matrix_rational(fg)
    large_M = hcat(M1, M2)
    @test rank(large_M) == minimum(size(large_M))
    solution = solve(large_M, flux_vector, side = :right)
    @test is_integer(solution[1])
    reconstructed_flux = flux_instance(fg, matrix(ZZ, [[solution[1]]]), solution[2:end,:])
    @test cohomology_class(qsm_g4_flux) == cohomology_class(reconstructed_flux)
    coho_R = cohomology_ring(ambient_space(qsm_model), completeness_check = false)
    gs = [gg.f for gg in gens(coho_R)]
    known_intersections = qsm_model.__attrs[:inter_dict]
    kbar_poly = polynomial(cohomology_class(anticanonical_bundle(ambient_space(qsm_model)))).f
    non_zero_entries = collect(filter(x -> x[2] != 0, known_intersections))
    sampled_dict = Dict()
    while length(sampled_dict) < min(100, length(non_zero_entries))
      i = rand(our_rng, 1:length(non_zero_entries))
      sampled_dict[non_zero_entries[i][1]] = non_zero_entries[i][2]
    end
    for (k,v) in sampled_dict
      desired_class = CohomologyClass(ambient_space(qsm_model), coho_R(gs[k[1]] * gs[k[2]] * gs[k[3]] * gs[k[4]] * kbar_poly), true)
      @test v == integrate(desired_class, completeness_check = false)
    end
  end
end
