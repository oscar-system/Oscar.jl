using Random
our_rng = Random.Xoshiro(1234)

#@testset "Blowups of global Tate models" begin
#id_i5_s = ideal([tate_polynomial(t_i5_s)]);
#tas = ambient_space(t_i5_s);
#irr_i5_s = irrelevant_ideal(tas);
#lin_i5_s = ideal_of_linear_relations(tas);
#id_fin,  = _blowup_global_sequence(id_i5_s, [[8, 9, 6], [2, 3, 1], [3, 4], [2, 4]], irr_i5_s, sri_i5_s, lin_i5_s)
#@test string(gens(id_fin)[end]) == "-b_4_1*b_2_1*a1p*z - b_4_1*b_2_2 - b_4_1*b_2_3*b_1_3^2*a3p*z^3 + b_4_2*b_3_2*b_2_1^2*b_1_1 + b_4_2*b_3_2*b_2_1^2*b_1_3*a2p*z^2 + b_4_2*b_3_2*b_2_1*b_2_3*b_1_3^3*a4p*z^4 + b_4_2*b_3_2*b_2_3^2*b_1_3^5*a6p*z^6"
#end

#@testset "Fibers" begin
#  inters = analyze_fibers(t_i5_s, [[7, 8, 6], [2, 3, 1], [3, 4], [2, 4]])
#  @test string(inters[1][2][1][2][1]) == "ideal(e_2*b_2_1 - b_1_1, b_3_1*a1p*z - b_3_2*b_1_1^2, b_4_1*a1p*z - b_4_2*b_3_2*b_2_1*b_1_1, b_4_1*b_1_1 - b_4_2*b_3_1*b_2_1, b_4_1*e_2 - b_4_2*b_3_1, e_4*b_4_2 - e_2, e_4*b_4_1 - b_3_1, y, x, v, b_1_3, b_1_2, e_1, b_2_3, b_2_2, e_3)"
#end

B3 = projective_space(NormalToricVariety, 3)
W = toric_line_bundle(2 * torusinvariant_prime_divisors(B3)[1])
w = generic_section(W; rng=our_rng)
Kbar = anticanonical_bundle(B3)
a10 = generic_section(Kbar; rng=our_rng)
a21 = generic_section(Kbar^2 * W^(-1); rng=our_rng)
a32 = generic_section(Kbar^3 * W^(-2); rng=our_rng)
a43 = generic_section(Kbar^4 * W^(-3); rng=our_rng)
a65 = 0
t = global_tate_model(B3, [a10, a21 * w, a32 * w^2, a43 * w^3, a65 * w^5])
add_resolution!(
  t,
  [["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]],
  ["e1", "e4", "e2", "e3", "s"],
)
toric_centers = [["x", "y", "x1"], ["y", "toric_e1"]]
toric_exceptionals = ["toric_e1", "toric_e2"]
add_resolution!(t, toric_centers, toric_exceptionals)
add_weighted_resolution!(t, [(["x", "y", "w"], [1, 1, 1])], ["weighted_e"])
explicit_model_sections(t)["w"] = w

# Check stable model data is preserved, while caches and resolution data are invalidated.
defining_classes(t)["W"] = toric_divisor_class(W)
set_attribute!(t, :model_description => "attribute audit")
set_attribute!(t, :chern_classes => nothing)
set_attribute!(t, :exceptional_divisor_indices, Int[])
section_classes = classes_of_model_sections(t)
tate_discriminant = discriminant(t)
original_resolutions = deepcopy(resolutions(t))
original_weighted_resolutions = deepcopy(weighted_resolutions(t))
original_exceptional_indices = copy(exceptional_divisor_indices(t))
t_res = resolve(t, 1)

@testset "Custom blowup of a global Tate model" begin
  @test ambient_space(t_res) isa CoveredScheme{QQField}
  @test defining_classes(t_res) == defining_classes(t)
  @test model_description(t_res) == "attribute audit"
  @test !has_attribute(t_res, :chern_classes)
  @test classes_of_model_sections(t_res) === section_classes
  @test discriminant(t_res) === tate_discriminant
  @test !has_attribute(t_res, :weierstrass_model)
  @test resolutions(t) == original_resolutions
  @test weighted_resolutions(t) == original_weighted_resolutions
  @test exceptional_divisor_indices(t) == original_exceptional_indices
  @test_throws ArgumentError exceptional_divisor_indices(t_res)
  @test !has_attribute(t_res, :resolutions)
  @test !has_attribute(t_res, :weighted_resolutions)
  @test_throws ArgumentError("No resolutions known for this model") resolve(t_res, 1)
  @test_throws ArgumentError("A resolution must contain at least one blowup") add_resolution!(
    t, Vector{Vector{String}}(), String[]
  )
end

# Check that planned toric blowups are reused and agree with successive transforms.
@testset "Batched toric resolution" begin
  R = coordinate_ring(ambient_space(t))
  single_model = blow_up(
    t, ideal([R[5], R[6], R[1]]); coordinate_name="single_e"
  )
  @test ambient_space(single_model) isa NormalToricVariety
  @test exceptional_divisor_indices(single_model) == [
    ngens(coordinate_ring(ambient_space(single_model)))
  ]
  @test discriminant(single_model) === tate_discriminant
  @test !has_attribute(single_model, :resolutions)

  blowups, non_toric_center = Oscar._toric_blow_up_plan(
    t, toric_centers, toric_exceptionals
  )
  @test isnothing(non_toric_center)

  transform = tate_polynomial(t)
  transform = foldl((f, b) -> strict_transform(b, f), blowups; init=transform)
  @test Oscar._strict_transform_toric_sequence(blowups, tate_polynomial(t)) == transform

  batched_model = resolve(t, 2)
  @test evaluate(tate_polynomial(batched_model), gens(parent(transform))) == transform
  @test coordinate_names(ambient_space(batched_model)) ==
    coordinate_names(domain(last(blowups)))
  @test exceptional_divisor_indices(batched_model) ==
    index_of_exceptional_ray.(blowups)
  @test classes_of_model_sections(batched_model) === section_classes
  @test discriminant(batched_model) === tate_discriminant
  @test !has_attribute(batched_model, :exceptional_classes)
  @test length(exceptional_classes(batched_model)) == length(blowups)
  @test has_attribute(batched_model, :exceptional_classes)
  @test !has_attribute(batched_model, :weierstrass_model)

  combined_blow_down = get_attribute(batched_model, :blow_down_morphism)
  @test domain(combined_blow_down) === ambient_space(batched_model)
  @test codomain(combined_blow_down) === ambient_space(t)
  composed_blow_down = foldl(*, reverse(blowups))
  @test matrix(lattice_homomorphism(combined_blow_down)) ==
    matrix(lattice_homomorphism(composed_blow_down))

  B2 = projective_space(NormalToricVariety, 2)
  weierstrass = weierstrass_model(
    B2; completeness_check=false, rng=Random.Xoshiro(5678)
  )
  add_resolution!(weierstrass, toric_centers, toric_exceptionals)
  weierstrass_blowups, _ = Oscar._toric_blow_up_plan(
    weierstrass, toric_centers, toric_exceptionals
  )
  transform = weierstrass_polynomial(weierstrass)
  transform = foldl((f, b) -> strict_transform(b, f), weierstrass_blowups; init=transform)
  batched_weierstrass = resolve(weierstrass, 1)
  @test evaluate(
    weierstrass_polynomial(batched_weierstrass), gens(parent(transform))
  ) == transform

  mixed_centers = [toric_centers[1], ["y", "toric_e1*x"]]
  mixed_exceptionals = ["toric_e1", "mixed_e2"]
  mixed_blowups, non_toric_center = Oscar._toric_blow_up_plan(
    t, mixed_centers, mixed_exceptionals
  )
  @test length(mixed_blowups) == 1
  @test !isnothing(non_toric_center)

  add_resolution!(t, mixed_centers, mixed_exceptionals)
  mixed_model = resolve(t, 3)
  @test ambient_space(mixed_model) isa CoveredScheme{QQField}
  @test !has_attribute(mixed_model, :resolutions)
  @test classes_of_model_sections(mixed_model) === section_classes
  @test discriminant(mixed_model) === tate_discriminant
end
