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
explicit_model_sections(t)["w"] = w
t_res = resolve(t, 1)

@testset "Custom blowup of a global Tate model" begin
  @test typeof(ambient_space(t_res)) == CoveredScheme{QQField}
end
