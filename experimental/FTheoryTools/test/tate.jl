#############################################################
# 1: Global Tate models over concrete base space
#############################################################

# sample base
base = sample_toric_variety()

# sections, used to test error messages
sec_a1 = generic_section(anticanonical_bundle(projective_space(NormalToricVariety,3)))
sec_a2 = generic_section(anticanonical_bundle(base)^2)
sec_a3 = generic_section(anticanonical_bundle(base)^3)
sec_a4 = generic_section(anticanonical_bundle(base)^4)
sec_a6 = generic_section(anticanonical_bundle(base)^6)

# construct and test one tate model
t = global_tate_model(base; completeness_check = false)

@testset "Attributes of global Tate models over concrete base space" begin
  @test parent(tate_section_a1(t)) == cox_ring(base_space(t))
  @test parent(tate_section_a2(t)) == cox_ring(base_space(t))
  @test parent(tate_section_a3(t)) == cox_ring(base_space(t))
  @test parent(tate_section_a4(t)) == cox_ring(base_space(t))
  @test parent(tate_section_a6(t)) == cox_ring(base_space(t))
  @test parent(tate_polynomial(t)) == cox_ring(ambient_space(t))
  @test parent(discriminant(t)) == cox_ring(base_space(t))
  @test dim(base_space(t)) == 3
  @test dim(ambient_space(t)) == 5
  @test base_fully_specified(t) == true
  @test base_fully_specified(t) == base_fully_specified(weierstrass_model(t))
  @test is_smooth(ambient_space(t)) == false
  @test toric_variety(calabi_yau_hypersurface(t)) == ambient_space(t)

  mktempdir() do path
    test_save_load_roundtrip(path, t) do loaded
      @test tate_polynomial(t) == tate_polynomial(loaded)
      @test tate_section_a1(t) == tate_section_a1(loaded)
      @test tate_section_a2(t) == tate_section_a2(loaded)
      @test tate_section_a3(t) == tate_section_a3(loaded)
      @test tate_section_a4(t) == tate_section_a4(loaded)
      @test tate_section_a6(t) == tate_section_a6(loaded)
      @test base_space(t) == base_space(loaded)
      @test ambient_space(t) == ambient_space(loaded)
    end
  end
end

@testset "Error messages in global Tate models over concrete base space" begin
  @test_throws ArgumentError global_tate_model(base, [sec_a2, sec_a3, sec_a4, sec_a6]; completeness_check = false)
  @test_throws ArgumentError global_tate_model(base, [sec_a1, sec_a2, sec_a3, sec_a4, sec_a6]; completeness_check = false)
end


#############################################################
# 2: Global Tate models over arbitrary base space
#############################################################

# rings needed for constructions
istar0_s_auxiliary_base_ring, (a1pp, a2pp, a3pp, a4pp, a6pp, vp, mp) = QQ["a1pp", "a2pp", "a3pp", "a4pp", "a6pp", "vp", "mp"];
tate_auxiliary_base_ring, (a1p, a2p, a3p, a4p, a6p, v) = QQ["a1p", "a2p", "a3p", "a4p", "a6p", "v"];

# construct Tate models over arbitrary base
t_istar0_s = global_tate_model(istar0_s_auxiliary_base_ring, [1 2 3 4 6 0 2; -1 -2 -2 -3 -4 1 -1], 3, [a1pp * vp^1, mp * vp^1 + a2pp * vp^2, a3pp * vp^2, mp^2 * vp^2 + a4pp * vp^3, a6pp * vp^4]);
t_i1 = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 0 -1 -1 -1 1], 3, [a1p * v^0, a2p * v^0, a3p * v^1, a4p * v^1, a6p * v^1]);
t_i2_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 0 -1 -1 -2 1], 3, [a1p * v^0, a2p * v^0, a3p * v^1, a4p * v^1, a6p * v^2]);
t_i2_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 -1 -1 -1 -2 1], 3, [a1p * v^0, a2p * v^1, a3p * v^1, a4p * v^1, a6p * v^2]);
t_i3_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 0 -2 -2 -3 1], 3, [a1p * v^0, a2p * v^0, a3p * v^2, a4p * v^2, a6p * v^3]);
t_i3_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 -1 -1 -2 -3 1], 3, [a1p * v^0, a2p * v^1, a3p * v^1, a4p * v^2, a6p * v^3]);
t_i4_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 0 -2 -2 -4 1], 3, [a1p * v^0, a2p * v^0, a3p * v^2, a4p * v^2, a6p * v^4]);
t_i4_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 -1 -2 -2 -4 1], 3, [a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^2, a6p * v^4]);
t_i5_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 0 -3 -3 -5 1], 3, [a1p * v^0, a2p * v^0, a3p * v^3, a4p * v^3, a6p * v^5]);
t_i5_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 -1 -2 -3 -5 1], 3, [a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5]);
t_i6_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 0 -3 -3 -6 1], 3, [a1p * v^0, a2p * v^0, a3p * v^3, a4p * v^3, a6p * v^6]);
t_i6_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 -1 -3 -3 -6 1], 3, [a1p * v^0, a2p * v^1, a3p * v^3, a4p * v^3, a6p * v^6]);
t_i7_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 0 -4 -4 -7 1], 3, [a1p * v^0, a2p * v^0, a3p * v^4, a4p * v^4, a6p * v^7]);
t_i7_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 -1 -3 -4 -7 1], 3, [a1p * v^0, a2p * v^1, a3p * v^3, a4p * v^4, a6p * v^7]);
t_ii = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -1 -1 -1 1], 3, [a1p * v^1, a2p * v^1, a3p * v^1, a4p * v^1, a6p * v^1]);
t_iii = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -1 -1 -2 1], 3, [a1p * v^1, a2p * v^1, a3p * v^1, a4p * v^1, a6p * v^2]);
t_iv_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -1 -2 -2 1], 3, [a1p * v^1, a2p * v^1, a3p * v^1, a4p * v^2, a6p * v^2]);
t_iv_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -1 -2 -3 1], 3, [a1p * v^1, a2p * v^1, a3p * v^1, a4p * v^2, a6p * v^3]);
t_istar0_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -2 -2 -3 1], 3, [a1p * v^1, a2p * v^1, a3p * v^2, a4p * v^2, a6p * v^3]);
t_istar0_ss = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -2 -2 -4 1], 3, [a1p * v^1, a2p * v^1, a3p * v^2, a4p * v^2, a6p * v^4]);
t_istar1_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -2 -3 -4 1], 3, [a1p * v^1, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^4]);
t_istar1_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -2 -3 -5 1], 3, [a1p * v^1, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5]);
t_istar2_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -3 -3 -5 1], 3, [a1p * v^1, a2p * v^1, a3p * v^3, a4p * v^3, a6p * v^5]);
t_istar2_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -3 -3 -6 1], 3, [a1p * v^1, a2p * v^1, a3p * v^3, a4p * v^3, a6p * v^6]);
t_istar3_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -3 -4 -6 1], 3, [a1p * v^1, a2p * v^1, a3p * v^3, a4p * v^4, a6p * v^6]);
t_istar3_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -3 -4 -7 1], 3, [a1p * v^1, a2p * v^1, a3p * v^3, a4p * v^4, a6p * v^7]);
t_istar4_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -4 -4 -7 1], 3, [a1p * v^1, a2p * v^1, a3p * v^4, a4p * v^4, a6p * v^7]);
t_istar4_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -1 -4 -4 -8 1], 3, [a1p * v^1, a2p * v^1, a3p * v^4, a4p * v^4, a6p * v^8]);
t_ivstar_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -2 -2 -3 -4 1], 3, [a1p * v^1, a2p * v^2, a3p * v^2, a4p * v^3, a6p * v^4]);
t_ivstar_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -2 -2 -3 -5 1], 3, [a1p * v^1, a2p * v^2, a3p * v^2, a4p * v^3, a6p * v^5]);
t_iiistar = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -2 -3 -3 -5 1], 3, [a1p * v^1, a2p * v^2, a3p * v^3, a4p * v^3, a6p * v^5]);
t_iistar = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -2 -3 -4 -5 1], 3, [a1p * v^1, a2p * v^2, a3p * v^3, a4p * v^4, a6p * v^5]);
t_nm = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; -1 -2 -3 -4 -6 1], 3, [a1p * v^1, a2p * v^2, a3p * v^3, a4p * v^4, a6p * v^6]);

@testset "Attributes of global Tate models over generic base space" begin
  @test parent(tate_section_a1(t_i5_s)) == cox_ring(base_space(t_i5_s))
  @test parent(tate_section_a2(t_i5_s)) == cox_ring(base_space(t_i5_s))
  @test parent(tate_section_a3(t_i5_s)) == cox_ring(base_space(t_i5_s))
  @test parent(tate_section_a4(t_i5_s)) == cox_ring(base_space(t_i5_s))
  @test parent(tate_section_a6(t_i5_s)) == cox_ring(base_space(t_i5_s))
  @test parent(tate_polynomial(t_i5_s)) == cox_ring(ambient_space(t_i5_s))
  @test parent(discriminant(t_i5_s)) == cox_ring(base_space(t_i5_s))
  @test dim(base_space(t_i5_s)) == 3
  @test dim(ambient_space(t_i5_s)) == 5
  @test base_fully_specified(t_i5_s) == false
  @test base_fully_specified(t_i5_s) == base_fully_specified(weierstrass_model(t_i5_s))
  @test is_smooth(ambient_space(t_i5_s)) == false
  @test toric_variety(calabi_yau_hypersurface(t_i5_s)) == ambient_space(t_i5_s)
end

@testset "Error messages in global Tate models over generic base space" begin
  @test_throws ArgumentError global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 -1 -2 -3 -5 1], 3, [a1p])
  @test_throws ArgumentError global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 -1 -2 -3 -5 1], 3, [a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, sec_a6])
  @test_throws ArgumentError global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 -1 -2 -3 -5 1], -1, [a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5])
  @test_throws ArgumentError global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0; 0 -1 -2 -3 -5 1], 7, [a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5])
end

@testset "Singular loci of global Tate models over generic base space" begin
  @test length(singular_loci(t_i1)) == 2
  @test length(singular_loci(t_i2_ns)) == 2
  @test length(singular_loci(t_i2_s)) == 2
  @test length(singular_loci(t_i3_ns)) == 2
  @test length(singular_loci(t_i3_s)) == 2
  @test length(singular_loci(t_i4_ns)) == 2
  @test length(singular_loci(t_i4_s)) == 2
  @test length(singular_loci(t_i5_ns)) == 2
  @test length(singular_loci(t_i5_s)) == 2
  @test length(singular_loci(t_i6_ns)) == 2
  @test length(singular_loci(t_i6_s)) == 2
  @test length(singular_loci(t_i7_ns)) == 2
  @test length(singular_loci(t_i7_s)) == 2
  @test length(singular_loci(t_ii)) == 2
  @test length(singular_loci(t_iii)) == 2
  @test length(singular_loci(t_iv_ns)) == 2
  @test length(singular_loci(t_iv_s)) == 2
  @test length(singular_loci(t_istar0_ns)) == 2
  @test length(singular_loci(t_istar0_ss)) == 2
  @test length(singular_loci(t_istar0_s)) == 2
  @test length(singular_loci(t_istar1_ns)) == 2
  @test length(singular_loci(t_istar1_s)) == 2
  @test length(singular_loci(t_istar2_ns)) == 2
  @test length(singular_loci(t_istar2_s)) == 2
  @test length(singular_loci(t_istar3_ns)) == 2
  @test length(singular_loci(t_istar3_s)) == 2
  @test length(singular_loci(t_istar4_ns)) == 2
  @test length(singular_loci(t_istar4_s)) == 2
  @test length(singular_loci(t_ivstar_ns)) == 2
  @test length(singular_loci(t_ivstar_s)) == 2
  @test length(singular_loci(t_iiistar)) == 2
  @test length(singular_loci(t_iistar)) == 2
  @test length(singular_loci(t_nm)) == 2
  @test singular_loci(t_i1)[1][2:3] == ((0, 0, 1), "I_1")
  @test singular_loci(t_i2_ns)[2][2:3] == ((0, 0, 2), "Non-split I_2")
  @test singular_loci(t_i2_s)[2][2:3] == ((0, 0, 2), "Split I_2")
  @test singular_loci(t_i3_ns)[2][2:3] == ((0, 0, 3), "Non-split I_3")
  @test singular_loci(t_i3_s)[2][2:3] == ((0, 0, 3), "Split I_3")
  @test singular_loci(t_i4_ns)[2][2:3] == ((0, 0, 4), "Non-split I_4")
  @test singular_loci(t_i4_s)[2][2:3] == ((0, 0, 4), "Split I_4")
  @test singular_loci(t_i5_ns)[2][2:3] == ((0, 0, 5), "Non-split I_5")
  @test singular_loci(t_i5_s)[2][2:3] == ((0, 0, 5), "Split I_5")
  @test singular_loci(t_i6_ns)[2][2:3] == ((0, 0, 6), "Non-split I_6")
  @test singular_loci(t_i6_s)[2][2:3] == ((0, 0, 6), "Split I_6")
  @test singular_loci(t_i7_ns)[2][2:3] == ((0, 0, 7), "Non-split I_7")
  @test singular_loci(t_i7_s)[2][2:3] == ((0, 0, 7), "Split I_7")
  @test singular_loci(t_ii)[2][2:3] == ((1, 1, 2), "II")
  @test singular_loci(t_iii)[2][2:3] == ((1, 2, 3), "III")
  @test singular_loci(t_iv_ns)[2][2:3] == ((2, 2, 4), "Non-split IV")
  @test singular_loci(t_iv_s)[2][2:3] == ((2, 2, 4), "Split IV")
  @test singular_loci(t_istar0_ns)[2][2:3] == ((2, 3, 6), "Non-split I^*_0")
  @test singular_loci(t_istar0_ss)[2][2:3] == ((2, 3, 6), "Semi-split I^*_0")
  @test singular_loci(t_istar0_s)[2][2:3] == ((2, 3, 6), "Split I^*_0")
  @test singular_loci(t_istar1_ns)[2][2:3] == ((2, 3, 7), "Non-split I^*_1")
  @test singular_loci(t_istar1_s)[2][2:3] == ((2, 3, 7), "Split I^*_1")
  @test singular_loci(t_istar2_ns)[2][2:3] == ((2, 3, 8), "Non-split I^*_2")
  @test singular_loci(t_istar2_s)[2][2:3] == ((2, 3, 8), "Split I^*_2")
  @test singular_loci(t_istar3_ns)[2][2:3] == ((2, 3, 9), "Non-split I^*_3")
  @test singular_loci(t_istar3_s)[2][2:3] == ((2, 3, 9), "Split I^*_3")
  @test singular_loci(t_istar4_ns)[2][2:3] == ((2, 3, 10), "Non-split I^*_4")
  @test singular_loci(t_istar4_s)[2][2:3] == ((2, 3, 10), "Split I^*_4")
  @test singular_loci(t_ivstar_ns)[2][2:3] == ((3, 4, 8), "Non-split IV^*")
  @test singular_loci(t_ivstar_s)[2][2:3] == ((3, 4, 8), "Split IV^*")
  @test singular_loci(t_iiistar)[2][2:3] == ((3, 5, 9), "III^*")
  @test singular_loci(t_iistar)[2][2:3] == ((4, 5, 10), "II^*")
  @test singular_loci(t_nm)[1][2:3] == ((4, 6, 12), "Non-minimal")
end

@testset "Blowups of global Tate models" begin
  id_i5_s = ideal([tate_polynomial(t_i5_s)]);
  tas = ambient_space(t_i5_s);
  irr_i5_s = irrelevant_ideal(tas);
  sri_i5_s = stanley_reisner_ideal(tas);
  lin_i5_s = ideal_of_linear_relations(tas);
  id_fin,  = _blowup_global_sequence(id_i5_s, [[7, 8, 6], [2, 3, 1], [3, 4], [2, 4]], irr_i5_s, sri_i5_s, lin_i5_s)
  @test string(gens(id_fin)[end]) == "-b_4_1*b_2_1*a1p*z - b_4_1*b_2_2 - b_4_1*b_2_3*b_1_3^2*a3p*z^3 + b_4_2*b_3_2*b_2_1^2*b_1_1 + b_4_2*b_3_2*b_2_1^2*b_1_3*a2p*z^2 + b_4_2*b_3_2*b_2_1*b_2_3*b_1_3^3*a4p*z^4 + b_4_2*b_3_2*b_2_3^2*b_1_3^5*a6p*z^6"
end

#@testset "Fibers" begin
#  inters = analyze_fibers(t_i5_s, [[7, 8, 6], [2, 3, 1], [3, 4], [2, 4]])
#  @test string(inters[1][2][1][2][1]) == "ideal(e_2*b_2_1 - b_1_1, b_3_1*a1p*z - b_3_2*b_1_1^2, b_4_1*a1p*z - b_4_2*b_3_2*b_2_1*b_1_1, b_4_1*b_1_1 - b_4_2*b_3_1*b_2_1, b_4_1*e_2 - b_4_2*b_3_1, e_4*b_4_2 - e_2, e_4*b_4_1 - b_3_1, y, x, v, b_1_3, b_1_2, e_1, b_2_3, b_2_2, e_3)"
#end
