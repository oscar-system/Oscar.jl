#############################################################
# 1: Global Tate models over concrete base space
#############################################################

using Random
our_rng = Random.Xoshiro(1234)

my_base = projective_space(NormalToricVariety, 3)
sec_a1 = generic_section(anticanonical_bundle(projective_space(NormalToricVariety,3)))
sec_a2 = generic_section(anticanonical_bundle(my_base)^2)
sec_a3 = generic_section(anticanonical_bundle(my_base)^3)
sec_a4 = generic_section(anticanonical_bundle(my_base)^4)
sec_a6 = generic_section(anticanonical_bundle(my_base)^6)
t = global_tate_model(my_base; completeness_check = false)

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
  @test is_base_space_fully_specified(t) == true
  @test is_base_space_fully_specified(t) == is_base_space_fully_specified(weierstrass_model(t))
  @test is_smooth(ambient_space(t)) == false
  @test toric_variety(calabi_yau_hypersurface(t)) == ambient_space(t)
  @test is_partially_resolved(t) == false
end

@testset "Error messages in global Tate models over concrete base space" begin
  @test_throws ArgumentError global_tate_model(my_base, [sec_a2, sec_a3, sec_a4, sec_a6]; completeness_check = false)
  @test_throws ArgumentError global_tate_model(my_base, [sec_a1, sec_a2, sec_a3, sec_a4, sec_a6]; completeness_check = false)
end
@test_skip @testset "Serialization Tests" begin
  B3 = projective_space(NormalToricVariety, 3)
  w = torusinvariant_prime_divisors(B3)[1]
  t2 = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)

  @testset "Saving and loading global Tate model over concrete base space" begin
    mktempdir() do path
      test_save_load_roundtrip(path, t2) do loaded
        @test tate_polynomial(t2) == tate_polynomial(loaded)
        @test tate_section_a1(t2) == tate_section_a1(loaded)
        @test tate_section_a2(t2) == tate_section_a2(loaded)
        @test tate_section_a3(t2) == tate_section_a3(loaded)
        @test tate_section_a4(t2) == tate_section_a4(loaded)
        @test tate_section_a6(t2) == tate_section_a6(loaded)
        @test base_space(t2) == base_space(loaded)
        @test ambient_space(t2) == ambient_space(loaded)
        @test is_base_space_fully_specified(t2) == is_base_space_fully_specified(loaded)
        for (key, value) in t2.__attrs
          if value isa String || value isa Vector{String} || value isa Bool
            @test t2.__attrs[key] == loaded.__attrs[key]
          end
        end
      end
    end
  end

  t2_copy = global_tate_model(my_base = projective_space(NormalToricVariety, 3); completeness_check = false)

  @testset "Saving and loading another global Tate model over concrete base space" begin
    mktempdir() do path
      test_save_load_roundtrip(path, t2_copy) do loaded
        @test tate_polynomial(t2_copy) == tate_polynomial(loaded)
        @test tate_section_a1(t2_copy) == tate_section_a1(loaded)
        @test tate_section_a2(t2_copy) == tate_section_a2(loaded)
        @test tate_section_a3(t2_copy) == tate_section_a3(loaded)
        @test tate_section_a4(t2_copy) == tate_section_a4(loaded)
        @test tate_section_a6(t2_copy) == tate_section_a6(loaded)
        @test base_space(t2_copy) == base_space(loaded)
        @test ambient_space(t2_copy) == ambient_space(loaded)
        @test is_base_space_fully_specified(t2_copy) == is_base_space_fully_specified(loaded)
        @test is_partially_resolved(t2_copy) == is_partially_resolved(loaded)
        @test explicit_model_sections(t2_copy) == explicit_model_sections(loaded)
        @test model_section_parametrization(t2_copy) == model_section_parametrization(loaded)
        @test defining_classes(t2_copy) == defining_classes(loaded)
      end
    end
  end
end



#############################################################
# 2: Global Tate models over arbitrary base space
#############################################################

# rings needed for constructions
istar0_s_auxiliary_base_ring, (a1pp, a2pp, a3pp, a4pp, a6pp, vp, wp, mp) = QQ[:a1pp, :a2pp, :a3pp, :a4pp, :a6pp, :vp, :wp, :mp];
tate_auxiliary_base_ring, (a1p, a2p, a3p, a4p, a6p, v, w) = QQ[:a1p, :a2p, :a3p, :a4p, :a6p, :v, :w];

# construct Tate models over arbitrary base
t_istar0_s = global_tate_model(istar0_s_auxiliary_base_ring, [1 2 3 4 6 0 0 2; -2 -4 -4 -6 -8 1 2 -2], 3, [a1pp * (vp^2 + wp)^1, mp * (vp^2 + wp)^1 + a2pp * (vp^2 + wp)^2, a3pp * (vp^2 + wp)^2, mp^2 * (vp^2 + wp)^2 + a4pp * (vp^2 + wp)^3, a6pp * (vp^2 + wp)^4]);
t_i1 = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -2 -2 -2 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^1, a4p * (v^2 + w)^1, a6p * (v^2 + w)^1]);
t_i2_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -2 -2 -4 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^1, a4p * (v^2 + w)^1, a6p * (v^2 + w)^2]);
t_i2_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -2 -2 -2 -4 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^1, a6p * (v^2 + w)^2]);
t_i3_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -4 -4 -6 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^2, a4p * (v^2 + w)^2, a6p * (v^2 + w)^3]);
t_i3_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -2 -2 -4 -6 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^2, a6p * (v^2 + w)^3]);
t_i4_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -4 -4 -8 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^2, a4p * (v^2 + w)^2, a6p * (v^2 + w)^4]);
t_i4_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -2 -4 -4 -8 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^1, a3p * (v^2 + w)^2, a4p * (v^2 + w)^2, a6p * (v^2 + w)^4]);
t_i5_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -6 -6 -10 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^5]);
t_i5_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -1 -2 -3 -5 1 2], 3, [a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5]);
t_i6_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -6 -6 -12 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^6]);
t_i6_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -2 -6 -6 -12 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^6]);
t_i7_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -8 -8 -14 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^4, a4p * (v^2 + w)^4, a6p * (v^2 + w)^7]);
t_i7_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -2 -6 -8 -14 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^4, a6p * (v^2 + w)^7]);
t_ii = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -2 -2 -2 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^1, a6p * (v^2 + w)^1]);
t_iii = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -2 -2 -4 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^1, a6p * (v^2 + w)^2]);
t_iv_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -2 -4 -4 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^2, a6p * (v^2 + w)^2]);
t_iv_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -2 -4 -6 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^2, a6p * (v^2 + w)^3]);
t_istar0_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -4 -4 -6 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^2, a4p * (v^2 + w)^2, a6p * (v^2 + w)^3]);
t_istar0_ss = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -4 -4 -8 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^2, a4p * (v^2 + w)^2, a6p * (v^2 + w)^4]);
t_istar1_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -4 -6 -8 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^2, a4p * (v^2 + w)^3, a6p * (v^2 + w)^4]);
t_istar1_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -4 -6 -10 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^2, a4p * (v^2 + w)^3, a6p * (v^2 + w)^5]);
t_istar2_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -6 -6 -10 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^5]);
t_istar2_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -6 -6 -12 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^6]);
t_istar3_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -6 -8 -12 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^4, a6p * (v^2 + w)^6]);
t_istar3_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -6 -8 -14 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^4, a6p * (v^2 + w)^7]);
t_istar4_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -8 -8 -14 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^4, a4p * (v^2 + w)^4, a6p * (v^2 + w)^7]);
t_istar4_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -8 -8 -16 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^4, a4p * (v^2 + w)^4, a6p * (v^2 + w)^8]);
t_ivstar_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -4 -4 -6 -8 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^2, a3p * (v^2 + w)^2, a4p * (v^2 + w)^3, a6p * (v^2 + w)^4]);
t_ivstar_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -4 -4 -6 -10 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^2, a3p * (v^2 + w)^2, a4p * (v^2 + w)^3, a6p * (v^2 + w)^5]);
t_iiistar = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -4 -6 -6 -10 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^2, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^5]);
t_iistar = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -4 -6 -8 -10 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^2, a3p * (v^2 + w)^3, a4p * (v^2 + w)^4, a6p * (v^2 + w)^5]);
t_nm = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -4 -6 -8 -12 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^2, a3p * (v^2 + w)^3, a4p * (v^2 + w)^4, a6p * (v^2 + w)^6]);

@testset "Attributes of global Tate models over generic base space" begin
  @test parent(tate_section_a1(t_i5_s)) == coordinate_ring(base_space(t_i5_s))
  @test parent(tate_section_a2(t_i5_s)) == coordinate_ring(base_space(t_i5_s))
  @test parent(tate_section_a3(t_i5_s)) == coordinate_ring(base_space(t_i5_s))
  @test parent(tate_section_a4(t_i5_s)) == coordinate_ring(base_space(t_i5_s))
  @test parent(tate_section_a6(t_i5_s)) == coordinate_ring(base_space(t_i5_s))
  @test parent(tate_polynomial(t_i5_s)) == coordinate_ring(ambient_space(t_i5_s))
  @test parent(discriminant(t_i5_s)) == coordinate_ring(base_space(t_i5_s))
  @test dim(base_space(t_i5_s)) == 3
  @test dim(ambient_space(t_i5_s)) == 5
  @test is_base_space_fully_specified(t_i5_s) == false
  @test is_base_space_fully_specified(t_i5_s) == is_base_space_fully_specified(weierstrass_model(t_i5_s))
end

@testset "Error messages in global Tate models over generic base space" begin
  @test_throws ArgumentError global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -1 -2 -3 -5 1 2], 3, [a1p])
  @test_throws ArgumentError global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -1 -2 -3 -5 1 2], 3, [a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, sec_a6])
  @test_throws ArgumentError global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -1 -2 -3 -5 1 2], -1, [a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5])
end

@testset "Singular loci of global Tate models over generic base space" begin
  @test [k[2:3] for k in singular_loci(t_i1; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 1), "I_1")]
  @test [k[2:3] for k in singular_loci(t_i2_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2")]
  @test [k[2:3] for k in singular_loci(t_i2_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Split I_2")]
  @test [k[2:3] for k in singular_loci(t_i3_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 3), "Non-split I_3")]
  @test [k[2:3] for k in singular_loci(t_i3_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 3), "Split I_3")]
  @test [k[2:3] for k in singular_loci(t_i4_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 4), "Non-split I_4")]
  @test [k[2:3] for k in singular_loci(t_i4_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 4), "Split I_4")]
  @test [k[2:3] for k in singular_loci(t_i5_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 5), "Non-split I_5")]
  @test [k[2:3] for k in singular_loci(t_i5_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 5), "Split I_5")]
  @test [k[2:3] for k in singular_loci(t_i6_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 6), "Non-split I_6")]
  @test [k[2:3] for k in singular_loci(t_i6_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 6), "Split I_6")]
  @test [k[2:3] for k in singular_loci(t_i7_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 7), "Non-split I_7")]
  @test [k[2:3] for k in singular_loci(t_i7_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 7), "Split I_7")]
  @test [k[2:3] for k in singular_loci(t_ii; rng = our_rng)] == [((0, 0, 1), "I_1"), ((1, 1, 2), "II")]
  @test [k[2:3] for k in singular_loci(t_iii; rng = our_rng)] == [((0, 0, 1), "I_1"), ((1, 2, 3), "III")]
  @test [k[2:3] for k in singular_loci(t_iv_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 2, 4), "Non-split IV")]
  @test [k[2:3] for k in singular_loci(t_iv_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 2, 4), "Split IV")]
  @test [k[2:3] for k in singular_loci(t_istar0_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 6), "Non-split I^*_0")]
  @test [k[2:3] for k in singular_loci(t_istar0_ss; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 6), "Semi-split I^*_0")]
  @test [k[2:3] for k in singular_loci(t_istar0_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 6), "Split I^*_0")]
  @test [k[2:3] for k in singular_loci(t_istar1_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 7), "Non-split I^*_1")]
  @test [k[2:3] for k in singular_loci(t_istar1_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 7), "Split I^*_1")]
  @test [k[2:3] for k in singular_loci(t_istar2_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 8), "Non-split I^*_2")]
  @test [k[2:3] for k in singular_loci(t_istar2_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 8), "Split I^*_2")]
  @test [k[2:3] for k in singular_loci(t_istar3_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 9), "Non-split I^*_3")]
  @test [k[2:3] for k in singular_loci(t_istar3_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 9), "Split I^*_3")]
  @test [k[2:3] for k in singular_loci(t_istar4_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 10), "Non-split I^*_4")]
  @test [k[2:3] for k in singular_loci(t_istar4_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((2, 3, 10), "Split I^*_4")]
  @test [k[2:3] for k in singular_loci(t_ivstar_ns; rng = our_rng)] == [((0, 0, 1), "I_1"), ((3, 4, 8), "Non-split IV^*")]
  @test [k[2:3] for k in singular_loci(t_ivstar_s; rng = our_rng)] == [((0, 0, 1), "I_1"), ((3, 4, 8), "Split IV^*")]
  @test [k[2:3] for k in singular_loci(t_iiistar; rng = our_rng)] == [((0, 0, 1), "I_1"), ((3, 5, 9), "III^*")]
  @test [k[2:3] for k in singular_loci(t_iistar; rng = our_rng)] == [((0, 0, 1), "I_1"), ((4, 5, 10), "II^*")]
  @test [k[2:3] for k in singular_loci(t_nm; rng = our_rng)] == [((0, 0, 1), "I_1"), ((4, 6, 12), "Non-minimal")]
end
