#############################################################
# Extra long singular loci tests
#############################################################

using Random
our_rng = Random.Xoshiro(1234)

B3 = projective_space(NormalToricVariety, 3)
Kbar = anticanonical_divisor(B3)
foah1_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.4", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah2_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.12", type = "weierstrass", base_space = B3, defining_classes = Dict("b7" => Kbar, "b9" => Kbar), completeness_check = false)
foah3_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.54", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah4_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.17", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah5_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.73", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah6_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.82", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah7_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.96", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah8_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.106", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah9_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.118", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah10_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.130", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah11_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.142", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah12_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.155", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah13_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.181", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah14_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.168", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah15_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.190", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah16_B3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.203", type = "weierstrass", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)

@testset "Singular loci for Weierstrass forms of models in F-theory on all toric hypersurfaces, defined over a concrete base" begin
  @test [k[2:3] for k in singular_loci(foah1_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1")]
  @test [k[2:3] for k in singular_loci(foah2_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1")]
  @test [k[2:3] for k in singular_loci(foah3_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1")]
  @test [k[2:3] for k in singular_loci(foah4_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2")]
  @test [k[2:3] for k in singular_loci(foah5_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1")]
  @test [k[2:3] for k in singular_loci(foah6_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2")]
  @test [k[2:3] for k in singular_loci(foah7_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1")]
  @test [k[2:3] for k in singular_loci(foah8_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 2), "Non-split I_2")]
  @test [k[2:3] for k in singular_loci(foah9_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2")]
  @test [k[2:3] for k in singular_loci(foah10_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 3), "Split I_3")]
  @test [k[2:3] for k in singular_loci(foah11_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 3), "Split I_3")]
  @test [k[2:3] for k in singular_loci(foah12_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 2), "Non-split I_2")]
  @test [k[2:3] for k in singular_loci(foah13_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 4), "Split I_4")]
  @test [k[2:3] for k in singular_loci(foah14_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 3), "Split I_3")]
  @test [k[2:3] for k in singular_loci(foah15_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 2), "Non-split I_2"), ((0, 0, 2), "Non-split I_2")]
  @test [k[2:3] for k in singular_loci(foah16_B3_weier; rng = our_rng)] == [((0, 0, 1), "I_1"), ((0, 0, 3), "Split I_3"), ((0, 0, 3), "Split I_3"), ((0, 0, 3), "Split I_3")]
end

# rings needed for constructions
istar0_s_auxiliary_base_ring, (a1pp, a2pp, a3pp, a4pp, a6pp, vp, wp, mp) = QQ[:a1pp, :a2pp, :a3pp, :a4pp, :a6pp, :vp, :wp, :mp]
tate_auxiliary_base_ring, (a1p, a2p, a3p, a4p, a6p, v, w) = QQ[:a1p, :a2p, :a3p, :a4p, :a6p, :v, :w]

# construct Tate models over arbitrary base
t_istar0_s = global_tate_model(istar0_s_auxiliary_base_ring, [1 2 3 4 6 0 0 2; -2 -4 -4 -6 -8 1 2 -2], 3, [a1pp * (vp^2 + wp)^1, mp * (vp^2 + wp)^1 + a2pp * (vp^2 + wp)^2, a3pp * (vp^2 + wp)^2, mp^2 * (vp^2 + wp)^2 + a4pp * (vp^2 + wp)^3, a6pp * (vp^2 + wp)^4])
t_i1 = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -2 -2 -2 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^1, a4p * (v^2 + w)^1, a6p * (v^2 + w)^1])
t_i2_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -2 -2 -4 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^1, a4p * (v^2 + w)^1, a6p * (v^2 + w)^2])
t_i2_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -2 -2 -2 -4 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^1, a6p * (v^2 + w)^2])
t_i3_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -4 -4 -6 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^2, a4p * (v^2 + w)^2, a6p * (v^2 + w)^3])
t_i3_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -2 -2 -4 -6 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^2, a6p * (v^2 + w)^3])
t_i4_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -4 -4 -8 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^2, a4p * (v^2 + w)^2, a6p * (v^2 + w)^4])
t_i4_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -2 -4 -4 -8 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^1, a3p * (v^2 + w)^2, a4p * (v^2 + w)^2, a6p * (v^2 + w)^4])
t_i5_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -6 -6 -10 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^5])
t_i5_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -1 -2 -3 -5 1 2], 3, [a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5])
t_i6_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -6 -6 -12 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^6])
t_i6_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -2 -6 -6 -12 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^6])
t_i7_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 0 -8 -8 -14 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^0, a3p * (v^2 + w)^4, a4p * (v^2 + w)^4, a6p * (v^2 + w)^7])
t_i7_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; 0 -2 -6 -8 -14 1 2], 3, [a1p * (v^2 + w)^0, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^4, a6p * (v^2 + w)^7])
t_ii = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -2 -2 -2 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^1, a6p * (v^2 + w)^1])
t_iii = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -2 -2 -4 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^1, a6p * (v^2 + w)^2])
t_iv_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -2 -4 -4 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^2, a6p * (v^2 + w)^2])
t_iv_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -2 -4 -6 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^1, a4p * (v^2 + w)^2, a6p * (v^2 + w)^3])
t_istar0_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -4 -4 -6 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^2, a4p * (v^2 + w)^2, a6p * (v^2 + w)^3])
t_istar0_ss = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -4 -4 -8 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^2, a4p * (v^2 + w)^2, a6p * (v^2 + w)^4])
t_istar1_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -4 -6 -8 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^2, a4p * (v^2 + w)^3, a6p * (v^2 + w)^4])
t_istar1_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -4 -6 -10 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^2, a4p * (v^2 + w)^3, a6p * (v^2 + w)^5])
t_istar2_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -6 -6 -10 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^5])
t_istar2_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -6 -6 -12 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^6])
t_istar3_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -6 -8 -12 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^4, a6p * (v^2 + w)^6])
t_istar3_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -6 -8 -14 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^3, a4p * (v^2 + w)^4, a6p * (v^2 + w)^7])
t_istar4_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -8 -8 -14 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^4, a4p * (v^2 + w)^4, a6p * (v^2 + w)^7])
t_istar4_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -2 -8 -8 -16 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^1, a3p * (v^2 + w)^4, a4p * (v^2 + w)^4, a6p * (v^2 + w)^8])
t_ivstar_ns = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -4 -4 -6 -8 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^2, a3p * (v^2 + w)^2, a4p * (v^2 + w)^3, a6p * (v^2 + w)^4])
t_ivstar_s = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -4 -4 -6 -10 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^2, a3p * (v^2 + w)^2, a4p * (v^2 + w)^3, a6p * (v^2 + w)^5])
t_iiistar = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -4 -6 -6 -10 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^2, a3p * (v^2 + w)^3, a4p * (v^2 + w)^3, a6p * (v^2 + w)^5])
t_iistar = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -4 -6 -8 -10 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^2, a3p * (v^2 + w)^3, a4p * (v^2 + w)^4, a6p * (v^2 + w)^5])
t_nm = global_tate_model(tate_auxiliary_base_ring, [1 2 3 4 6 0 0; -2 -4 -6 -8 -12 1 2], 3, [a1p * (v^2 + w)^1, a2p * (v^2 + w)^2, a3p * (v^2 + w)^3, a4p * (v^2 + w)^4, a6p * (v^2 + w)^6])

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
