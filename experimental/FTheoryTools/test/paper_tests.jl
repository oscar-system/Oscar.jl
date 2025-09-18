# The following code is used in the publication "FTheoryTools: Advancing Computational Capabilities for F-Theory Research".
# See BMT25 in the OSCAR references for more details on this publication.
# A preprint is available at https://arxiv.org/abs/2506.13849.
# Consequently, both the input and output should ideally remain unchanged.

using Random
our_rng = Random.Xoshiro(1234)

#######################################
# 3.1 Weierstrass Models
#######################################

B2 = projective_space(NormalToricVariety, 2)
x1, x2, x3 = gens(cox_ring(B2))
weier_f =
  1//48 * (
    -(28 * x1 * x2^5 + 169 * x3^6)^2 +
    24 * (2 * x1^3 * (x2 + x3)^9 + 13 * x1^2 * x2^4 * x3^6)
  )
weier_g =
  1//864 * (
    216 * x1^4 * x2^8 * x3^6 + (28 * x1 * x2^5 + 169 * x3^6)^3 -
    36 * (28 * x1 * x2^5 + 169 * x3^6) * (2 * x1^3 * (x2 + x3)^9 + 13 * x1^2 * x2^4 * x3^6)
  )
w = weierstrass_model(B2, weier_f, weier_g; completeness_check=false)

@testset "FTheoryToolsPaper Section 3.1 - Test 1" begin
  @test dim(ambient_space(w)) == 4
  @test length(explicit_model_sections(w)) == 2
end

w_generic = weierstrass_model(B2; completeness_check=false, rng=our_rng)

@testset "FTheoryToolsPaper Section 3.1 - Test 2" begin
  @test length(explicit_model_sections(w_generic)) == 2
end

#######################################
# 3.2 The Refined Tate Fiber Type
#######################################

@testset "FTheoryToolsPaper Section 3.2" begin
  @test length(singular_loci(w_generic)) == 1
  @test length(singular_loci(w)) == 2
end

#######################################
# 3.3 Global Tate Models
#######################################

a1 = 13 * x3^3
a2 = 7 * x1 * x2^5
a3 = x1^2 * x2^4 * x3^3
a4 = x1^3 * (x2 + x3)^9
a6 = zero(cox_ring(B2))
t = global_tate_model(B2, [a1, a2, a3, a4, a6])

@testset "FTheoryToolsPaper Section 3.3" begin
  @test hypersurface_equation(t) == tate_polynomial(t)
  @test weierstrass_section_f(weierstrass_model(t)) == weierstrass_section_f(w)
  @test weierstrass_section_g(weierstrass_model(t)) == weierstrass_section_g(w)
end

#######################################
# 3.4 Hypersurface Models
#######################################

fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2, 3, 1])
set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])
D1 = 2 * anticanonical_divisor_class(B2)
D2 = 3 * anticanonical_divisor_class(B2)
amb_ring, (x1, x2, x3, x, y, z) = polynomial_ring(QQ, ["x1", "x2", "x3", "x", "y", "z"])
p =
  x^3 + 7 * x1 * x2^5 * x^2 * z^2 + x1^3 * (x2 + x3)^9 * x * z^4 - y^2 -
  13 * x3^3 * x * y * z -
  x1^2 * x2^4 * x3^3 * y * z^3
h = hypersurface_model(B2, fiber_ambient_space, [D1, D2], p; completeness_check=false)

@testset "FTheoryToolsPaper Section 3.4" begin
  @test string(hypersurface_equation(h)) == string(tate_polynomial(t))
  @test rank(grading_group(cox_ring(ambient_space(h)))) ==
    rank(grading_group(cox_ring(ambient_space(t))))
end

#######################################
# 3.5 Resolution of Singularities
#######################################

B3 = projective_space(NormalToricVariety, 3)
Kbar = anticanonical_bundle(B3)
W = toric_line_bundle(torusinvariant_prime_divisors(B3)[1])
w = generic_section(W; rng=our_rng)
w = gens(cox_ring(B3))[1]
a10 = generic_section(Kbar; rng=our_rng)
a21 = generic_section(Kbar^2 * W^(-1); rng=our_rng)
a32 = generic_section(Kbar^3 * W^(-2); rng=our_rng)
a43 = generic_section(Kbar^4 * W^(-3); rng=our_rng)
a6 = zero(cox_ring(B3))
t = global_tate_model(B3, [a10, a21 * w, a32 * w^2, a43 * w^3, a6])
amb = ambient_space(t)
cox_ring(amb)
t1 = blow_up(t, ["x", "y", string(w)]; coordinate_name="e1")

@testset "FTheoryToolsPaper Section 3.5.1 Toric Resolution Part 1" begin
  @test length(singular_loci(t)) == 2
  @test_throws ArgumentError resolutions(t)
end

add_resolution!(
  t,
  [["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]],
  ["e1", "e4", "e2", "e3", "s"],
)
explicit_model_sections(t)["w"] = w
t_res = resolve(t, 1)
tate_polynomial(t_res)

@testset "FTheoryToolsPaper Section 3.5.1 Toric Resolution Part 2" begin
  @test is_partially_resolved(t_res)
end

W = toric_line_bundle(2 * torusinvariant_prime_divisors(B3)[1])
w = generic_section(W; rng=our_rng)
a10 = generic_section(Kbar; rng=our_rng)
a21 = generic_section(Kbar^2 * W^(-1); rng=our_rng)
a32 = generic_section(Kbar^3 * W^(-2); rng=our_rng)
a43 = generic_section(Kbar^4 * W^(-3); rng=our_rng)
a6 = zero(cox_ring(B3))
t2 = global_tate_model(B3, [a10, a21 * w, a32 * w^2, a43 * w^3, a6])

@testset "FTheoryToolsPaper Section 3.5.2 Non-Toric Resolution Part 1" begin
  @test length(singular_loci(t2)) == 2
end

t2_1 = blow_up(t, ["x", "y", string(w)]; coordinate_name="e1")
ambient_space(t2_1)

@testset "FTheoryToolsPaper Section 3.5.2 Non-Toric Resolution Part 2" begin
  @test typeof(ambient_space(t2_1)) == CoveredScheme{QQField}
end

add_resolution!(
  t2,
  [["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]],
  ["e1", "e4", "e2", "e3", "s"],
)
explicit_model_sections(t2)["w"] = w
t2_res = resolve(t2, 1)

@testset "FTheoryToolsPaper Section 3.5.2 Non-Toric Resolution Part 3" begin
  @test_throws ArgumentError tate_polynomial(t2_res)
  @test typeof(tate_ideal_sheaf(t2_res)) ==
    Oscar.StrictTransformIdealSheaf{CoveredScheme{QQField},AbsAffineScheme,Ideal,Map}
end

#######################################
# 4.1 Fundamentals of Literature Models
#######################################

W = torusinvariant_prime_divisors(B3)[1]
t = literature_model(;
  arxiv_id="1109.3454",
  equation="3.1",
  base_space=B3,
  defining_classes=Dict("w" => W),
  rng=our_rng,
)

@testset "FTheoryToolsPaper Section 4.1 Part 1" begin
  @test journal_name(t) == "Nucl. Phys. B"
  @test paper_authors(t) == ["Sven Krause", "Christoph Mayrhofer", "Timo Weigand"]
  @test arxiv_doi(t) == "10.48550/arXiv.1109.3454"
  @test journal_model_equation_number(t) == "3.1"
  @test journal_model_page(t) == "9"
  @test length(defining_classes(t)) == 1
  @test length(model_sections(t)) == 9
  @test length(model_section_parametrization(t)) == 4
  @test length(tunable_sections(t)) == 5
  @test length(classes_of_model_sections(t)) == 9
  @test length(classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes(t)) == 5
  @test length(explicit_model_sections(t)) == 9
  @test length(resolutions(t)) == 1
end

t_res = resolve(t, 1)
v = ambient_space(t_res)
cox_ring(v)
M = matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(v))
MT = transpose(M)

@testset "FTheoryToolsPaper Section 4.1 Part 2" begin
  @test rank(MT) == 7
end

W = 2 * torusinvariant_prime_divisors(B3)[1]
t = literature_model(;
  arxiv_id="1109.3454",
  equation="3.1",
  base_space=B3,
  defining_classes=Dict("w" => W),
  rng=our_rng,
)
t_res = resolve(t, 1)

@testset "FTheoryToolsPaper Section 4.1 Part 3" begin
  @test is_partially_resolved(t_res)
  @test typeof(ambient_space(t_res)) == CoveredScheme{QQField}
end

#######################################
# 4.2 Toric Hypersurface Fibrations
#######################################

display_all_literature_models(Dict("gauge_algebra" => ["su(3)", "su(2)", "u(1)"]))
foah11_B3 = literature_model(;
  arxiv_id="1408.4808",
  equation="3.142",
  type="hypersurface",
  base_space=B3,
  defining_classes=Dict("s7" => Kbar, "s9" => Kbar),
  rng=our_rng,
)
foah11_B3 = literature_model(
  34; base_space=B3, defining_classes=Dict("s7" => Kbar, "s9" => Kbar), rng=our_rng
)
C = algebraic_closure(QQ)
diagonal_matrix(root_of_unity(C, 3), 3)

@testset "FTheoryToolsPaper Section 4.2 Part 1" begin
  @test model_description(foah11_B3) ==
    "F-theory hypersurface model with fiber ambient space F_11"
  @test is_zero(zero_section(foah11_B3)) == false
  @test length(generating_sections(foah11_B3)) == 1
  @test dim(gauge_algebra(foah11_B3)) == 12
  @test global_gauge_group_quotient(foah11_B3)[1][1] ==
    "diagonal_matrix(root_of_unity(C,3),3)"
end

foah11_B3_weier = weierstrass_model(foah11_B3)
R = cox_ring(base_space(foah11_B3_weier))
tuned_model = tune(foah11_B3_weier, Dict("s5" => zero(R)))

@testset "FTheoryToolsPaper Section 4.2 Part 2" begin
  @test length(birational_literature_models(foah11_B3)) == 1
  @test length(singular_loci(foah11_B3_weier)) == 3
  @test length(tunable_sections(foah11_B3_weier)) == 6
  @test length(singular_loci(tuned_model)) == 4
end

#######################################
# 5.1 Computation of G4-Fluxes
#######################################

qsm_model = literature_model(;
  arxiv_id="1903.00009", model_parameters=Dict("k" => 283), rng=our_rng
)
cox_ring(ambient_space(qsm_model))
cohomology_ring(ambient_space(qsm_model))
g4_gens = chosen_g4_flux_gens(qsm_model)

@testset "FTheoryToolsPaper Section 5.1 Part 1" begin
  @test ngens(stanley_reisner_ideal(ambient_space(qsm_model))) == 20
  @test betti_number(ambient_space(qsm_model), 4) == 25
  @test length(basis_of_h22_ambient(qsm_model; completeness_check=false)) == 25
  @test length(g4_gens) == 25
  @test parent(polynomial(cohomology_class(g4_gens[1]))) ==
    cohomology_ring(ambient_space(qsm_model))
end

fg = special_flux_family(qsm_model; completeness_check=false, rng=our_rng)
mat_int = matrix_integral(fg)
mat_rat = matrix_rational(fg)

@testset "FTheoryToolsPaper Section 5.1 Part 2" begin
  @test size(mat_int) == (25, 17)
  @test size(mat_rat) == (25, 0)
end

g4 = random_flux_instance(fg; rng=our_rng)
g4_2 = random_flux(qsm_model; completeness_check=false, rng=our_rng)
fg_not_breaking = special_flux_family(
  qsm_model; not_breaking=true, completeness_check=false, rng=our_rng
)

@testset "FTheoryToolsPaper Section 5.1 Part 3" begin
  @test size(matrix_integral(fg_not_breaking)) == (25, 1)
  @test size(matrix_rational(fg_not_breaking)) == (25, 0)
end

mat_int = matrix_integral(fg_not_breaking)
g4_sample = sum(mat_int[k] * g4_gens[k] for k in 1:length(g4_gens))
cohomology_class(g4_sample)

@testset "FTheoryToolsPaper Section 5.1 Part 4" begin
  @test breaks_non_abelian_gauge_group(g4_sample) == false
  @test kbar3(qsm_model) == 30
end

divs = torusinvariant_prime_divisors(ambient_space(qsm_model))
e1 = cohomology_class(divs[13])
e2 = cohomology_class(divs[10])
e4 = cohomology_class(divs[12])
u = cohomology_class(divs[11])
v = cohomology_class(divs[8])
pb_Kbar = cohomology_class(sum(divs[1:7]))
g4_class =
  (-3)//kbar3(qsm_model) *
  (5 * e1 * e4 + pb_Kbar * (-3 * e1 - 2 * e2 - 6 * e4 + pb_Kbar - 4 * u + v))
g4_sample2 = g4_flux(qsm_model, g4_class; completeness_check=false, consistency_check=false)
cohomology_class(g4_sample2)

@testset "FTheoryToolsPaper Section 5.1 Part 5" begin
  @test qsm_flux(qsm_model) == g4_sample2
  @test is_trivial(3 * cohomology_class(g4_sample) - cohomology_class(g4_sample2))
end

d3_tadpole_constraint(fg_not_breaking; rng=our_rng)
g4_exp = flux_instance(fg_not_breaking, [3], [])

@testset "FTheoryToolsPaper Section 5.1 Part 6" begin
  @test d3_tadpole_constraint(g4_exp) == 30
  @test d3_tadpole_constraint(g4_exp) == -1//12 * 3^2 + 123//4
end

#######################################
# 5.2 G4-Fluxes in a Complex Model
#######################################

t = literature_model(; arxiv_id="1511.03209", rng=our_rng)
t_res = resolve(t, 1)
amb = ambient_space(t_res)
cohomology_ring(amb; completeness_check=false)
g4_amb_candidates = chosen_g4_flux_gens(t_res; completeness_check=false)
fg_quant = special_flux_family(t_res; completeness_check=false, rng=our_rng)
fg_quant_no_break = special_flux_family(
  t_res; not_breaking=true, completeness_check=false, rng=our_rng
)

@testset "FTheoryToolsPaper Section 5.2" begin
  @test betti_number(amb, 4) == 1109
  @test length(g4_amb_candidates) == 629
  size(matrix_integral(fg_quant)) == (629, 224)
  size(matrix_rational(fg_quant)) == (629, 127)
  size(matrix_integral(fg_quant_no_break)) == (629, 1)
  size(matrix_rational(fg_quant_no_break)) == (629, 127)
end
