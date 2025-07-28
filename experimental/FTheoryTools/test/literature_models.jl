#############################################################
# 1: Literature Tate model over concrete base
#############################################################

using Random
our_rng = Random.Xoshiro(1234)

B3 = projective_space(NormalToricVariety, 3)
Kbar = anticanonical_divisor_class(B3)
w = torusinvariant_prime_divisors(B3)[1]
t1 = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, model_sections = Dict("w" => w), completeness_check = false)
D = classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes(t1)

@testset "Test defining data for literature Tate model over concrete base" begin
  @test parent(tate_section_a1(t1)) == cox_ring(base_space(t1))
  @test parent(tate_section_a2(t1)) == cox_ring(base_space(t1))
  @test parent(tate_section_a3(t1)) == cox_ring(base_space(t1))
  @test parent(tate_section_a4(t1)) == cox_ring(base_space(t1))
  @test parent(tate_section_a6(t1)) == cox_ring(base_space(t1))
  @test parent(tate_polynomial(t1)) == cox_ring(ambient_space(t1))
  @test parent(discriminant(t1)) == cox_ring(base_space(t1))
  @test length(singular_loci(t1; rng = our_rng)) == 2
  @test dim(base_space(t1)) == 3
  @test dim(ambient_space(t1)) == 5
  @test is_base_space_fully_specified(t1) == true
  @test is_base_space_fully_specified(t1) == is_base_space_fully_specified(weierstrass_model(t1))
  @test is_smooth(ambient_space(t1)) == false
  @test toric_variety(calabi_yau_hypersurface(t1)) == ambient_space(t1)
  @test sum(D["a43"].*[Kbar, toric_divisor_class(w)]) == classes_of_model_sections(t1)["a43"]
end

@testset "Test meta data for literature Tate model over concrete base" begin
  @test arxiv_id(t1) == "1109.3454"
  @test arxiv_doi(t1) == "10.48550/arXiv.1109.3454"
  @test arxiv_link(t1) == "https://arxiv.org/abs/1109.3454v2"
  @test arxiv_model_equation_number(t1) ==  "3.1"
  @test arxiv_model_page(t1) == "10"
  @test arxiv_model_section(t1) == "3"
  @test arxiv_version(t1) ==  "2"
  @test length(generating_sections(t1)) == 1
  @test journal_doi(t1) == "10.1016/j.nuclphysb.2011.12.013"
  @test journal_link(t1) == "https://www.sciencedirect.com/science/article/pii/S0550321311007115"
  @test journal_model_equation_number(t1) == "3.1"
  @test journal_model_page(t1) == "9"
  @test journal_model_section(t1) == "3"
  @test journal_name(t1) == "Nucl. Phys. B"
  @test journal_pages(t1) == "1–47"
  @test journal_volume(t1) == "858"
  @test journal_year(t1) ==  "2012"
  @test literature_identifier(t1) ==  "1109_3454"
  @test model_description(t1) == "SU(5)xU(1) restricted Tate model"
  @test paper_authors(t1) == ["Sven Krause", "Christoph Mayrhofer", "Timo Weigand"]
  @test paper_buzzwords(t1) == ["GUT model", "Tate", "U(1)", "SU(5)"]
  @test paper_description(t1) == "SU(5)xU(1) restricted Tate model"
  @test paper_title(t1) == "\$G_4\$ flux, chiral matter and singularity resolution in F-theory compactifications"
  @test resolutions(t1) == [([["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]], ["e1", "e4", "e2", "e3", "s"])]
  @test length(resolution_generating_sections(t1)) == 1
  @test length(resolution_zero_sections(t1)) == 1
  @test weighted_resolutions(t1) == [([(["x", "y", "w"], [1, 1, 1]), (["x", "y", "w"], [1, 2, 1]), (["x", "y", "w"], [2, 2, 1]), (["x", "y", "w"], [2, 3, 1]), (["x", "y"], [1, 1])], ["e1", "e4", "e2", "e3", "s"])]
  @test length(weighted_resolution_generating_sections(t1)) == 1
  @test length(weighted_resolution_zero_sections(t1)) == 1
end

@testset "Test error messages for literature Tate model over concrete base" begin
  @test_throws ArgumentError associated_literature_models(t1)
  @test_throws ArgumentError journal_report_numbers(t1)
  @test_throws ArgumentError model_parameters(t1)
  @test_throws ArgumentError birational_literature_models(t1)
end

t2 = resolve(t1, 1)

@testset "Test resolving literature Tate model over concrete base" begin
  @test is_smooth(ambient_space(t2)) == false
  @test is_partially_resolved(t2) == true
  @test base_space(t1) == base_space(t2)
  @test length(exceptional_divisor_indices(t2)) == length(exceptional_classes(t2))
  @test length(exceptional_divisor_indices(t2)) == length(exceptional_divisor_indices(t1)) + 5
  @test length(exceptional_classes(t2)) == length(exceptional_classes(t1)) + 5
end

add_resolution(t1, [["x", "y"], ["y", "s", "w"], ["s", "e4"], ["s", "e3"], ["s", "e1"]], ["s", "w", "e3", "e1", "e2"])

@testset "Test adding new resolution to literature Tate model over concrete base" begin
  @test length(resolutions(t1)) == 2
end



#############################################################
# 2: Literature Weierstrass model over concrete base
#############################################################

B2 = projective_space(NormalToricVariety, 2)
b = torusinvariant_prime_divisors(B2)[1]
w1 = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, defining_classes = Dict("b" => b), completeness_check = false)

@testset "Test defining data for literature Weierstrass model over concrete base" begin
  @test parent(weierstrass_section_f(w1)) == cox_ring(base_space(w1))
  @test parent(weierstrass_section_g(w1)) == cox_ring(base_space(w1))
  @test parent(weierstrass_polynomial(w1)) == cox_ring(ambient_space(w1))
  @test parent(discriminant(w1)) == cox_ring(base_space(w1))
  @test length(singular_loci(w1; rng = our_rng)) == 1
  @test dim(base_space(w1)) == 2
  @test dim(ambient_space(w1)) == 4
  @test is_base_space_fully_specified(w1) == true
  @test is_smooth(ambient_space(w1)) == false
  @test toric_variety(calabi_yau_hypersurface(w1)) == ambient_space(w1)
end

@testset "Test meta data for literature Weierstrass model over concrete base" begin
  @test arxiv_id(w1) == "1208.2695"
  @test arxiv_doi(w1) == "10.48550/arXiv.1208.2695"
  @test arxiv_link(w1) == "https://arxiv.org/abs/1208.2695v2"
  @test arxiv_model_equation_number(w1) ==  "B.19"
  @test arxiv_model_page(w1) == "34"
  @test arxiv_model_section(w1) == "B"
  @test arxiv_version(w1) == "2"
  @test birational_literature_models(w1) == ["1208_2695-1"]
  @test length(generating_sections(w1)) == 1
  @test journal_doi(w1) == "10.1007/JHEP10(2012)128"
  @test journal_link(w1) == "https://link.springer.com/article/10.1007/JHEP10(2012)128"
  @test journal_model_equation_number(w1) == "B.19"
  @test journal_model_page(w1) == "34"
  @test journal_model_section(w1) == "B"
  # FIXME: This should likely be a range? But then, maybe we do not care?
  @test journal_pages(w1) == "128"
  @test journal_report_numbers(w1) == ["UCSB-MATH-2012-27", "MIT-CTP-4388"]
  @test journal_volume(w1) == "10"
  @test journal_year(w1) == "2012"
  @test literature_identifier(w1) == "1208_2695-2"
  @test model_description(w1) == "U(1) Weierstrass model"
  @test paper_authors(w1) == ["David R. Morrison", "Daniel S. Park"]
  @test paper_buzzwords(w1) == ["U(1)", "Mordell–Weil", "rational sections"]
  @test paper_description(w1) == "General construction of U(1) model"
  @test paper_title(w1) ==  "F-Theory and the Mordell-Weil Group of Elliptically-Fibered Calabi-Yau Threefolds"
end

@testset "Test error messages for literature Weierstrass model over concrete base" begin
  @test_throws ArgumentError model_parameters(w1)
  @test_throws ArgumentError associated_literature_models(w1)
  @test_throws ArgumentError resolutions(w1)
  @test_throws ArgumentError resolution_generating_sections(w1)
  @test_throws ArgumentError resolution_zero_sections(w1)
  @test_throws ArgumentError weighted_resolutions(w1)
  @test_throws ArgumentError weighted_resolution_generating_sections(w1)
  @test_throws ArgumentError weighted_resolution_zero_sections(w1)
end



#############################################################
# 3: Literature Tate model over arbitrary base
#############################################################

t3 = literature_model(arxiv_id = "1109.3454", equation = "3.1")

@testset "Basic tests for literature Tate model over arbitrary base" begin
  @test parent(tate_section_a1(t3)) == cox_ring(base_space(t3))
  @test parent(tate_section_a2(t3)) == cox_ring(base_space(t3))
  @test parent(tate_section_a3(t3)) == cox_ring(base_space(t3))
  @test parent(tate_section_a4(t3)) == cox_ring(base_space(t3))
  @test parent(tate_section_a6(t3)) == cox_ring(base_space(t3))
  @test parent(tate_polynomial(t3)) == cox_ring(ambient_space(t3))
  @test parent(discriminant(t3)) == cox_ring(base_space(t3))
  @test length(singular_loci(t3; rng = our_rng)) == 2
  @test dim(base_space(t3)) == 3
  @test dim(ambient_space(t3)) == 5
  @test is_base_space_fully_specified(t3) == false
  @test is_base_space_fully_specified(t3) == is_base_space_fully_specified(weierstrass_model(t3))
end

@testset "Test meta data for literature Tate model over arbitrary base" begin
  @test arxiv_id(t3) == "1109.3454"
  @test arxiv_doi(t3) == "10.48550/arXiv.1109.3454"
  @test arxiv_link(t3) == "https://arxiv.org/abs/1109.3454v2"
  @test arxiv_model_equation_number(t3) ==  "3.1"
  @test arxiv_model_page(t3) == "10"
  @test arxiv_model_section(t3) == "3"
  @test arxiv_version(t3) ==  "2"
  @test length(generating_sections(t3)) == 1
  @test journal_doi(t3) == "10.1016/j.nuclphysb.2011.12.013"
  @test journal_link(t3) == "https://www.sciencedirect.com/science/article/pii/S0550321311007115"
  @test journal_model_equation_number(t3) == "3.1"
  @test journal_model_page(t3) == "9"
  @test journal_model_section(t3) == "3"
  @test journal_pages(t3) == "1–47"
  @test journal_volume(t3) == "858"
  @test journal_year(t3) ==  "2012"
  @test literature_identifier(t3) ==  "1109_3454"
  @test model_description(t3) == "SU(5)xU(1) restricted Tate model"
  @test paper_authors(t3) == ["Sven Krause", "Christoph Mayrhofer", "Timo Weigand"]
  @test paper_buzzwords(t3) == ["GUT model", "Tate", "U(1)", "SU(5)"]
  @test paper_description(t3) == "SU(5)xU(1) restricted Tate model"
  @test paper_title(t3) == "\$G_4\$ flux, chiral matter and singularity resolution in F-theory compactifications"
  @test resolutions(t3) == [([["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]], ["e1", "e4", "e2", "e3", "s"])]
  @test length(resolution_generating_sections(t3)) == 1
  @test length(resolution_zero_sections(t3)) == 1
  @test weighted_resolutions(t3) == [([(["x", "y", "w"], [1, 1, 1]), (["x", "y", "w"], [1, 2, 1]), (["x", "y", "w"], [2, 2, 1]), (["x", "y", "w"], [2, 3, 1]), (["x", "y"], [1, 1])], ["e1", "e4", "e2", "e3", "s"])]
  @test length(weighted_resolution_generating_sections(t3)) == 1
  @test length(weighted_resolution_zero_sections(t3)) == 1
end

@testset "Test error messages for literature Tate model over arbitrary base" begin
  @test_throws ArgumentError associated_literature_models(t3)
  @test_throws ArgumentError journal_report_numbers(t3)
  @test_throws ArgumentError model_parameters(t3)
  @test_throws ArgumentError birational_literature_models(t3)
  @test_throws ArgumentError literature_model(arxiv_id = "1212.2949", equation = "3.2")
end


t4 = literature_model(arxiv_id = "1212.2949", equation = "3.2", model_parameters = Dict("k" => 5))
t5 = literature_model(arxiv_id = "1212.2949", equation = "3.42", model_parameters = Dict("k" => 4))
t6 = literature_model(arxiv_id = "1212.2949", equation = "4.1", model_parameters = Dict("k" => 3))
t7 = literature_model(arxiv_id = "1212.2949", equation = "4.23", model_parameters = Dict("k" => 4))
t8  = literature_model(arxiv_id = "1212.2949", equation = "5.1")
t9 = literature_model(arxiv_id = "1212.2949", equation = "5.7")
t10 = literature_model(arxiv_id = "1212.2949", equation = "5.13")

@testset "Test more Tate models (from paper Tate form on Steroids) in our database, by constructing them over arbitrary bases" begin
  @test is_base_space_fully_specified(t4) == false
  @test is_partially_resolved(t5) == false
  @test length(resolutions(t6)) == 1
  @test length(paper_authors(t7)) == 2
  @test model_description(t8) == "E6 Tate model"
  @test model_description(t9) == "E7 Tate model"
  @test model_description(t10) == "E8 Tate model"
end


#############################################################
# 4: Literature Tate model over concrete base
#############################################################

ζ0 = torusinvariant_prime_divisors(B3)[1]
t4b = literature_model(arxiv_id = "1212.2949", equation = "3.2", model_parameters = Dict("k" => 2), base_space = B3, model_sections = Dict("ζ0" => toric_line_bundle(ζ0)))
t5b = literature_model(arxiv_id = "1212.2949", equation = "3.42", model_parameters = Dict("k" => 2), base_space = B3, model_sections = Dict("ζ0" => toric_divisor_class(ζ0)))
t6b = literature_model(arxiv_id = "1212.2949", equation = "4.1", model_parameters = Dict("k" => 2), base_space = B3, defining_classes = Dict("ζ0" => toric_line_bundle(ζ0)))
t7b = literature_model(arxiv_id = "1212.2949", equation = "4.23", model_parameters = Dict("k" => 2), base_space = B3, defining_classes = Dict("ζ0" => toric_divisor_class(ζ0)))
t8b  = literature_model(arxiv_id = "1212.2949", equation = "5.1", base_space = B3, defining_classes = Dict("ζ0" => ζ0))
t9b = literature_model(arxiv_id = "1212.2949", equation = "5.7", base_space = B3, defining_classes = Dict("ζ0" => ζ0))
t10b = literature_model(arxiv_id = "1212.2949", equation = "5.13", base_space = B3, defining_classes = Dict("ζ0" => ζ0))

@testset "Test more Tate models (from paper Tate form on Steroids) in our database, by constructing them over concrete bases" begin
  @test is_base_space_fully_specified(t4b) == true
  @test is_partially_resolved(t5b) == false
  @test length(resolutions(t6b)) == 1
  @test length(paper_authors(t7b)) == 2
  @test model_description(t8b) == "E6 Tate model"
  @test model_description(t9b) == "E7 Tate model"
  @test model_description(t10b) == "E8 Tate model"
end


#############################################################
# 5: Literature Weierstrass model over arbitrary base
#############################################################

w2 = literature_model(arxiv_id = "1208.2695", equation = "B.19", completeness_check = false)

@testset "Test defining data for literature Weierstrass model over arbitrary base" begin
  @test parent(weierstrass_section_f(w2)) == cox_ring(base_space(w2))
  @test parent(weierstrass_section_g(w2)) == cox_ring(base_space(w2))
  @test parent(weierstrass_polynomial(w2)) == cox_ring(ambient_space(w2))
  @test parent(discriminant(w2)) == cox_ring(base_space(w2))
  @test length(singular_loci(w2; rng = our_rng)) == 1
  @test dim(base_space(w2)) == 2
  @test dim(ambient_space(w2)) == 4
  @test is_base_space_fully_specified(w2) == false
end

@testset "Test meta data for literature Weierstrass model over arbitrary base" begin
  @test arxiv_id(w2) == "1208.2695"
  @test arxiv_doi(w2) == "10.48550/arXiv.1208.2695"
  @test arxiv_link(w2) == "https://arxiv.org/abs/1208.2695v2"
  @test arxiv_model_equation_number(w2) ==  "B.19"
  @test arxiv_model_page(w2) == "34"
  @test arxiv_model_section(w2) == "B"
  @test arxiv_version(w2) == "2"
  @test birational_literature_models(w2) == ["1208_2695-1"]
  @test length(generating_sections(w2)) == 1
  @test journal_doi(w2) == "10.1007/JHEP10(2012)128"
  @test journal_link(w2) == "https://link.springer.com/article/10.1007/JHEP10(2012)128"
  @test journal_model_equation_number(w2) == "B.19"
  @test journal_model_page(w2) == "34"
  @test journal_model_section(w2) == "B"
  # FIXME: This should likely be a range? But then, maybe we do not care?
  @test journal_pages(w2) == "128"
  @test journal_report_numbers(w2) == ["UCSB-MATH-2012-27", "MIT-CTP-4388"]
  @test journal_volume(w2) == "10"
  @test journal_year(w2) == "2012"
  @test literature_identifier(w2) == "1208_2695-2"
  @test model_description(w2) == "U(1) Weierstrass model"
  @test paper_authors(w2) == ["David R. Morrison", "Daniel S. Park"]
  @test paper_buzzwords(w2) == ["U(1)", "Mordell–Weil", "rational sections"]
  @test paper_description(w2) == "General construction of U(1) model"
  @test paper_title(w2) ==  "F-Theory and the Mordell-Weil Group of Elliptically-Fibered Calabi-Yau Threefolds"
end

@testset "Test error messages for literature Weierstrass model over arbitrary base" begin
  @test_throws ArgumentError model_parameters(w2)
  @test_throws ArgumentError associated_literature_models(w2)
  @test_throws ArgumentError resolutions(w2)
  @test_throws ArgumentError resolution_generating_sections(w2)
  @test_throws ArgumentError resolution_zero_sections(w2)
  @test_throws ArgumentError weighted_resolutions(w2)
  @test_throws ArgumentError weighted_resolution_generating_sections(w2)
  @test_throws ArgumentError weighted_resolution_zero_sections(w2)
end


w3 = literature_model(arxiv_id = "1507.05954", equation = "A.1")


@testset "Test more Weierstrass literature models in our database over arbitrary base" begin
  @test model_description(w3) == "U(1)xU(1) Weierstrass model"
end


#############################################################
# 6: Literature Weierstrass model over arbitrary base
#############################################################

B2 = projective_space(NormalToricVariety, 2)
b = torusinvariant_prime_divisors(B2)[2]
w4 = literature_model(arxiv_id = "1208.2695", equation = "B.19", completeness_check = false, base_space = B2, defining_classes = Dict("b" =>b))
b2 = anticanonical_divisor(B2)
w5 = literature_model(arxiv_id = "1507.05954", equation = "A.1", completeness_check = false, base_space = B2, defining_classes = Dict("s8" => b2, "a1" => b, "a2" => b, "a3" => b))

@testset "Test defining data for literature Weierstrass model over concrete base" begin
  @test length(singular_loci(w4; rng = our_rng)) == 1
  @test dim(base_space(w4)) == 2
  @test dim(ambient_space(w4)) == 4
  @test is_base_space_fully_specified(w4) == true
  @test model_description(w5) == "U(1)xU(1) Weierstrass model"
end



#############################################################
# 7: Test other literature model constructor
#############################################################

B2 = projective_space(NormalToricVariety, 2)
b = torusinvariant_prime_divisors(B2)[1]
w6 = literature_model(3, base_space = B2, defining_classes = Dict("b" => b), completeness_check = false)

@testset "Test defining data for literature model defined by model index" begin
  @test length(singular_loci(w6; rng = our_rng)) == 1
  @test dim(base_space(w6)) == 2
  @test dim(ambient_space(w6)) == 4
  @test is_base_space_fully_specified(w6) == true
  @test model_description(w6) == "U(1) Weierstrass model"
end

@testset "Test error messages for literature Weierstrass model over arbitrary base" begin
  @test_throws ArgumentError literature_model(-1)
  @test_throws ArgumentError literature_model(0)
  @test_throws ArgumentError literature_model(205)
end

h = literature_model(arxiv_id = "1507.05954", equation = "3.4")

@testset "Test for literature hypersurface model over arbitrary base" begin
  @test parent(hypersurface_equation(h)) == coordinate_ring(ambient_space(h))
  @test dim(base_space(h)) == 2
  @test is_smooth(fiber_ambient_space(h)) == true
  @test symbols(cox_ring(fiber_ambient_space(h))) == [:u, :v, :w]
  @test is_base_space_fully_specified(h) == false
  @test is_partially_resolved(h) == false
  @test string.(zero_section(h)) == ["0", "-b1", "a1"]
end



#############################################################################
# 8: Test models from F-theory on all toric hypersurfaces over arbitrary base
#############################################################################

foah1 = literature_model(arxiv_id = "1408.4808", equation = "3.4", type = "hypersurface")
foah2 = literature_model(arxiv_id = "1408.4808", equation = "3.12", type = "hypersurface")
foah3 = literature_model(arxiv_id = "1408.4808", equation = "3.54", type = "hypersurface")
foah4 = literature_model(arxiv_id = "1408.4808", equation = "3.17", type = "hypersurface")
foah5 = literature_model(arxiv_id = "1408.4808", equation = "3.73", type = "hypersurface")
foah6 = literature_model(arxiv_id = "1408.4808", equation = "3.82", type = "hypersurface")
foah7 = literature_model(arxiv_id = "1408.4808", equation = "3.96", type = "hypersurface")
foah8 = literature_model(arxiv_id = "1408.4808", equation = "3.106", type = "hypersurface")
foah9 = literature_model(arxiv_id = "1408.4808", equation = "3.118", type = "hypersurface")
foah10 = literature_model(arxiv_id = "1408.4808", equation = "3.130", type = "hypersurface")
foah11 = literature_model(arxiv_id = "1408.4808", equation = "3.142", type = "hypersurface")
foah12 = literature_model(arxiv_id = "1408.4808", equation = "3.155", type = "hypersurface")
foah13 = literature_model(arxiv_id = "1408.4808", equation = "3.181", type = "hypersurface")
foah14 = literature_model(arxiv_id = "1408.4808", equation = "3.168", type = "hypersurface")
foah15 = literature_model(arxiv_id = "1408.4808", equation = "3.190", type = "hypersurface")
foah16 = literature_model(arxiv_id = "1408.4808", equation = "3.203", type = "hypersurface")

@testset "Test hypersurface form of models in F-theory on all toric hypersurfaces, defined over arbitrary base" begin
  @test dim(base_space(foah1)) == 3
  @test dim(base_space(foah2)) == 3
  @test dim(base_space(foah3)) == 3
  @test dim(base_space(foah4)) == 3
  @test dim(base_space(foah5)) == 3
  @test dim(base_space(foah6)) == 3
  @test dim(base_space(foah7)) == 3
  @test dim(base_space(foah8)) == 3
  @test dim(base_space(foah9)) == 3
  @test dim(base_space(foah10)) == 3
  @test dim(base_space(foah11)) == 3
  @test dim(base_space(foah12)) == 3
  @test dim(base_space(foah13)) == 3
  @test dim(base_space(foah14)) == 3
  @test dim(base_space(foah15)) == 3
  @test dim(base_space(foah16)) == 3
  @test dim(ambient_space(foah1)) == 5
  @test dim(ambient_space(foah2)) == 5
  @test dim(ambient_space(foah3)) == 5
  @test dim(ambient_space(foah4)) == 5
  @test dim(ambient_space(foah5)) == 5
  @test dim(ambient_space(foah6)) == 5
  @test dim(ambient_space(foah7)) == 5
  @test dim(ambient_space(foah8)) == 5
  @test dim(ambient_space(foah9)) == 5
  @test dim(ambient_space(foah10)) == 5
  @test dim(ambient_space(foah11)) == 5
  @test dim(ambient_space(foah12)) == 5
  @test dim(ambient_space(foah13)) == 5
  @test dim(ambient_space(foah14)) == 5
  @test dim(ambient_space(foah15)) == 5
  @test dim(ambient_space(foah16)) == 5
  @test is_base_space_fully_specified(foah1) == false
  @test is_base_space_fully_specified(foah2) == false
  @test is_base_space_fully_specified(foah3) == false
  @test is_base_space_fully_specified(foah4) == false
  @test is_base_space_fully_specified(foah5) == false
  @test is_base_space_fully_specified(foah6) == false
  @test is_base_space_fully_specified(foah7) == false
  @test is_base_space_fully_specified(foah8) == false
  @test is_base_space_fully_specified(foah9) == false
  @test is_base_space_fully_specified(foah10) == false
  @test is_base_space_fully_specified(foah11) == false
  @test is_base_space_fully_specified(foah12) == false
  @test is_base_space_fully_specified(foah13) == false
  @test is_base_space_fully_specified(foah14) == false
  @test is_base_space_fully_specified(foah15) == false
  @test is_base_space_fully_specified(foah16) == false
  @test model_description(foah1) == "F-theory hypersurface model with fiber ambient space F_1"
  @test model_description(foah2) == "F-theory hypersurface model with fiber ambient space F_2"
  @test model_description(foah3) == "F-theory hypersurface model with fiber ambient space F_3"
  @test model_description(foah4) == "F-theory hypersurface model with fiber ambient space F_4"
  @test model_description(foah5) == "F-theory hypersurface model with fiber ambient space F_5"
  @test model_description(foah6) == "F-theory hypersurface model with fiber ambient space F_6"
  @test model_description(foah7) == "F-theory hypersurface model with fiber ambient space F_7"
  @test model_description(foah8) == "F-theory hypersurface model with fiber ambient space F_8"
  @test model_description(foah9) == "F-theory hypersurface model with fiber ambient space F_9"
  @test model_description(foah10) == "F-theory hypersurface model with fiber ambient space F_10"
  @test model_description(foah11) == "F-theory hypersurface model with fiber ambient space F_11"
  @test model_description(foah12) == "F-theory hypersurface model with fiber ambient space F_12"
  @test model_description(foah13) == "F-theory hypersurface model with fiber ambient space F_13"
  @test model_description(foah14) == "F-theory hypersurface model with fiber ambient space F_14"
  @test model_description(foah15) == "F-theory hypersurface model with fiber ambient space F_15"
  @test model_description(foah16) == "F-theory hypersurface model with fiber ambient space F_16"
  @test haskey(explicit_model_sections(foah6), "s9") == false
  @test dim(gauge_algebra(foah6)) == 4
  @test length(global_gauge_group_quotient(foah6)) == 2
  @test dim(gauge_algebra(foah8)) == 7
  @test length(global_gauge_group_quotient(foah8)) == 3
  @test dim(gauge_algebra(foah9)) == 5
  @test length(global_gauge_group_quotient(foah9)) == 3
  @test dim(gauge_algebra(foah11)) == 12
  @test length(global_gauge_group_quotient(foah11)) == 3
  @test dim(gauge_algebra(foah12)) == 8
  @test length(global_gauge_group_quotient(foah12)) == 4
  @test dim(gauge_algebra(foah13)) == 21
  @test length(global_gauge_group_quotient(foah13)) == 3
  @test dim(gauge_algebra(foah14)) == 15
  @test length(global_gauge_group_quotient(foah14)) == 4
  @test dim(gauge_algebra(foah15)) == 13
  @test length(global_gauge_group_quotient(foah15)) == 5
  @test dim(gauge_algebra(foah16)) == 24
  @test length(global_gauge_group_quotient(foah16)) == 3
end



#############################################################################
# 9: Test models from F-theory on all toric hypersurfaces over concrete base
#############################################################################

B3 = projective_space(NormalToricVariety, 3)
Kbar = anticanonical_divisor(B3)
foah1_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.4", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah2_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.12", type = "hypersurface", base_space = B3, defining_classes = Dict("b7" => Kbar, "b9" => Kbar), completeness_check = false)
foah3_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.54", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah4_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.17", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah5_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.73", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah6_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.82", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah7_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.96", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah8_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.106", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah9_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.118", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah10_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.130", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah11_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.142", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah12_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.155", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah13_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.181", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah14_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.168", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah15_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.190", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)
foah16_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.203", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar), completeness_check = false)

@testset "Test hypersurface form of models in F-theory on all toric hypersurfaces, defined over concrete base" begin
  @test dim(base_space(foah1_B3)) == 3
  @test dim(base_space(foah2_B3)) == 3
  @test dim(base_space(foah3_B3)) == 3
  @test dim(base_space(foah4_B3)) == 3
  @test dim(base_space(foah5_B3)) == 3
  @test dim(base_space(foah6_B3)) == 3
  @test dim(base_space(foah7_B3)) == 3
  @test dim(base_space(foah8_B3)) == 3
  @test dim(base_space(foah9_B3)) == 3
  @test dim(base_space(foah10_B3)) == 3
  @test dim(base_space(foah11_B3)) == 3
  @test dim(base_space(foah12_B3)) == 3
  @test dim(base_space(foah13_B3)) == 3
  @test dim(base_space(foah14_B3)) == 3
  @test dim(base_space(foah15_B3)) == 3
  @test dim(base_space(foah16_B3)) == 3
  @test dim(ambient_space(foah1_B3)) == 5
  @test dim(ambient_space(foah2_B3)) == 5
  @test dim(ambient_space(foah3_B3)) == 5
  @test dim(ambient_space(foah4_B3)) == 5
  @test dim(ambient_space(foah5_B3)) == 5
  @test dim(ambient_space(foah6_B3)) == 5
  @test dim(ambient_space(foah7_B3)) == 5
  @test dim(ambient_space(foah8_B3)) == 5
  @test dim(ambient_space(foah9_B3)) == 5
  @test dim(ambient_space(foah10_B3)) == 5
  @test dim(ambient_space(foah11_B3)) == 5
  @test dim(ambient_space(foah12_B3)) == 5
  @test dim(ambient_space(foah13_B3)) == 5
  @test dim(ambient_space(foah14_B3)) == 5
  @test dim(ambient_space(foah15_B3)) == 5
  @test dim(ambient_space(foah16_B3)) == 5
  @test is_base_space_fully_specified(foah1_B3) == true
  @test is_base_space_fully_specified(foah2_B3) == true
  @test is_base_space_fully_specified(foah3_B3) == true
  @test is_base_space_fully_specified(foah4_B3) == true
  @test is_base_space_fully_specified(foah5_B3) == true
  @test is_base_space_fully_specified(foah6_B3) == true
  @test is_base_space_fully_specified(foah7_B3) == true
  @test is_base_space_fully_specified(foah8_B3) == true
  @test is_base_space_fully_specified(foah9_B3) == true
  @test is_base_space_fully_specified(foah10_B3) == true
  @test is_base_space_fully_specified(foah11_B3) == true
  @test is_base_space_fully_specified(foah12_B3) == true
  @test is_base_space_fully_specified(foah13_B3) == true
  @test is_base_space_fully_specified(foah14_B3) == true
  @test is_base_space_fully_specified(foah15_B3) == true
  @test is_base_space_fully_specified(foah16_B3) == true
  @test model_description(foah1_B3) == "F-theory hypersurface model with fiber ambient space F_1"
  @test model_description(foah2_B3) == "F-theory hypersurface model with fiber ambient space F_2"
  @test model_description(foah3_B3) == "F-theory hypersurface model with fiber ambient space F_3"
  @test model_description(foah4_B3) == "F-theory hypersurface model with fiber ambient space F_4"
  @test model_description(foah5_B3) == "F-theory hypersurface model with fiber ambient space F_5"
  @test model_description(foah6_B3) == "F-theory hypersurface model with fiber ambient space F_6"
  @test model_description(foah7_B3) == "F-theory hypersurface model with fiber ambient space F_7"
  @test model_description(foah8_B3) == "F-theory hypersurface model with fiber ambient space F_8"
  @test model_description(foah9_B3) == "F-theory hypersurface model with fiber ambient space F_9"
  @test model_description(foah10_B3) == "F-theory hypersurface model with fiber ambient space F_10"
  @test model_description(foah11_B3) == "F-theory hypersurface model with fiber ambient space F_11"
  @test model_description(foah12_B3) == "F-theory hypersurface model with fiber ambient space F_12"
  @test model_description(foah13_B3) == "F-theory hypersurface model with fiber ambient space F_13"
  @test model_description(foah14_B3) == "F-theory hypersurface model with fiber ambient space F_14"
  @test model_description(foah15_B3) == "F-theory hypersurface model with fiber ambient space F_15"
  @test model_description(foah16_B3) == "F-theory hypersurface model with fiber ambient space F_16"
  @test parent(explicit_model_sections(foah1_B3)["s7"]) == cox_ring(base_space(foah1_B3))
  @test parent(explicit_model_sections(foah2_B3)["b7"]) == cox_ring(base_space(foah2_B3))
  @test parent(explicit_model_sections(foah3_B3)["s7"]) == cox_ring(base_space(foah3_B3))
  @test parent(explicit_model_sections(foah4_B3)["d4"]) == cox_ring(base_space(foah4_B3))
  @test parent(explicit_model_sections(foah5_B3)["s7"]) == cox_ring(base_space(foah5_B3))
  @test parent(explicit_model_sections(foah6_B3)["s7"]) == cox_ring(base_space(foah6_B3))
  @test parent(explicit_model_sections(foah7_B3)["s7"]) == cox_ring(base_space(foah7_B3))
  @test parent(explicit_model_sections(foah8_B3)["s7"]) == cox_ring(base_space(foah8_B3))
  @test parent(explicit_model_sections(foah9_B3)["s7"]) == cox_ring(base_space(foah9_B3))
  @test parent(explicit_model_sections(foah10_B3)["s5"]) == cox_ring(base_space(foah10_B3))
  @test parent(explicit_model_sections(foah11_B3)["s5"]) == cox_ring(base_space(foah11_B3))
  @test parent(explicit_model_sections(foah12_B3)["s7"]) == cox_ring(base_space(foah12_B3))
  @test parent(explicit_model_sections(foah13_B3)["s1"]) == cox_ring(base_space(foah13_B3))
  @test parent(explicit_model_sections(foah14_B3)["s7"]) == cox_ring(base_space(foah14_B3))
  @test parent(explicit_model_sections(foah15_B3)["s7"]) == cox_ring(base_space(foah15_B3))
  @test parent(explicit_model_sections(foah16_B3)["s7"]) == cox_ring(base_space(foah16_B3))
  @test string(hypersurface_equation_parametrization(foah1_B3)) == "s1*u^3 + s2*u^2*v + s3*u*v^2 + s4*v^3 + s5*u^2*w + s6*u*v*w + s7*v^2*w + s8*u*w^2 + s9*v*w^2 + s10*w^3"
end



##########################################################################################################
# 10: Test Weierstrass counterparts of models from F-theory on all toric hypersurfaces over arbitrary base
##########################################################################################################

foah1_weier = literature_model(arxiv_id = "1408.4808", equation = "3.4", type = "weierstrass")
foah2_weier = literature_model(arxiv_id = "1408.4808", equation = "3.12", type = "weierstrass")
foah3_weier = literature_model(arxiv_id = "1408.4808", equation = "3.54", type = "weierstrass")
foah4_weier = literature_model(arxiv_id = "1408.4808", equation = "3.17", type = "weierstrass")
foah5_weier = literature_model(arxiv_id = "1408.4808", equation = "3.73", type = "weierstrass")
foah6_weier = literature_model(arxiv_id = "1408.4808", equation = "3.82", type = "weierstrass")
foah7_weier = literature_model(arxiv_id = "1408.4808", equation = "3.96", type = "weierstrass")
foah8_weier = literature_model(arxiv_id = "1408.4808", equation = "3.106", type = "weierstrass")
foah9_weier = literature_model(arxiv_id = "1408.4808", equation = "3.118", type = "weierstrass")
foah10_weier = literature_model(arxiv_id = "1408.4808", equation = "3.130", type = "weierstrass")
foah11_weier = literature_model(arxiv_id = "1408.4808", equation = "3.142", type = "weierstrass")
foah12_weier = literature_model(arxiv_id = "1408.4808", equation = "3.155", type = "weierstrass")
foah13_weier = literature_model(arxiv_id = "1408.4808", equation = "3.181", type = "weierstrass")
foah14_weier = literature_model(arxiv_id = "1408.4808", equation = "3.168", type = "weierstrass")
foah15_weier = literature_model(arxiv_id = "1408.4808", equation = "3.190", type = "weierstrass")
foah16_weier = weierstrass_model(foah16)

@testset "Test Weierstrass form of models in F-theory on all toric hypersurfaces, defined over arbitrary base" begin
  @test dim(base_space(foah1_weier)) == 3
  @test dim(base_space(foah2_weier)) == 3
  @test dim(base_space(foah3_weier)) == 3
  @test dim(base_space(foah4_weier)) == 3
  @test dim(base_space(foah5_weier)) == 3
  @test dim(base_space(foah6_weier)) == 3
  @test dim(base_space(foah7_weier)) == 3
  @test dim(base_space(foah8_weier)) == 3
  @test dim(base_space(foah9_weier)) == 3
  @test dim(base_space(foah10_weier)) == 3
  @test dim(base_space(foah11_weier)) == 3
  @test dim(base_space(foah12_weier)) == 3
  @test dim(base_space(foah13_weier)) == 3
  @test dim(base_space(foah14_weier)) == 3
  @test dim(base_space(foah15_weier)) == 3
  @test dim(base_space(foah16_weier)) == 3
  @test dim(ambient_space(foah1_weier)) == 5
  @test dim(ambient_space(foah2_weier)) == 5
  @test dim(ambient_space(foah3_weier)) == 5
  @test dim(ambient_space(foah4_weier)) == 5
  @test dim(ambient_space(foah5_weier)) == 5
  @test dim(ambient_space(foah6_weier)) == 5
  @test dim(ambient_space(foah7_weier)) == 5
  @test dim(ambient_space(foah8_weier)) == 5
  @test dim(ambient_space(foah9_weier)) == 5
  @test dim(ambient_space(foah10_weier)) == 5
  @test dim(ambient_space(foah11_weier)) == 5
  @test dim(ambient_space(foah12_weier)) == 5
  @test dim(ambient_space(foah13_weier)) == 5
  @test dim(ambient_space(foah14_weier)) == 5
  @test dim(ambient_space(foah15_weier)) == 5
  @test dim(ambient_space(foah16_weier)) == 5
  @test is_base_space_fully_specified(foah1_weier) == false
  @test is_base_space_fully_specified(foah2_weier) == false
  @test is_base_space_fully_specified(foah3_weier) == false
  @test is_base_space_fully_specified(foah4_weier) == false
  @test is_base_space_fully_specified(foah5_weier) == false
  @test is_base_space_fully_specified(foah6_weier) == false
  @test is_base_space_fully_specified(foah7_weier) == false
  @test is_base_space_fully_specified(foah8_weier) == false
  @test is_base_space_fully_specified(foah9_weier) == false
  @test is_base_space_fully_specified(foah10_weier) == false
  @test is_base_space_fully_specified(foah11_weier) == false
  @test is_base_space_fully_specified(foah12_weier) == false
  @test is_base_space_fully_specified(foah13_weier) == false
  @test is_base_space_fully_specified(foah14_weier) == false
  @test is_base_space_fully_specified(foah15_weier) == false
  @test is_base_space_fully_specified(foah16_weier) == false
  @test model_description(foah1_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_1"
  @test model_description(foah2_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_2"
  @test model_description(foah3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_3"
  @test model_description(foah4_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_4"
  @test model_description(foah5_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_5"
  @test model_description(foah6_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_6"
  @test model_description(foah7_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_7"
  @test model_description(foah8_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_8"
  @test model_description(foah9_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_9"
  @test model_description(foah10_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_10"
  @test model_description(foah11_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_11"
  @test model_description(foah12_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_12"
  @test model_description(foah13_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_13"
  @test model_description(foah14_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_14"
  @test model_description(foah15_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_15"
  @test model_description(foah16_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_16"
end



########################################################################################################
# 11: Test Weierstrass counterparts of models from F-theory on all toric hypersurfaces over concrete base
########################################################################################################

B3 = projective_space(NormalToricVariety, 3)
Kbar = anticanonical_divisor(B3)
foah1_B3_weier = weierstrass_model(foah1_B3)
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

@testset "Test Weierstrass form of models in F-theory on all toric hypersurfaces, defined over concrete base" begin
  @test dim(base_space(foah1_B3_weier)) == 3
  @test dim(base_space(foah2_B3_weier)) == 3
  @test dim(base_space(foah3_B3_weier)) == 3
  @test dim(base_space(foah4_B3_weier)) == 3
  @test dim(base_space(foah5_B3_weier)) == 3
  @test dim(base_space(foah6_B3_weier)) == 3
  @test dim(base_space(foah7_B3_weier)) == 3
  @test dim(base_space(foah8_B3_weier)) == 3
  @test dim(base_space(foah9_B3_weier)) == 3
  @test dim(base_space(foah10_B3_weier)) == 3
  @test dim(base_space(foah11_B3_weier)) == 3
  @test dim(base_space(foah12_B3_weier)) == 3
  @test dim(base_space(foah13_B3_weier)) == 3
  @test dim(base_space(foah14_B3_weier)) == 3
  @test dim(base_space(foah15_B3_weier)) == 3
  @test dim(base_space(foah16_B3_weier)) == 3
  @test dim(ambient_space(foah1_B3_weier)) == 5
  @test dim(ambient_space(foah2_B3_weier)) == 5
  @test dim(ambient_space(foah3_B3_weier)) == 5
  @test dim(ambient_space(foah4_B3_weier)) == 5
  @test dim(ambient_space(foah5_B3_weier)) == 5
  @test dim(ambient_space(foah6_B3_weier)) == 5
  @test dim(ambient_space(foah7_B3_weier)) == 5
  @test dim(ambient_space(foah8_B3_weier)) == 5
  @test dim(ambient_space(foah9_B3_weier)) == 5
  @test dim(ambient_space(foah10_B3_weier)) == 5
  @test dim(ambient_space(foah11_B3_weier)) == 5
  @test dim(ambient_space(foah12_B3_weier)) == 5
  @test dim(ambient_space(foah13_B3_weier)) == 5
  @test dim(ambient_space(foah14_B3_weier)) == 5
  @test dim(ambient_space(foah15_B3_weier)) == 5
  @test dim(ambient_space(foah16_B3_weier)) == 5
  @test is_base_space_fully_specified(foah1_B3_weier) == true
  @test is_base_space_fully_specified(foah2_B3_weier) == true
  @test is_base_space_fully_specified(foah3_B3_weier) == true
  @test is_base_space_fully_specified(foah4_B3_weier) == true
  @test is_base_space_fully_specified(foah5_B3_weier) == true
  @test is_base_space_fully_specified(foah6_B3_weier) == true
  @test is_base_space_fully_specified(foah7_B3_weier) == true
  @test is_base_space_fully_specified(foah8_B3_weier) == true
  @test is_base_space_fully_specified(foah9_B3_weier) == true
  @test is_base_space_fully_specified(foah10_B3_weier) == true
  @test is_base_space_fully_specified(foah11_B3_weier) == true
  @test is_base_space_fully_specified(foah12_B3_weier) == true
  @test is_base_space_fully_specified(foah13_B3_weier) == true
  @test is_base_space_fully_specified(foah14_B3_weier) == true
  @test is_base_space_fully_specified(foah15_B3_weier) == true
  @test is_base_space_fully_specified(foah16_B3_weier) == true
  @test model_description(foah1_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_1"
  @test model_description(foah2_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_2"
  @test model_description(foah3_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_3"
  @test model_description(foah4_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_4"
  @test model_description(foah5_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_5"
  @test model_description(foah6_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_6"
  @test model_description(foah7_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_7"
  @test model_description(foah8_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_8"
  @test model_description(foah9_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_9"
  @test model_description(foah10_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_10"
  @test model_description(foah11_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_11"
  @test model_description(foah12_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_12"
  @test model_description(foah13_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_13"
  @test model_description(foah14_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_14"
  @test model_description(foah15_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_15"
  @test model_description(foah16_B3_weier) == "F-theory weierstrass model dual to hypersurface model with fiber ambient space F_16"
  @test parent(explicit_model_sections(foah1_B3_weier)["s7"]) == cox_ring(base_space(foah1_B3_weier))
  @test parent(explicit_model_sections(foah2_B3_weier)["b7"]) == cox_ring(base_space(foah2_B3_weier))
  @test parent(explicit_model_sections(foah3_B3_weier)["s7"]) == cox_ring(base_space(foah3_B3_weier))
  @test parent(explicit_model_sections(foah4_B3_weier)["d4"]) == cox_ring(base_space(foah4_B3_weier))
  @test parent(explicit_model_sections(foah5_B3_weier)["s7"]) == cox_ring(base_space(foah5_B3_weier))
  @test parent(explicit_model_sections(foah6_B3_weier)["s7"]) == cox_ring(base_space(foah6_B3_weier))
  @test parent(explicit_model_sections(foah7_B3_weier)["s7"]) == cox_ring(base_space(foah7_B3_weier))
  @test parent(explicit_model_sections(foah8_B3_weier)["s7"]) == cox_ring(base_space(foah8_B3_weier))
  @test parent(explicit_model_sections(foah9_B3_weier)["s7"]) == cox_ring(base_space(foah9_B3_weier))
  @test parent(explicit_model_sections(foah10_B3_weier)["s5"]) == cox_ring(base_space(foah10_B3_weier))
  @test parent(explicit_model_sections(foah11_B3_weier)["s5"]) == cox_ring(base_space(foah11_B3_weier))
  @test parent(explicit_model_sections(foah12_B3_weier)["s7"]) == cox_ring(base_space(foah12_B3_weier))
  @test parent(explicit_model_sections(foah13_B3_weier)["s1"]) == cox_ring(base_space(foah13_B3_weier))
  @test parent(explicit_model_sections(foah14_B3_weier)["s7"]) == cox_ring(base_space(foah14_B3_weier))
  @test parent(explicit_model_sections(foah15_B3_weier)["s7"]) == cox_ring(base_space(foah15_B3_weier))
  @test parent(explicit_model_sections(foah16_B3_weier)["s7"]) == cox_ring(base_space(foah16_B3_weier))
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
