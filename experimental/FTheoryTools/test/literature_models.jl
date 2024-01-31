#############################################################
# 1: Literature Tate model over concrete base
#############################################################

B3 = projective_space(NormalToricVariety, 3)
w = torusinvariant_prime_divisors(B3)[1]
t1 = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, model_sections = Dict("w" => w), completeness_check = false)

@testset "Test defining data for literature Tate model over concrete base" begin
  @test parent(tate_section_a1(t1)) == cox_ring(base_space(t1))
  @test parent(tate_section_a2(t1)) == cox_ring(base_space(t1))
  @test parent(tate_section_a3(t1)) == cox_ring(base_space(t1))
  @test parent(tate_section_a4(t1)) == cox_ring(base_space(t1))
  @test parent(tate_section_a6(t1)) == cox_ring(base_space(t1))
  @test parent(tate_polynomial(t1)) == cox_ring(ambient_space(t1))
  @test parent(discriminant(t1)) == cox_ring(base_space(t1))
  @test length(singular_loci(t1)) == 2
  @test dim(base_space(t1)) == 3
  @test dim(ambient_space(t1)) == 5
  @test is_base_space_fully_specified(t1) == true
  @test is_base_space_fully_specified(t1) == is_base_space_fully_specified(weierstrass_model(t1))
  @test is_smooth(ambient_space(t1)) == false
  @test toric_variety(calabi_yau_hypersurface(t1)) == ambient_space(t1)
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
  @test journal_pages(t1) == "1–47"
  @test journal_volume(t1) == "858"
  @test journal_year(t1) ==  "2012"
  @test literature_identifier(t1) ==  "1109_3454"
  @test model_description(t1) == "SU(5)xU(1) restricted Tate model"
  @test paper_authors(t1) == ["Sven Krause", "Christoph Mayrhofer", "Timo Weigand"]
  @test paper_buzzwords(t1) == ["GUT model", "Tate", "U(1)", "SU(5)"]
  @test paper_description(t1) == "SU(5)xU(1) restricted Tate model"
  @test paper_title(t1) == "\$G_4\$ flux, chiral matter and singularity resolution in F-theory compactifications"
  @test resolutions(t1) == [[[["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]], ["e1", "e4", "e2", "e3", "s"]]]
  @test length(resolution_generating_sections(t1)) == 1
  @test length(resolution_zero_sections(t1)) == 1
  @test weighted_resolutions(t1) == [[[[["x", "y", "w"], [1, 1, 1]], [["x", "y", "w"], [1, 2, 1]], [["x", "y", "w"], [2, 2, 1]], [["x", "y", "w"], [2, 3, 1]], [["x", "y"], [1, 1]]], ["e1", "e4", "e2", "e3", "s"]]]
  @test length(weighted_resolution_generating_sections(t1)) == 1
  @test length(weighted_resolution_zero_sections(t1)) == 1
end

@testset "Test error messages for literature Tate model over concrete base" begin
  @test_throws ArgumentError associated_literature_models(t1)
  @test_throws ArgumentError journal_report_numbers(t1)
  @test_throws ArgumentError model_parameters(t1)
  @test_throws ArgumentError related_literature_models(t1)
end

set_model_description(t1, "Testing...")

@testset "Test modifying the model description for literature Tate model over concrete base" begin
  @test model_description(t1) == "Testing..."
end

t2 = resolve(t1, 1)

@testset "Test resolving literature Tate model over concrete base" begin
  @test is_smooth(ambient_space(t2)) == false
  @test is_partially_resolved(t2) == true
  @test base_space(t1) == base_space(t2)
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
w1 = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, model_sections = Dict("b" => b), completeness_check = false)

@testset "Test defining data for literature Weierstrass model over concrete base" begin
  @test parent(weierstrass_section_f(w1)) == cox_ring(base_space(w1))
  @test parent(weierstrass_section_g(w1)) == cox_ring(base_space(w1))
  @test parent(weierstrass_polynomial(w1)) == cox_ring(ambient_space(w1))
  @test parent(discriminant(w1)) == cox_ring(base_space(w1))
  @test length(singular_loci(w1)) == 1
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
  @test associated_literature_models(w1) == ["1208_2695-1"]
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
  @test_throws ArgumentError related_literature_models(w1)
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
  @test length(singular_loci(t3)) == 2
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
  @test resolutions(t3) == [[[["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]], ["e1", "e4", "e2", "e3", "s"]]]
  @test length(resolution_generating_sections(t3)) == 1
  @test length(resolution_zero_sections(t3)) == 1
  @test weighted_resolutions(t3) == [[[[["x", "y", "w"], [1, 1, 1]], [["x", "y", "w"], [1, 2, 1]], [["x", "y", "w"], [2, 2, 1]], [["x", "y", "w"], [2, 3, 1]], [["x", "y"], [1, 1]]], ["e1", "e4", "e2", "e3", "s"]]]
  @test length(weighted_resolution_generating_sections(t3)) == 1
  @test length(weighted_resolution_zero_sections(t3)) == 1
end

@testset "Test error messages for literature Tate model over arbitrary base" begin
  @test_throws ArgumentError associated_literature_models(t3)
  @test_throws ArgumentError journal_report_numbers(t3)
  @test_throws ArgumentError model_parameters(t3)
  @test_throws ArgumentError related_literature_models(t3)
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
t4b = literature_model(arxiv_id = "1212.2949", equation = "3.2", model_parameters = Dict("k" => 2), base_space = B3, model_sections = Dict("ζ0" => ζ0))
t5b = literature_model(arxiv_id = "1212.2949", equation = "3.42", model_parameters = Dict("k" => 2), base_space = B3, model_sections = Dict("ζ0" => ζ0))
t6b = literature_model(arxiv_id = "1212.2949", equation = "4.1", model_parameters = Dict("k" => 2), base_space = B3, model_sections = Dict("ζ0" => ζ0))
t7b = literature_model(arxiv_id = "1212.2949", equation = "4.23", model_parameters = Dict("k" => 2), base_space = B3, model_sections = Dict("ζ0" => ζ0))
t8b  = literature_model(arxiv_id = "1212.2949", equation = "5.1", base_space = B3, model_sections = Dict("ζ0" => ζ0))
t9b = literature_model(arxiv_id = "1212.2949", equation = "5.7", base_space = B3, model_sections = Dict("ζ0" => ζ0))
t10b = literature_model(arxiv_id = "1212.2949", equation = "5.13", base_space = B3, model_sections = Dict("ζ0" => ζ0))

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
  @test length(singular_loci(w2)) == 1
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
  @test associated_literature_models(w2) == ["1208_2695-1"]
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
  @test_throws ArgumentError related_literature_models(w2)
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
w4 = literature_model(arxiv_id = "1208.2695", equation = "B.19", completeness_check = false, base_space = B2, model_sections = Dict("b" =>b))
b2 = anticanonical_divisor(B2)
w5 = literature_model(arxiv_id = "1507.05954", equation = "A.1", completeness_check = false, base_space = B2, model_sections = Dict("s8" => b2, "a1" => b, "a2" => b, "a3" => b))

@testset "Test defining data for literature Weierstrass model over concrete base" begin
  @test length(singular_loci(w4)) == 1
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
w6 = literature_model(3, base_space = B2, model_sections = Dict("b" => b), completeness_check = false)

@testset "Test defining data for literature model defined by model index" begin
  @test length(singular_loci(w6)) == 1
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

@testset "Test for literature hypersurface model over arbitary base" begin
  @test parent(hypersurface_equation(h)) == coordinate_ring(ambient_space(h))
  @test dim(base_space(h)) == 2
  @test is_smooth(fiber_ambient_space(h)) == true
  @test [string(g) for g in gens(cox_ring(fiber_ambient_space(h)))] == ["u", "v", "w"]
  @test is_base_space_fully_specified(h) == false
  @test is_partially_resolved(h) == false
  @test string.(zero_section(h)) == ["0", "-b1", "a1"]
end
