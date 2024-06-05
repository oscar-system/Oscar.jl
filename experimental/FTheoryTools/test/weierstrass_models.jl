#############################################################
# 1: Weierstrass models over concrete base space
#############################################################

base = sample_toric_variety()
sec_f = generic_section(anticanonical_bundle(projective_space(NormalToricVariety,3))^4)
sec_g = generic_section(anticanonical_bundle(base)^6)
w = weierstrass_model(base; completeness_check = false)

@testset "Attributes of Weierstrass models over concrete base spaces" begin
  @test parent(weierstrass_section_f(w)) == cox_ring(base_space(w))
  @test parent(weierstrass_section_g(w)) == cox_ring(base_space(w))
  @test parent(weierstrass_polynomial(w)) == cox_ring(ambient_space(w))
  @test parent(discriminant(w)) == cox_ring(base_space(w))
  @test dim(base_space(w)) == 3
  @test dim(ambient_space(w)) == 5
  @test is_smooth(ambient_space(w)) == false
  @test toric_variety(calabi_yau_hypersurface(w)) == ambient_space(w)
  @test length(singular_loci(w)) == 1
  @test is_base_space_fully_specified(w) == true
  @test is_base_space_fully_specified(w) == is_base_space_fully_specified(weierstrass_model(w))
  @test is_partially_resolved(w) == false
end

@testset "Error messages in Weierstrass models over concrete base spaces" begin
  @test_throws ArgumentError weierstrass_model(base, sec_f, sec_g; completeness_check = false)
end

Kbar = anticanonical_bundle(base_space(w))
my_choice = Dict("f" => basis_of_global_sections(Kbar^4)[1])
w2 = tune(w, my_choice; completeness_check = false)

B2 = projective_space(NormalToricVariety, 2)
b = torusinvariant_prime_divisors(B2)[1]
w3 = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, defining_classes = Dict("b" => b), completeness_check = false)

@testset "Saving and loading Weierstrass models over concrete base space" begin
  mktempdir() do path
    test_save_load_roundtrip(path, w3) do loaded
      @test weierstrass_polynomial(w3) == weierstrass_polynomial(loaded)
      @test weierstrass_section_f(w3) == weierstrass_section_f(loaded)
      @test weierstrass_section_g(w3) == weierstrass_section_g(loaded)
      @test base_space(w3) == base_space(loaded)
      @test ambient_space(w3) == ambient_space(loaded)
      @test is_base_space_fully_specified(w3) == is_base_space_fully_specified(loaded)
      @test is_partially_resolved(w3) == is_partially_resolved(loaded)
    end
  end
end

x1, x2, x3, x, y, z = gens(parent(weierstrass_polynomial(w3)))
new_weierstrass_polynomial = x^3 - y^2 - x3^12 * x * z^4
tuned_w3 = tune(w3, new_weierstrass_polynomial)

@testset "Tuning of a Weierstrass model over a concrete toric base" begin
  @test base_space(w2) == base_space(w)
  @test weierstrass_section_f(w2) == my_choice["f"]
  @test weierstrass_section_f(w2) != weierstrass_section_f(w)
  @test weierstrass_section_g(w2) == weierstrass_section_g(w)
  @test hypersurface_equation(tuned_w3) == new_weierstrass_polynomial
  @test w3 == tune(w3, weierstrass_polynomial(w3))
  @test base_space(tuned_w3) == base_space(w3)
  @test fiber_ambient_space(tuned_w3) == fiber_ambient_space(w3)
end

other_Kbar = anticanonical_bundle(projective_space(NormalToricVariety, 3))

@testset "Error messages from tuning Weierstrass models over concrete toric base" begin
  @test_throws ArgumentError tune(w, Dict("f" => basis_of_global_sections(Kbar^2)[1]))
  @test_throws ArgumentError tune(w, Dict("f" => basis_of_global_sections(other_Kbar^4)[1]))
  @test_throws ArgumentError tune(w3, basis_of_global_sections(other_Kbar^4)[1])
  @test_throws ArgumentError tune(w3, x1)
end


#############################################################
# 2: Weierstrass models over generic base space
#############################################################

auxiliary_base_ring, (f, g, Kbar, u) = QQ["f", "g", "Kbar", "u"]
auxiliary_base_grading = [4 6 1 0; 0 0 0 1]
w4 = weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, 3, f, g)

@testset "Attributes of Weierstrass models over generic base space" begin
  @test parent(weierstrass_section_f(w4)) == coordinate_ring(base_space(w4))
  @test parent(weierstrass_section_g(w4)) == coordinate_ring(base_space(w4))
  @test parent(weierstrass_polynomial(w4)) == coordinate_ring(ambient_space(w4))
  @test dim(base_space(w4)) == 3
  @test dim(ambient_space(w4)) == 5
end

@testset "Error messages in Weierstrass models over generic base space" begin
  @test_throws ArgumentError weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, 3, f, sec_f)
  @test_throws ArgumentError weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, 0, f, g)
  @test_throws ArgumentError tune(w3, basis_of_global_sections(other_Kbar)[1])
  @test_throws ArgumentError tune(w4, Dict("f" => basis_of_global_sections(other_Kbar)[1]))
end
