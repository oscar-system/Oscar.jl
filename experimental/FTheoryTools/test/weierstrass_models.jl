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


#############################################################
# 2: Weierstrass models over generic base space
#############################################################

auxiliary_base_ring, (f, g, Kbar, u) = QQ["f", "g", "Kbar", "u"]
auxiliary_base_grading = [4 6 1 0; 0 0 0 1]
w2 = weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, 3, f, g)

@testset "Attributes of Weierstrass models over generic base space" begin
  @test parent(weierstrass_section_f(w2)) == coordinate_ring(base_space(w2))
  @test parent(weierstrass_section_g(w2)) == coordinate_ring(base_space(w2))
  @test parent(weierstrass_polynomial(w2)) == coordinate_ring(ambient_space(w2))
  @test dim(base_space(w2)) == 3
  @test dim(ambient_space(w2)) == 5
end

@testset "Error messages in Weierstrass models over generic base space" begin
  @test_throws ArgumentError weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, 3, f, sec_f)
  @test_throws ArgumentError weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, 0, f, g)
end
