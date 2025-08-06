#############################################################
# 1: Weierstrass models over concrete base space
#############################################################

using Random
our_rng = Random.Xoshiro(1234)

my_base = projective_space(NormalToricVariety, 3)
sec_f = generic_section(anticanonical_bundle(projective_space(NormalToricVariety,3))^4)
sec_g = generic_section(anticanonical_bundle(my_base)^6)
w = weierstrass_model(my_base; completeness_check = false)

@testset "Attributes of Weierstrass models over concrete base spaces" begin
  @test parent(weierstrass_section_f(w)) == cox_ring(base_space(w))
  @test parent(weierstrass_section_g(w)) == cox_ring(base_space(w))
  @test parent(weierstrass_polynomial(w)) == cox_ring(ambient_space(w))
  @test parent(discriminant(w)) == cox_ring(base_space(w))
  @test dim(base_space(w)) == 3
  @test dim(ambient_space(w)) == 5
  @test is_smooth(ambient_space(w)) == false
  @test toric_variety(calabi_yau_hypersurface(w)) == ambient_space(w)
  @test length(singular_loci(w; rng = our_rng)) == 1
  @test is_base_space_fully_specified(w) == true
  @test is_partially_resolved(w) == false
end

@testset "Error messages in Weierstrass models over concrete base spaces" begin
  @test_throws ArgumentError weierstrass_model(my_base, sec_f, sec_g; completeness_check = false)
end

B2 = projective_space(NormalToricVariety, 2)
b = torusinvariant_prime_divisors(B2)[1]
w3 = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, defining_classes = Dict("b" => b), completeness_check = false)

@testset "Saving and loading Weierstrass model over concrete base space" begin
  mktempdir() do path
    test_save_load_roundtrip(path, w3) do loaded
      @test weierstrass_polynomial(w3) == weierstrass_polynomial(loaded)
      @test weierstrass_section_f(w3) == weierstrass_section_f(loaded)
      @test weierstrass_section_g(w3) == weierstrass_section_g(loaded)
      @test base_space(w3) == base_space(loaded)
      @test ambient_space(w3) == ambient_space(loaded)
      @test is_base_space_fully_specified(w3) == is_base_space_fully_specified(loaded)
      @test explicit_model_sections(w3) == explicit_model_sections(loaded)
      @test model_section_parametrization(w3) == model_section_parametrization(loaded)
      @test defining_classes(w3) == defining_classes(loaded)
      for (key, value) in w3.__attrs
        if value isa String || value isa Vector{String} || value isa Bool
          @test w3.__attrs[key] == loaded.__attrs[key]
        end
      end
    end
  end
end

w3_copy = weierstrass_model(projective_space(NormalToricVariety, 3); completeness_check = false)

@testset "Saving and loading another Weierstrass model over concrete base space" begin
  mktempdir() do path
    test_save_load_roundtrip(path, w3_copy) do loaded
      @test weierstrass_polynomial(w3_copy) == weierstrass_polynomial(loaded)
      @test weierstrass_section_f(w3_copy) == weierstrass_section_f(loaded)
      @test weierstrass_section_g(w3_copy) == weierstrass_section_g(loaded)
      @test base_space(w3_copy) == base_space(loaded)
      @test ambient_space(w3_copy) == ambient_space(loaded)
      @test is_base_space_fully_specified(w3_copy) == is_base_space_fully_specified(loaded)
      @test is_partially_resolved(w3_copy) == is_partially_resolved(loaded)
      @test explicit_model_sections(w3_copy) == explicit_model_sections(loaded)
      @test model_section_parametrization(w3_copy) == model_section_parametrization(loaded)
      @test defining_classes(w3_copy) == defining_classes(loaded)
    end
  end
end


#############################################################
# 2: Weierstrass models over generic base space
#############################################################

auxiliary_base_ring, (f, g, Kbar, u) = QQ[:f, :g, :Kbar, :u]
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
end
