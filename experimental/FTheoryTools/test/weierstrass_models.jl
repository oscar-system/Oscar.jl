#############################################################
# 1: Weierstrass models over concrete base space
#############################################################

using Random
our_rng = Random.Xoshiro(1234)

# Base space and sections
P3 = projective_space(NormalToricVariety, 3)
section_f = generic_section(
  anticanonical_bundle(projective_space(NormalToricVariety, 3))^4; rng=our_rng
)
section_g = generic_section(anticanonical_bundle(P3)^6; rng=our_rng)

# Weierstrass model over P3
weierstrass_P3 = weierstrass_model(P3; completeness_check=false, rng=our_rng)

@testset "Attributes of Weierstrass models over concrete base spaces" begin
  @test parent(weierstrass_section_f(weierstrass_P3)) ==
    coordinate_ring(base_space(weierstrass_P3))
  @test parent(weierstrass_section_g(weierstrass_P3)) ==
    coordinate_ring(base_space(weierstrass_P3))
  @test parent(weierstrass_polynomial(weierstrass_P3)) ==
    coordinate_ring(ambient_space(weierstrass_P3))
  @test parent(discriminant(weierstrass_P3)) == coordinate_ring(base_space(weierstrass_P3))
  @test dim(base_space(weierstrass_P3)) == 3
  @test dim(ambient_space(weierstrass_P3)) == 5
  @test is_smooth(ambient_space(weierstrass_P3)) == false
  @test toric_variety(calabi_yau_hypersurface(weierstrass_P3)) ==
    ambient_space(weierstrass_P3)
  @test length(singular_loci(weierstrass_P3; rng=our_rng)) == 1
  @test is_base_space_fully_specified(weierstrass_P3) == true
  @test is_partially_resolved(weierstrass_P3) == false
end

@testset "Error messages in Weierstrass models over concrete base spaces" begin
  @test_throws ArgumentError weierstrass_model(
    P3, section_f, section_g; completeness_check=false
  )
end

# Literature example over P2
P2 = projective_space(NormalToricVariety, 2)
divisor_b = torusinvariant_prime_divisors(P2)[1]

literature_model_P2 = literature_model(;
  arxiv_id="1208.2695",
  equation="B.19",
  base_space=P2,
  defining_classes=Dict("b" => divisor_b),
  completeness_check=false,
  rng=our_rng,
)

@testset "Saving and loading Weierstrass model over concrete base space" begin
  mktempdir() do path
    test_save_load_roundtrip(path, literature_model_P2) do loaded
      @test weierstrass_polynomial(literature_model_P2) == weierstrass_polynomial(loaded)
      @test weierstrass_section_f(literature_model_P2) == weierstrass_section_f(loaded)
      @test weierstrass_section_g(literature_model_P2) == weierstrass_section_g(loaded)
      @test base_space(literature_model_P2) == base_space(loaded)
      @test ambient_space(literature_model_P2) == ambient_space(loaded)
      @test is_base_space_fully_specified(literature_model_P2) ==
        is_base_space_fully_specified(loaded)
      @test explicit_model_sections(literature_model_P2) == explicit_model_sections(loaded)
      @test model_section_parametrization(literature_model_P2) ==
        model_section_parametrization(loaded)
      @test defining_classes(literature_model_P2) == defining_classes(loaded)
      for (key, value) in literature_model_P2.__attrs
        if value isa String || value isa Vector{String} || value isa Bool
          @test literature_model_P2.__attrs[key] == loaded.__attrs[key]
        end
      end
    end
  end
end

# Another example over P3
another_weierstrass_P3 = weierstrass_model(P3; completeness_check=false, rng=our_rng)

@testset "Saving and loading another Weierstrass model over concrete base space" begin
  mktempdir() do path
    test_save_load_roundtrip(path, another_weierstrass_P3) do loaded
      @test weierstrass_polynomial(another_weierstrass_P3) == weierstrass_polynomial(loaded)
      @test weierstrass_section_f(another_weierstrass_P3) == weierstrass_section_f(loaded)
      @test weierstrass_section_g(another_weierstrass_P3) == weierstrass_section_g(loaded)
      @test base_space(another_weierstrass_P3) == base_space(loaded)
      @test ambient_space(another_weierstrass_P3) == ambient_space(loaded)
      @test is_base_space_fully_specified(another_weierstrass_P3) ==
        is_base_space_fully_specified(loaded)
      @test is_partially_resolved(another_weierstrass_P3) == is_partially_resolved(loaded)
      @test explicit_model_sections(another_weierstrass_P3) ==
        explicit_model_sections(loaded)
      @test model_section_parametrization(another_weierstrass_P3) ==
        model_section_parametrization(loaded)
      @test defining_classes(another_weierstrass_P3) == defining_classes(loaded)
    end
  end
end

#############################################################
# 2: Weierstrass models over generic base space
#############################################################

auxiliary_base_ring, (f, g, Kbar, u) = QQ[:f, :g, :Kbar, :u]
auxiliary_base_grading = [4 6 1 0; 0 0 0 1]

weierstrass_generic = weierstrass_model(
  auxiliary_base_ring, auxiliary_base_grading, 3, f, g
)

@testset "Attributes of Weierstrass models over generic base space" begin
  @test parent(weierstrass_section_f(weierstrass_generic)) ==
    coordinate_ring(base_space(weierstrass_generic))
  @test parent(weierstrass_section_g(weierstrass_generic)) ==
    coordinate_ring(base_space(weierstrass_generic))
  @test parent(weierstrass_polynomial(weierstrass_generic)) ==
    coordinate_ring(ambient_space(weierstrass_generic))
  @test dim(base_space(weierstrass_generic)) == 3
  @test dim(ambient_space(weierstrass_generic)) == 5
end

@testset "Error messages in Weierstrass models over generic base space" begin
  @test_throws ArgumentError weierstrass_model(
    auxiliary_base_ring, auxiliary_base_grading, 3, f, section_f
  )
  @test_throws ArgumentError weierstrass_model(
    auxiliary_base_ring, auxiliary_base_grading, 0, f, g
  )
end
