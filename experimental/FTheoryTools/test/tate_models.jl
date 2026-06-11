#############################################################
# 1: Global Tate models over concrete base space
#############################################################

using Random
our_rng = Random.Xoshiro(1234)

# Base space and Tate sections
P3 = projective_space(NormalToricVariety, 3)
section_a1 = generic_section(
  anticanonical_bundle(projective_space(NormalToricVariety, 3)); rng=our_rng
)
section_a2 = generic_section(anticanonical_bundle(P3)^2; rng=our_rng)
section_a3 = generic_section(anticanonical_bundle(P3)^3; rng=our_rng)
section_a4 = generic_section(anticanonical_bundle(P3)^4; rng=our_rng)
section_a6 = generic_section(anticanonical_bundle(P3)^6; rng=our_rng)

# Global Tate model over P3
tate_P3 = global_tate_model(P3; completeness_check=false, rng=our_rng)

@testset "Attributes of global Tate models over concrete base space" begin
  @test parent(tate_section_a1(tate_P3)) == coordinate_ring(base_space(tate_P3))
  @test parent(tate_section_a2(tate_P3)) == coordinate_ring(base_space(tate_P3))
  @test parent(tate_section_a3(tate_P3)) == coordinate_ring(base_space(tate_P3))
  @test parent(tate_section_a4(tate_P3)) == coordinate_ring(base_space(tate_P3))
  @test parent(tate_section_a6(tate_P3)) == coordinate_ring(base_space(tate_P3))
  @test parent(tate_polynomial(tate_P3)) == coordinate_ring(ambient_space(tate_P3))
  @test parent(discriminant(tate_P3)) == coordinate_ring(base_space(tate_P3))
  @test dim(base_space(tate_P3)) == 3
  @test dim(ambient_space(tate_P3)) == 5
  @test is_base_space_fully_specified(tate_P3) == true
  @test is_base_space_fully_specified(tate_P3) ==
    is_base_space_fully_specified(weierstrass_model(tate_P3))
  @test is_smooth(ambient_space(tate_P3)) == false
  @test toric_variety(calabi_yau_hypersurface(tate_P3)) == ambient_space(tate_P3)
  @test is_partially_resolved(tate_P3) == false
end

@testset "Error messages in global Tate models over concrete base space" begin
  @test_throws ArgumentError global_tate_model(
    P3, [section_a2, section_a3, section_a4, section_a6]; completeness_check=false
  )
  @test_throws ArgumentError global_tate_model(
    P3,
    [section_a1, section_a2, section_a3, section_a4, section_a6];
    completeness_check=false,
  )
end

@test_skip @testset "Serialization Tests" begin
  P3b = projective_space(NormalToricVariety, 3)
  divisor_w = torusinvariant_prime_divisors(P3b)[1]

  literature_tate_P3 = literature_model(;
    arxiv_id="1109.3454",
    equation="3.1",
    base_space=P3b,
    defining_classes=Dict("w" => divisor_w),
    completeness_check=false,
    rng=our_rng,
  )

  @testset "Saving and loading global Tate model over concrete base space" begin
    mktempdir() do path
      test_save_load_roundtrip(path, literature_tate_P3) do loaded
        @test tate_polynomial(literature_tate_P3) == tate_polynomial(loaded)
        @test tate_section_a1(literature_tate_P3) == tate_section_a1(loaded)
        @test tate_section_a2(literature_tate_P3) == tate_section_a2(loaded)
        @test tate_section_a3(literature_tate_P3) == tate_section_a3(loaded)
        @test tate_section_a4(literature_tate_P3) == tate_section_a4(loaded)
        @test tate_section_a6(literature_tate_P3) == tate_section_a6(loaded)
        @test base_space(literature_tate_P3) == base_space(loaded)
        @test ambient_space(literature_tate_P3) == ambient_space(loaded)
        @test is_base_space_fully_specified(literature_tate_P3) ==
          is_base_space_fully_specified(loaded)
        for (key, value) in literature_tate_P3.__attrs
          if value isa String || value isa Vector{String} || value isa Bool
            @test literature_tate_P3.__attrs[key] == loaded.__attrs[key]
          end
        end
      end
    end
  end

  another_tate_P3 = global_tate_model(P3; completeness_check=false, rng=our_rng)

  @testset "Saving and loading another global Tate model over concrete base space" begin
    mktempdir() do path
      test_save_load_roundtrip(path, another_tate_P3) do loaded
        @test tate_polynomial(another_tate_P3) == tate_polynomial(loaded)
        @test tate_section_a1(another_tate_P3) == tate_section_a1(loaded)
        @test tate_section_a2(another_tate_P3) == tate_section_a2(loaded)
        @test tate_section_a3(another_tate_P3) == tate_section_a3(loaded)
        @test tate_section_a4(another_tate_P3) == tate_section_a4(loaded)
        @test tate_section_a6(another_tate_P3) == tate_section_a6(loaded)
        @test base_space(another_tate_P3) == base_space(loaded)
        @test ambient_space(another_tate_P3) == ambient_space(loaded)
        @test is_base_space_fully_specified(another_tate_P3) ==
          is_base_space_fully_specified(loaded)
        @test is_partially_resolved(another_tate_P3) == is_partially_resolved(loaded)
        @test explicit_model_sections(another_tate_P3) == explicit_model_sections(loaded)
        @test model_section_parametrization(another_tate_P3) ==
          model_section_parametrization(loaded)
        @test defining_classes(another_tate_P3) == defining_classes(loaded)
      end
    end
  end
end

#############################################################
# 2: Global Tate models over arbitrary base space
#############################################################

# rings needed for constructions
aux_ring_star0, (a1pp, a2pp, a3pp, a4pp, a6pp, vp, wp, mp) = QQ[
  :a1pp, :a2pp, :a3pp, :a4pp, :a6pp, :vp, :wp, :mp
]
aux_ring_tate, (a1p, a2p, a3p, a4p, a6p, v, w) = QQ[:a1p, :a2p, :a3p, :a4p, :a6p, :v, :w]

# construct Tate models over arbitrary base
tate_generic_I1 = global_tate_model(
  aux_ring_tate,
  [1 2 3 4 6 0 0; 0 0 -2 -2 -2 1 2],
  3,
  [a1p, a2p, a3p * (v^2 + w), a4p * (v^2 + w), a6p * (v^2 + w)],
)

tate_generic_I5s = global_tate_model(
  aux_ring_tate,
  [1 2 3 4 6 0 0; 0 -1 -2 -3 -5 1 2],
  3,
  [a1p, a2p * v, a3p * v^2, a4p * v^3, a6p * v^5],
)

tate_generic_nonminimal = global_tate_model(
  aux_ring_tate,
  [1 2 3 4 6 0 0; -2 -4 -6 -8 -12 1 2],
  3,
  [
    a1p * (v^2 + w),
    a2p * (v^2 + w)^2,
    a3p * (v^2 + w)^3,
    a4p * (v^2 + w)^4,
    a6p * (v^2 + w)^6,
  ],
)

@testset "Attributes of global Tate models over generic base space" begin
  @test parent(tate_section_a1(tate_generic_I5s)) ==
    coordinate_ring(base_space(tate_generic_I5s))
  @test parent(tate_section_a2(tate_generic_I5s)) ==
    coordinate_ring(base_space(tate_generic_I5s))
  @test parent(tate_section_a3(tate_generic_I5s)) ==
    coordinate_ring(base_space(tate_generic_I5s))
  @test parent(tate_section_a4(tate_generic_I5s)) ==
    coordinate_ring(base_space(tate_generic_I5s))
  @test parent(tate_section_a6(tate_generic_I5s)) ==
    coordinate_ring(base_space(tate_generic_I5s))
  @test parent(tate_polynomial(tate_generic_I5s)) ==
    coordinate_ring(ambient_space(tate_generic_I5s))
  @test parent(discriminant(tate_generic_I5s)) ==
    coordinate_ring(base_space(tate_generic_I5s))
  @test dim(base_space(tate_generic_I5s)) == 3
  @test dim(ambient_space(tate_generic_I5s)) == 5
  @test is_base_space_fully_specified(tate_generic_I5s) == false
  @test is_base_space_fully_specified(tate_generic_I5s) ==
    is_base_space_fully_specified(weierstrass_model(tate_generic_I5s))
end

@testset "Error messages in global Tate models over generic base space" begin
  @test_throws ArgumentError global_tate_model(
    aux_ring_tate, [1 2 3 4 6 0 0; 0 -1 -2 -3 -5 1 2], 3, [a1p]
  )
  @test_throws ArgumentError global_tate_model(
    aux_ring_tate,
    [1 2 3 4 6 0 0; 0 -1 -2 -3 -5 1 2],
    3,
    [a1p, a2p * v, a3p * v^2, a4p * v^3, section_a6],
  )
  @test_throws ArgumentError global_tate_model(
    aux_ring_tate,
    [1 2 3 4 6 0 0; 0 -1 -2 -3 -5 1 2],
    -1,
    [a1p, a2p * v, a3p * v^2, a4p * v^3, a6p * v^5],
  )
end

@testset "Singular loci of global Tate models over generic base space" begin
  @test [k[2:3] for k in singular_loci(tate_generic_I1; rng=our_rng)] ==
    [((0, 0, 1), "I_1"), ((0, 0, 1), "I_1")]
  @test [k[2:3] for k in singular_loci(tate_generic_nonminimal; rng=our_rng)] ==
    [((0, 0, 1), "I_1"), ((4, 6, 12), "Non-minimal")]
end
