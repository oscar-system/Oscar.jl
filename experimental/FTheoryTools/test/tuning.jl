#############################################################
# 1: Tuning Weierstrass models over concrete base space
#############################################################

using Random
our_rng = Random.Xoshiro(1234)

my_base = projective_space(NormalToricVariety, 3)
sec_f = generic_section(
  anticanonical_bundle(projective_space(NormalToricVariety, 3))^4; rng=our_rng
)
sec_g = generic_section(anticanonical_bundle(my_base)^6; rng=our_rng)
w = weierstrass_model(my_base; completeness_check=false, rng=our_rng)
Kbar = anticanonical_bundle(base_space(w))
my_choice = Dict("f" => basis_of_global_sections(Kbar^4)[1])
w2 = tune(w, my_choice; completeness_check=false)

# The tests associated to the below code did not actually test the tune functionality of Weierstrass models,
# but instead the tune functionality of abstract F-theory models, inherited by Weierstrass models
# This functionality has been removed for the time being, because it did not correspond to a proper tuning
# These tests have been retained for the (potential) future date when we reintroduce this functionality
#
# x1, x2, x3, x, y, z = gens(parent(weierstrass_polynomial(w3)))
# new_weierstrass_polynomial = x^3 - y^2 - x3^12 * x * z^4
# tuned_w3 = tune(w3, new_weierstrass_polynomial)

@testset "Tuning of a Weierstrass model over a concrete toric base" begin
  @test base_space(w2) == base_space(w)
  @test weierstrass_section_f(w2) == my_choice["f"]
  @test weierstrass_section_f(w2) != weierstrass_section_f(w)
  @test weierstrass_section_g(w2) == weierstrass_section_g(w)
  # @test hypersurface_equation(tuned_w3) == new_weierstrass_polynomial
  # @test w3 == tune(w3, weierstrass_polynomial(w3))
  # @test base_space(tuned_w3) == base_space(w3)
  # @test fiber_ambient_space(tuned_w3) == fiber_ambient_space(w3)
end

other_Kbar = anticanonical_bundle(projective_space(NormalToricVariety, 3))

@testset "Error messages from tuning Weierstrass models over concrete toric base" begin
  @test_throws ArgumentError tune(w, Dict("f" => basis_of_global_sections(Kbar^2)[1]))
  @test_throws ArgumentError tune(w, Dict("f" => basis_of_global_sections(other_Kbar^4)[1]))
  # Removed, see above
  # @test_throws ArgumentError tune(w3, basis_of_global_sections(other_Kbar^4)[1])
  # @test_throws ArgumentError tune(w3, x1)
end

#############################################################
# 2: Tuning Weierstrass models over generic base spaces
#############################################################

auxiliary_base_ring, (f, g, Kbar, u) = QQ[:f, :g, :Kbar, :u]
auxiliary_base_grading = [4 6 1 0; 0 0 0 1]
w4 = weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, 3, f, g)

@testset "Error messages in Weierstrass models over generic base space" begin
  # @test_throws ArgumentError tune(w3, basis_of_global_sections(other_Kbar)[1]) # Removed, see above
  @test_throws ArgumentError tune(w4, Dict("f" => basis_of_global_sections(other_Kbar)[1]))
end

#############################################################
# 3: Tuning global Tate models over concrete base space
#############################################################

my_base = projective_space(NormalToricVariety, 3)
sec_a1 = generic_section(
  anticanonical_bundle(projective_space(NormalToricVariety, 3)); rng=our_rng
)
sec_a2 = generic_section(anticanonical_bundle(my_base)^2; rng=our_rng)
sec_a3 = generic_section(anticanonical_bundle(my_base)^3; rng=our_rng)
sec_a4 = generic_section(anticanonical_bundle(my_base)^4; rng=our_rng)
sec_a6 = generic_section(anticanonical_bundle(my_base)^6; rng=our_rng)
t = global_tate_model(my_base; completeness_check=false, rng=our_rng)

Kbar = anticanonical_bundle(base_space(t))
my_choice = Dict("a1" => basis_of_global_sections(Kbar)[1])
my_choice["a2"] = basis_of_global_sections(Kbar^2)[1]
my_choice["a3"] = basis_of_global_sections(Kbar^3)[1]
my_choice["a4"] = basis_of_global_sections(Kbar^4)[1]
my_choice["a6"] = basis_of_global_sections(Kbar^6)[1]

t3 = tune(t, my_choice; completeness_check=false)

# The tests associated to the below code did not actually test the tune functionality of Tate models,
# but instead the tune functionality of abstract F-theory models, inherited by Tate models
# This functionality has been removed for the time being, because it did not correspond to a proper tuning
# These tests have been retained for the (potential) future date when we reintroduce this functionality
#
# x1, x2, x3, x4, x, y, z = gens(parent(tate_polynomial(t2)))
# new_tate_polynomial = x^3 - y^2 - x * y * z * x4^4
# tuned_t2 = tune(t2, new_tate_polynomial)

@testset "Tuning of a Tate model over a concrete toric space" begin
  @test base_space(t3) == base_space(t)
  @test tate_section_a1(t3) == my_choice["a1"]
  @test tate_section_a1(t3) != tate_section_a1(t)
  @test tate_section_a2(t3) == my_choice["a2"]
  @test tate_section_a2(t3) != tate_section_a2(t)
  @test tate_section_a3(t3) == my_choice["a3"]
  @test tate_section_a3(t3) != tate_section_a3(t)
  @test tate_section_a4(t3) == my_choice["a4"]
  @test tate_section_a4(t3) != tate_section_a4(t)
  @test tate_section_a6(t3) == my_choice["a6"]
  @test tate_section_a6(t3) != tate_section_a6(t)
  # Removed, see above
  # @test t2 == tune(t2, tate_polynomial(t2))
  # @test hypersurface_equation(tuned_t2) == new_tate_polynomial
  # @test base_space(tuned_t2) == base_space(t2)
  # @test fiber_ambient_space(tuned_t2) == fiber_ambient_space(t2)
end

other_Kbar = anticanonical_bundle(projective_space(NormalToricVariety, 3))

@testset "Error messages from tuning Tate models" begin
  @test_throws ArgumentError tune(t, Dict("a1" => basis_of_global_sections(Kbar^2)[1]))
  @test_throws ArgumentError tune(t, Dict("a1" => basis_of_global_sections(other_Kbar)[1]))
  # Removed, see above
  # @test_throws ArgumentError tune(t2, basis_of_global_sections(other_Kbar)[1])
  # @test_throws ArgumentError tune(t2, x1)
end

#############################################################
# 4: Tuning global Tate models over arbitrary base space
#############################################################

# rings needed for constructions
tate_auxiliary_base_ring, (a1p, a2p, a3p, a4p, a6p, v, w) = QQ[
  :a1p, :a2p, :a3p, :a4p, :a6p, :v, :w
];
t_i5_s = global_tate_model(
  tate_auxiliary_base_ring,
  [1 2 3 4 6 0 0; 0 -1 -2 -3 -5 1 2],
  3,
  [a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5],
);

@testset "Error messages in global Tate models over generic base space" begin
  @test_throws ArgumentError tune(
    t_i5_s, Dict("a2" => basis_of_global_sections(other_Kbar)[1])
  )
  # @test_throws ArgumentError tune(t_i5_s, basis_of_global_sections(other_Kbar)[1]) # Removed, see above
end

#############################################################
# 5: Tuning hypersurface models over concrete bases
#############################################################

B3 = projective_space(NormalToricVariety, 3)
ambient_space_of_fiber = weighted_projective_space(NormalToricVariety, [2, 3, 1])
set_coordinate_names(ambient_space_of_fiber, ["x", "y", "z"])
D1 = 2 * anticanonical_divisor_class(B3)
D2 = 3 * anticanonical_divisor_class(B3)
D3 = trivial_divisor_class(B3)
p = "x^3 - 2*y^2 + x1^16*x*z^4 + x2^24*z^6 + 13*x3^4*x*y*z"
h = hypersurface_model(
  B3, ambient_space_of_fiber, [D1, D2, D3], p; completeness_check=false
)

# The tests below did not actually test the tune functionality of hypersurface models,
# but instead the tune functionality of abstract F-theory models, inherited by hypersurface models
# This functionality has been removed for the time being, because it did not correspond to a proper tuning
# These tests have been retained for the (potential) future date when we reintroduce this functionality
#
# (x1, x2, x3, x4, x, y, z) = gens(parent(hypersurface_equation(h)))
# new_poly = x^3 - y^2 + 3*x1^16*x*z^4 - 7*x2^24*z^6
# h2 = tune(h, new_poly)

# @testset "Tune hypersurface model" begin
#   @test h == tune(h, hypersurface_equation(h))
#   @test hypersurface_equation(h2) == new_poly
#   @test hypersurface_equation(h2) != hypersurface_equation(h)
#   @test fiber_ambient_space(h2) == fiber_ambient_space(h)
# end

# These tests have also been removed, see above
# @testset "Errors from tuning hypersurface models" begin
#   @test_throws ArgumentError tune(h, new_poly^2)
#   @test_throws ArgumentError tune(h, hypersurface_equation(h3))
# end

#############################################################
# 6: Tuning hypersurface models over generic base space
#############################################################

auxiliary_base_vars = ["a1", "a21", "a32", "a43", "a65", "w"]
auxiliary_base_grading = [1 2 3 4 6 0; 0 -1 -2 -3 -5 1]
D1 = [4, 0]
D2 = [6, 0]
D3 = [0, 0]
d = 3
ambient_space_of_fiber_2 = weighted_projective_space(NormalToricVariety, [2, 3, 1])
set_coordinate_names(ambient_space_of_fiber_2, ["x", "y", "z"])
auxiliary_ambient_ring, (a1, a21, a32, a43, a65, w, x, y, z) = QQ[
  :a1, :a21, :a32, :a43, :a65, :w, :x, :y, :z
]
p =
  x^3 - y^2 - x * y * z * a1 + x^2 * z^2 * a21 * w - y * z^3 * a32 * w^2 +
  x * z^4 * a43 * w^3 + z^6 * a65 * w^5
h4 = hypersurface_model(
  auxiliary_base_vars, auxiliary_base_grading, d, ambient_space_of_fiber_2, [D1, D2, D3], p
)

@testset "Error messages in hypersurface models over not fully specified base spaces" begin
  # @test_throws ArgumentError tune(h4, hypersurface_equation(h4)) # Removed, see above
end
