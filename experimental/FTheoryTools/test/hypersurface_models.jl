#############################################################
# 1: Hypersurface models over concrete bases
#############################################################

base = projective_space(NormalToricVariety, 2)
h1 = hypersurface_model(base; completeness_check = false)

@testset "Attributes and properties of hypersurface models over concrete base space and fiber ambient space P231" begin
  @test parent(hypersurface_equation(h1)) == cox_ring(ambient_space(h1))
  @test dim(base_space(h1)) == dim(base)
  @test is_smooth(fiber_ambient_space(h1)) == false
  @test is_simplicial(fiber_ambient_space(h1)) == true
  @test [string(g) for g in gens(cox_ring(fiber_ambient_space(h1)))] == ["x", "y", "z"]
  @test toric_variety(calabi_yau_hypersurface(h1)) == ambient_space(h1)
  @test is_base_space_fully_specified(h1) == true
  @test is_partially_resolved(h1) == false
end

ambient_space_of_fiber = projective_space(NormalToricVariety, 2)
set_coordinate_names(ambient_space_of_fiber, ["x", "y", "z"])
D1 = 2 * anticanonical_divisor_class(base)
D2 = 3 * anticanonical_divisor_class(base)
h2 = hypersurface_model(base, ambient_space_of_fiber, D1, D2; completeness_check = false)

@testset "Attributes and properties of hypersurface models over concrete base space and fiber ambient space P2" begin
  @test parent(hypersurface_equation(h2)) == cox_ring(ambient_space(h2))
  @test dim(base_space(h2)) == dim(base)
  @test fiber_ambient_space(h2) == ambient_space_of_fiber
  @test is_smooth(fiber_ambient_space(h2)) == true
  @test [string(g) for g in gens(cox_ring(fiber_ambient_space(h2)))] == ["x", "y", "z"]
  @test toric_variety(calabi_yau_hypersurface(h2)) == ambient_space(h2)
  @test is_base_space_fully_specified(h2) == true
end

@testset "Error messages in hypersurface models over concrete base spaces" begin
  @test_throws ArgumentError hypersurface_model(affine_space(NormalToricVariety, 3))
  @test_throws ArgumentError weierstrass_model(h1)
  @test_throws ArgumentError global_tate_model(h1)
  @test_throws ArgumentError discriminant(h1)
  @test_throws ArgumentError singular_loci(h1)
  @test_throws ArgumentError weierstrass_model(h2)
  @test_throws ArgumentError global_tate_model(h2)
  @test_throws ArgumentError discriminant(h2)
  @test_throws ArgumentError singular_loci(h2)
end

@testset "Saving and loading hypersurface models over concrete base space" begin
  mktempdir() do path
    test_save_load_roundtrip(path, h2) do loaded
      @test hypersurface_equation(h2) == hypersurface_equation(loaded)
      @test base_space(h2) == base_space(loaded)
      @test ambient_space(h2) == ambient_space(loaded)
      @test fiber_ambient_space(h2) == fiber_ambient_space(loaded)
      @test is_base_space_fully_specified(h2) == is_base_space_fully_specified(loaded)
      @test is_partially_resolved(h2) == is_partially_resolved(loaded)
    end
  end
end

h3 = hypersurface_model_over_projective_space(2)
R = parent(hypersurface_equation(h3))
new_poly = gens(R)[4]^3
h4 = tune(h3, new_poly)

@testset "Tune hypersurface model" begin
  @test h3 == tune(h3, hypersurface_equation(h3))
  @test hypersurface_equation(h4) == new_poly
  @test hypersurface_equation(h4) != hypersurface_equation(h3)
  @test fiber_ambient_space(h4) == fiber_ambient_space(h3)
end

@testset "Errors from tuning hypersurface models" begin
  @test_throws ArgumentError tune(h3, new_poly^2)
  @test_throws ArgumentError tune(h3, hypersurface_equation(h1))
end

# Currently, none of the hypersurface models in our database has corresponding Weierstrass/Tate models.
# This code thus only tests if the code works, but the assignment is mathematically speaking wrong.
w_model = weierstrass_model_over_projective_space(2)
gt_model = global_tate_model_over_projective_space(2)
set_weierstrass_model(h3, w_model)
set_global_tate_model(h3, gt_model)

@testset "Test assignment of Weierstrass and global Tate models to hypersurface models" begin
  @test weierstrass_model(h3) == w_model
  @test global_tate_model(h3) == gt_model
end

B2 = projective_space(NormalToricVariety, 2)
b = torusinvariant_prime_divisors(B2)[1]
h5 = literature_model(arxiv_id = "1208.2695", equation = "B.5", base_space = B2, model_sections = Dict("b" => b))

@testset "Attributes and properties of hypersurface literature models over concrete base space" begin
  @test parent(hypersurface_equation(h5)) == cox_ring(ambient_space(h5))
  @test dim(base_space(h5)) == 2
  @test is_smooth(fiber_ambient_space(h5)) == false
  @test is_simplicial(fiber_ambient_space(h5)) == true
  @test [string(g) for g in gens(cox_ring(fiber_ambient_space(h5)))] == ["u", "w", "v"]
  @test is_base_space_fully_specified(h5) == true
  @test is_partially_resolved(h5) == false
end



#############################################################
# 2: Hypersurface models over generic base space
#############################################################

auxiliary_base_vars = ["a1", "a21", "a32", "a43", "a65", "w"]
auxiliary_base_grading = [1 2 3 4 6 0; 0 -1 -2 -3 -5 1]
D1 = [4,0]
D2 = [6,0]
d = 3
ambient_space_of_fiber_2 = weighted_projective_space(NormalToricVariety, [2,3,1])
set_coordinate_names(ambient_space_of_fiber_2, ["x", "y", "z"])
auxiliary_ambient_ring, (a1, a21, a32, a43, a65, w, x, y, z)  = QQ["a1", "a21", "a32", "a43", "a65", "w", "x", "y", "z"]
p = x^3 - y^2 - x * y * z * a1 + x^2 * z^2 * a21 * w - y * z^3 * a32 * w^2 + x * z^4 * a43 * w^3 + z^6 * a65 * w^5
h6 = hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, ambient_space_of_fiber_2, D1, D2, p)

@testset "Attributes and properties of hypersurface models over concrete base space and fiber ambient space P2" begin
  @test parent(hypersurface_equation(h6)) == coordinate_ring(ambient_space(h6))
  @test dim(base_space(h6)) == d
  @test fiber_ambient_space(h6) == ambient_space_of_fiber_2
  @test is_smooth(fiber_ambient_space(h6)) == false
  @test is_simplicial(fiber_ambient_space(h6)) == true
  @test [string(g) for g in gens(cox_ring(fiber_ambient_space(h6)))] == ["x", "y", "z"]
  @test is_base_space_fully_specified(h6) == false
  @test is_partially_resolved(h6) == false
end

@testset "Error messages in hypersurface models over not fully specified base spaces" begin
  @test_throws ArgumentError hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, -1, ambient_space_of_fiber_2, D1, D2, p)
  @test_throws ArgumentError tune(h6, hypersurface_equation(h6))
end

h7 = literature_model(arxiv_id = "1208.2695", equation = "B.5")

@testset "Attributes and properties of hypersurface models over concrete base space and fiber ambient space P2" begin
  @test parent(hypersurface_equation(h7)) == coordinate_ring(ambient_space(h7))
  @test dim(base_space(h7)) == 2
  @test is_smooth(fiber_ambient_space(h7)) == false
  @test is_simplicial(fiber_ambient_space(h7)) == true
  @test [string(g) for g in gens(cox_ring(fiber_ambient_space(h7)))] == ["u", "w", "v"]
  @test is_base_space_fully_specified(h7) == false
  @test is_partially_resolved(h7) == false
end
