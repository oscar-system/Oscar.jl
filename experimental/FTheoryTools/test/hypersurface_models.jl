#############################################################
# 1: Hypersurface models over concrete bases
#############################################################

B3 = projective_space(NormalToricVariety, 3)
ambient_space_of_fiber = weighted_projective_space(NormalToricVariety, [2,3,1])
set_coordinate_names(ambient_space_of_fiber, ["x", "y", "z"])
D1 = 2 * anticanonical_divisor_class(B3)
D2 = 3 * anticanonical_divisor_class(B3)
D3 = trivial_divisor_class(B3)
p = "x^3 - 2*y^2 + x1^16*x*z^4 + x2^24*z^6 + 13*x3^4*x*y*z"
h = hypersurface_model(B3, ambient_space_of_fiber, [D1, D2, D3], p; completeness_check = false)

@testset "Attributes and properties of hypersurface models over concrete base space and fiber ambient space P2" begin
  @test parent(hypersurface_equation(h)) == cox_ring(ambient_space(h))
  @test dim(base_space(h)) == dim(B3)
  @test fiber_ambient_space(h) == ambient_space_of_fiber
  @test is_smooth(fiber_ambient_space(h)) == false
  @test [string(g) for g in gens(cox_ring(fiber_ambient_space(h)))] == ["x", "y", "z"]
  @test toric_variety(calabi_yau_hypersurface(h)) == ambient_space(h)
  @test is_base_space_fully_specified(h) == true
end

@testset "Error messages in hypersurface models over concrete base spaces" begin
  @test_throws ArgumentError weierstrass_model(h)
  @test_throws ArgumentError global_tate_model(h)
  @test_throws ArgumentError discriminant(h)
  @test_throws ArgumentError singular_loci(h)
end

@testset "Saving and loading hypersurface models over concrete base space" begin
  mktempdir() do path
    test_save_load_roundtrip(path, h) do loaded
      @test hypersurface_equation(h) == hypersurface_equation(loaded)
      @test base_space(h) == base_space(loaded)
      @test ambient_space(h) == ambient_space(loaded)
      @test fiber_ambient_space(h) == fiber_ambient_space(loaded)
      @test is_base_space_fully_specified(h) == is_base_space_fully_specified(loaded)
      @test is_partially_resolved(h) == is_partially_resolved(loaded)
    end
  end
end

(x1, x2, x3, x4, x, y, z) = gens(parent(hypersurface_equation(h)))
new_poly = x^3 - y^2 + 3*x1^16*x*z^4 - 7*x2^24*z^6
h2 = tune(h, new_poly)

@testset "Tune hypersurface model" begin
  @test h == tune(h, hypersurface_equation(h))
  @test hypersurface_equation(h2) == new_poly
  @test hypersurface_equation(h2) != hypersurface_equation(h)
  @test fiber_ambient_space(h2) == fiber_ambient_space(h)
end

B2 = projective_space(NormalToricVariety, 2)
b = torusinvariant_prime_divisors(B2)[1]
h3 = literature_model(arxiv_id = "1208.2695", equation = "B.5", base_space = B2, model_sections = Dict("b" => b))

@testset "Errors from tuning hypersurface models" begin
  @test_throws ArgumentError tune(h, new_poly^2)
  @test_throws ArgumentError tune(h, hypersurface_equation(h3))
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

@testset "Attributes and properties of hypersurface literature models over concrete base space" begin
  @test parent(hypersurface_equation(h3)) == cox_ring(ambient_space(h3))
  @test dim(base_space(h3)) == 2
  @test is_smooth(fiber_ambient_space(h3)) == false
  @test is_simplicial(fiber_ambient_space(h3)) == true
  @test [string(g) for g in gens(cox_ring(fiber_ambient_space(h3)))] == ["u", "w", "v"]
  @test is_base_space_fully_specified(h3) == true
  @test is_partially_resolved(h3) == false
end



#############################################################
# 2: Hypersurface models over generic base space
#############################################################

auxiliary_base_vars = ["a1", "a21", "a32", "a43", "a65", "w"]
auxiliary_base_grading = [1 2 3 4 6 0; 0 -1 -2 -3 -5 1]
D1 = [4,0]
D2 = [6,0]
D3 = [0,0]
d = 3
ambient_space_of_fiber_2 = weighted_projective_space(NormalToricVariety, [2,3,1])
set_coordinate_names(ambient_space_of_fiber_2, ["x", "y", "z"])
auxiliary_ambient_ring, (a1, a21, a32, a43, a65, w, x, y, z)  = QQ["a1", "a21", "a32", "a43", "a65", "w", "x", "y", "z"]
p = x^3 - y^2 - x * y * z * a1 + x^2 * z^2 * a21 * w - y * z^3 * a32 * w^2 + x * z^4 * a43 * w^3 + z^6 * a65 * w^5
h4 = hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, ambient_space_of_fiber_2, [D1, D2, D3], p)

@testset "Attributes and properties of hypersurface models over concrete base space and fiber ambient space P2" begin
  @test parent(hypersurface_equation(h4)) == coordinate_ring(ambient_space(h4))
  @test dim(base_space(h4)) == d
  @test fiber_ambient_space(h4) == ambient_space_of_fiber_2
  @test is_smooth(fiber_ambient_space(h4)) == false
  @test is_simplicial(fiber_ambient_space(h4)) == true
  @test [string(g) for g in gens(cox_ring(fiber_ambient_space(h4)))] == ["x", "y", "z"]
  @test is_base_space_fully_specified(h4) == false
  @test is_partially_resolved(h4) == false
end

@testset "Error messages in hypersurface models over not fully specified base spaces" begin
  @test_throws ArgumentError hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, -1, ambient_space_of_fiber_2, [D1, D2, D3], p)
  @test_throws ArgumentError tune(h4, hypersurface_equation(h4))
end

h5 = literature_model(arxiv_id = "1208.2695", equation = "B.5")

@testset "Attributes and properties of hypersurface models over concrete base space and fiber ambient space P2" begin
  @test parent(hypersurface_equation(h5)) == coordinate_ring(ambient_space(h5))
  @test dim(base_space(h5)) == 2
  @test is_smooth(fiber_ambient_space(h5)) == false
  @test is_simplicial(fiber_ambient_space(h5)) == true
  @test [string(g) for g in gens(cox_ring(fiber_ambient_space(h5)))] == ["u", "w", "v"]
  @test is_base_space_fully_specified(h5) == false
  @test is_partially_resolved(h5) == false
end
