# Check that chosen Cox equations yield certified geometric presentations.
@testset "Toric hypersurface and complete-intersection presentations" begin
  X = projective_space(NormalToricVariety, 2)
  x, y, z = gens(cox_ring(X))
  H = hypersurface_in_toric_variety(X, x^2 - y * z)
  @test H isa HypersurfaceInToricVariety
  @test H isa AbsCompleteIntersectionInToricVariety
  @test H isa AbsCoveredScheme{QQField}
  @test sprint(show, MIME("text/plain"), H) ==
    "Hypersurface in a normal toric variety"
  @test Oscar.underlying_scheme(H) === covered_scheme(H)
  @test base_ring(H) === base_ring(covered_scheme(H))
  @test default_covering(H) === default_covering(covered_scheme(H))
  @test length(affine_charts(H)) == length(affine_charts(covered_scheme(H)))
  @test has_expected_codimension(H)
  @test is_smooth(H)
  @test is_empty(H) == is_empty(covered_scheme(H))
  @test is_irreducible(H) == is_irreducible(covered_scheme(H))
  @test is_integral(H) == is_integral(covered_scheme(H))
  @test defining_equation(H) == x^2 - y * z
  @test defining_divisor_class(H) == toric_divisor_class(X, degree(x^2))
  @test defining_divisor_classes(H) == [toric_divisor_class(X, degree(x^2))]
  @test fundamental_class(H) ==
    rational_equivalence_class(only(defining_divisor_classes(H)))
  @test fundamental_class(H) == rational_equivalence_class(H)
  @test covered_scheme(H) === covered_scheme(H)
  inc = Oscar.closed_embedding(H)
  @test inc isa Oscar.CoveredClosedEmbedding
  @test domain(inc) === covered_scheme(H)
  @test codomain(inc) === X
  @test image_ideal(inc) == ideal_sheaf(X, defining_ideal(H))
  @test defining_ideal(H) == ideal(cox_ring(X), x^2 - y * z)

  C = complete_intersection_in_toric_variety(X, [x, y])
  @test C isa CompleteIntersectionInToricVariety
  @test C isa AbsCoveredScheme{QQField}
  @test Oscar.underlying_scheme(C) === covered_scheme(C)
  @test is_reduced(C) == is_reduced(covered_scheme(C))
  @test is_empty(C) == is_empty(covered_scheme(C))
  @test length(defining_equations(C)) == 2
  C2 = complete_intersection_in_toric_variety(X, [x, x + y])
  @test defining_ideal(C) == defining_ideal(C2)
  @test defining_equations(C) != defining_equations(C2)
  @test has_expected_codimension(C)
  @test dim(C) == 0
  @test_throws ArgumentError complete_intersection_in_toric_variety(X, [x, x])
  @test_throws MethodError CompleteIntersectionInToricVariety(
    X, [x], [toric_divisor_class(X, degree(x^3))]
  )
  @test_throws ArgumentError hypersurface_in_toric_variety(X, x + y^2)
  @test_throws ArgumentError hypersurface_in_toric_variety(X, zero(cox_ring(X)))

  singular = hypersurface_in_toric_variety(X, x * y)
  @test !is_smooth(singular)
  @test is_smooth(covered_scheme(H))
end

@testset "Toric complete-intersection input validation" begin
  X = projective_space(NormalToricVariety, 2)
  Y = projective_space(NormalToricVariety, 1)
  x, y, z = gens(cox_ring(X))
  u = first(gens(cox_ring(Y)))

  @test_throws ArgumentError hypersurface_in_toric_variety(X, u)
  @test_throws ArgumentError complete_intersection_in_toric_variety(X, [x, x])
  @test_throws ArgumentError complete_intersection_in_toric_variety(X, [x, y, z])
  @test_throws ArgumentError complete_intersection_in_toric_variety(X, typeof(x)[])
  @test_throws ArgumentError hypersurface_in_toric_variety(X, x + y^2)
  @test_throws ArgumentError hypersurface_in_toric_variety(X, zero(cox_ring(X)))

  # Skipping a certification check must not manufacture a positive result.
  Z = complete_intersection_in_toric_variety(
    X, [x, x]; check_codimension=false
  )
  @test !has_expected_codimension(Z)
  @test !is_empty(Z) # The codimension is wrong even though the scheme is nonempty.
  @test_throws ArgumentError complete_intersection_in_toric_variety(
    X, [x, x]; check_codimension=true
  )

  # Empty zero loci are valid presentations, including with certification.
  P = normal_toric_variety(cube(2))
  x1, x2, y1, y2 = gens(cox_ring(P))
  E = complete_intersection_in_toric_variety(P, [x1, x2])
  @test is_empty(E)
  @test is_empty(E) == is_empty(covered_scheme(E))
  @test has_expected_codimension(E)
  @test is_trivial(fundamental_class(E))

  # A nowhere-vanishing constant section defines an empty hypersurface.
  Hempty = hypersurface_in_toric_variety(P, one(cox_ring(P)))
  @test is_empty(Hempty)
  @test has_expected_codimension(Hempty)
  @test is_trivial(fundamental_class(Hempty))

  # The supplied equation and its degree determine the presentation.
  H = hypersurface_in_toric_variety(X, x * y)
  @test defining_equation(H) == x * y
  @test !is_smooth(H)
end

@testset "Toric equation accessors do not expose cached polynomials" begin
  X = projective_space(NormalToricVariety, 2)
  x, y, z = gens(cox_ring(X))
  H = hypersurface_in_toric_variety(X, x)
  class_before = fundamental_class(H)
  @test is_reduced(H)

  # Mutate polynomial elements, not merely the returned vector.
  fs = defining_equations(H)
  set_exponent_vector!(fs[1], 1, [2, 0, 0])
  @test defining_equation(H) == x
  f = defining_equation(H)
  set_exponent_vector!(f, 1, [2, 0, 0])
  @test defining_equation(H) == x
  @test defining_ideal(H) == ideal(cox_ring(X), x)
  @test fundamental_class(H) == class_before
  @test is_reduced(H)

  C = complete_intersection_in_toric_variety(X, [x, y])
  gs = defining_equations(C)
  set_exponent_vector!(gs[1], 1, [2, 0, 0])
  @test defining_equations(C) == [x, y]
  @test defining_ideal(C) == ideal(cox_ring(X), [x, y])
end

@testset "Toric complete-intersection ambient fan validation" begin
  # A pure fan of one-dimensional rays in ambient dimension two.
  X = normal_toric_variety(
    incidence_matrix([[1], [2], [3]]),
    [[1, 0], [0, 1], [-1, -1]],
  )
  x = first(gens(cox_ring(X)))
  @test_throws ArgumentError hypersurface_in_toric_variety(
    X, x; check_codimension=false
  )
  @test_throws ArgumentError complete_intersection_in_toric_variety(
    X, [x]; check_codimension=false
  )
  @test_throws ArgumentError hypersurface_in_toric_variety(X, x)

  # Full-dimensional cones alone are insufficient if the fan is not pure.
  Y = normal_toric_variety(
    incidence_matrix([[1, 2], [3]]),
    [[1, 0], [0, 1], [-1, -1]],
  )
  y = first(gens(cox_ring(Y)))
  @test_throws ArgumentError complete_intersection_in_toric_variety(
    Y, [y]; check_codimension=false
  )
end

@testset "Toric hypersurface reducedness" begin
  X = projective_space(NormalToricVariety, 2)
  x = first(gens(cox_ring(X)))
  H = hypersurface_in_toric_variety(X, x^2)
  @test H isa HypersurfaceInToricVariety
  @test H isa AbsCoveredScheme{QQField}
  @test !is_reduced(H)
  @test has_expected_codimension(H)
  @test defining_divisor_class(H) == toric_divisor_class(X, degree(x^2))
  @test fundamental_class(H) == rational_equivalence_class(defining_divisor_class(H))
  @test dim(H) == dim(covered_scheme(H))
  @test is_smooth(H) == is_smooth(covered_scheme(H))
end
