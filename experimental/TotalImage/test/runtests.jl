# Test projective, recursive, closed, and affine images translated from TotalImage.m2.
using Oscar.TotalImage

@testset "Total image of a rational map" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  cremona = [y*z, x*z, x*y]

  @testset "target ring and scheme-theoretic image" begin
    target = Oscar.TotalImage._target(cremona)
    @test ngens(target) == 3
    source = ideal(R, [zero(R)])
    @test is_zero(Oscar.TotalImage._scheme_theoretic_image(cremona, source, target))

    R2, (s, t) = polynomial_ring(QQ, [:s, :t])
    veronese = [s^2, s*t, t^2]
    target2 = Oscar.TotalImage._target(veronese)
    b0, b1, b2 = gens(target2)
    image2 = Oscar.TotalImage._scheme_theoretic_image(
      veronese, ideal(R2, [zero(R2)]), target2
    )
    @test image2 == ideal(target2, [b0*b2 - b1^2])
  end

  image, exceptional_loci = partial_image(cremona)
  @test is_zero(image)
  @test length(exceptional_loci) == 3
  @test all(I -> dim(I) == 2, exceptional_loci)

  T, _ = polynomial_ring(QQ, [:p, :q, :r])
  custom_image, custom_exception = partial_image(
    cremona, ideal(R, [zero(R)]), T; decompose=false
  )
  @test base_ring(custom_image) === T
  @test dim(custom_exception) == 2

  tree = total_image(cremona)
  @test tree isa ConstructibleSetTree
  @test length(tree) == 10
  @test length(edges(tree)) == 9
  @test sort(dim.(vertices(tree))) == [1, 1, 1, 1, 1, 1, 2, 2, 2, 3]
  @test !is_closed_image(cremona)
  @test is_contained([1, 1, 1], tree)
  @test !is_contained([1, 1, 0], tree)
  @test is_contained([1, 0, 0], tree)
  @test sprint(show, tree) isa String
  @test sprint(show, "text/plain", tree) isa String
  @test_throws ArgumentError is_contained([0, 0, 0], tree)
  @test_throws ArgumentError is_contained([1, 1], tree)

  @testset "surjective cubics" begin
    cubics = [x^3 + y^2*z, y^3 + x*z^2, z^3 + x^2*y]
    cubic_image, cubic_exception = partial_image(cubics)
    @test is_zero(cubic_image)
    @test isempty(cubic_exception)
  end

  removable = total_image([x^2, x*y])
  @test length(edges(removable)) == 1
  @test sort(dim.(vertices(removable))) == [1, 2]
  _, removable_exception = partial_image([x^2, x*y])
  @test !isempty(removable_exception)

  _, removable_exception2 = partial_image([x*z, y*z, z^2])
  @test !isempty(removable_exception2)

  @testset "tree accessors" begin
    simple = total_image([z^2, x^2, x*y])
    @test simple isa ConstructibleSetTree
    @test length(simple) == 3
    @test length(edges(simple)) == 2
    @test sort(dim.(vertices(simple))) == [1, 2, 3]
    @test vertices(simple) == Oscar.TotalImage.image_ideals(simple)
  end

  cubic = ideal(R, [y^2*z - x^3 - 2*x^2*z - z^3])
  @test is_closed_image([x, y - z], cubic)

  nodal_cubic = ideal(R, [y^2*z - x^3 - x^2*z])
  nodal_image = total_image([x, y], nodal_cubic)
  @test sort(dim.(vertices(nodal_image))) == [1, 1, 2]

  conic = ideal(R, [x*y - z^2])
  conic_projection = total_image([y, z], conic)
  @test sort(dim.(vertices(conic_projection))) == [1, 2]
  @test is_contained([1, 0], conic_projection)
  @test is_contained([1, 1], conic_projection)
  @test !is_contained([0, 1], conic_projection)

  S, (a, b, c, d) = polynomial_ring(QQ, [:a, :b, :c, :d])
  cone_image = total_image([a, b, c], ideal(S, [a*b - c^2]))
  @test isempty(edges(cone_image))
  @test dim(only(vertices(cone_image))) == 2
  @test is_contained([1, 1, 1], cone_image)
  @test !is_contained([1, 1, 0], cone_image)

  A, (u, v) = polynomial_ring(QQ, [:u, :v])
  hyperbola_image = total_image([u], ideal(A, [u*v - 1]))
  @test is_affine(hyperbola_image)
  @test length(edges(hyperbola_image)) == 1
  @test sort(dim.(vertices(hyperbola_image))) == [0, 1]
  @test is_contained([2], hyperbola_image)
  @test !is_contained([0], hyperbola_image)

  @testset "affine polynomial map" begin
    B, (r, s, t) = polynomial_ring(QQ, [:r, :s, :t])
    affine_image = total_image([r*s, r*t, s^2, t^2]; affine=true)
    @test sort(dim.(vertices(affine_image))) == [0, 2, 3]
    @test is_contained([1, 1, 1, 1], affine_image)
    @test is_contained([0, 0, 0, 0], affine_image)
    @test !is_contained([1, 0, 0, 0], affine_image)
  end

  @testset "classical embeddings" begin
    P1, (s, t) = polynomial_ring(QQ, [:s, :t])
    veronese = [s^2, s*t, t^2]
    @test is_closed_image(veronese)
    veronese_image = total_image(veronese)
    @test is_contained([1, 1, 1], veronese_image)
    @test !is_contained([1, 0, 1], veronese_image)

    Q1, (u, v) = polynomial_ring(QQ, [:u, :v])
    twisted_cubic = [u^3, u^2*v, u*v^2, v^3]
    @test is_closed_image(twisted_cubic)
    twisted_cubic_image = total_image(twisted_cubic)
    @test is_contained([1, 1, 1, 1], twisted_cubic_image)
    @test !is_contained([1, 0, 1, 0], twisted_cubic_image)
  end
end
