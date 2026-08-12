using Oscar.RationalMapImage

@testset "RationalMapImage" begin

  @testset "generate_target_ring" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    L = [y*z, x*z, x*y]
    T = Oscar.RationalMapImage._generate_target_ring(L)
    @test ngens(T) == 3
    @test symbols(T) == [:b_0, :b_1, :b_2]
  end

  @testset "scheme_theoretic_image" begin
    # Dominant map: image is all of P^2
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    L = [y*z, x*z, x*y]
    T = Oscar.RationalMapImage._generate_target_ring(L)
    X = ideal(R, typeof(zero(R))[])
    img = Oscar.RationalMapImage._scheme_theoretic_image(L, X, T)
    @test is_zero(img)

    # Map from P^1 to P^2 (twisted cubic restricted)
    R2, (s, t) = polynomial_ring(QQ, [:s, :t])
    L2 = [s^2, s*t, t^2]
    T2 = Oscar.RationalMapImage._generate_target_ring(L2)
    X2 = ideal(R2, typeof(zero(R2))[])
    img2 = Oscar.RationalMapImage._scheme_theoretic_image(L2, X2, T2)
    # Image is the conic b_0*b_2 - b_1^2 = 0
    @test dim(img2) == 2
  end

  @testset "partial_image - surjective map" begin
    # Random cubics in P^2: map should be surjective (empty exceptional locus)
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    # Use fixed polynomials that give surjection onto P^2
    f1 = x^3 + y^2*z
    f2 = y^3 + x*z^2
    f3 = z^3 + x^2*y
    L = [f1, f2, f3]
    (img, exc) = partial_image(L)
    @test is_zero(img)
    @test isempty(exc)
  end

  @testset "partial_image - removable base locus 1" begin
    # L = [x^2, x*y]: base locus is {x=0} which is removable
    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    L = [x^2, x*y]
    (img, exc) = partial_image(L)
    @test !isempty(exc)
  end

  @testset "partial_image- removable base locus 2" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    L = [x*z, y*z, z^2]
    (img, exc) = partial_image(L)
    @test !isempty(exc)
  end

  @testset "total_image - simple non-trivial" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    L = [z^2, x^2, x*y]
    (N, E) = total_image(L; clean = true)
    @test length(E) == 2
    @test sort(dim.(N)) == [1, 2, 3]
  end

  @testset "total_image - Cremona transformation" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    L = [y*z, x*z, x*y]
    (N, E) = total_image(L; clean = true)
    @test !is_closed_image(L)
    @test sort(dim.(N)) == [1, 1, 1, 1, 1, 1, 2, 2, 2, 3]
  end

  @testset "total_image - projection of nodal cubic" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    X = ideal(R, [y^2*z - x^3 - x^2*z])
    L = [x, y]
    (N, E) = total_image(L, X; clean = true)
    @test sort(dim.(N)) == [1, 1, 2]
  end

  @testset "total_image - smooth cubic is closed" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    X = ideal(R, [y^2*z - x^3 - 2*x^2*z - z^3])
    L = [x, y - z]
    @test is_closed_image(L, X)
  end

  @testset "total_image - conic projection misses infinity" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    X = ideal(R, [x*y - z^2])
    L = [y, z]
    (N, E) = total_image(L, X; clean = true)
    @test sort(dim.(N)) == [1, 2]
  end

  @testset "total_image - affine map" begin
    R, (u, v, t) = polynomial_ring(QQ, [:u, :v, :t])
    L = [u*v, u*t, v^2, t^2]
    (N, E) = total_image(L; affine = true, clean = true)
    @test sort(dim.(N)) == [2, 3]
  end

  @testset "total_image - singular cone projection" begin
    R, (x, y, z, w) = polynomial_ring(QQ, [:x, :y, :z, :w])
    X = ideal(R, [x*y - z^2])
    L = [x, y, z]
    (N, E) = total_image(L, X; clean = true)
    @test sort(dim.(N)) == [2]
  end

  @testset "is_closed_image" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    # Cremona is NOT closed
    @test !is_closed_image([y*z, x*z, x*y])
    # Projection of smooth cubic onto P^1 IS closed
    X = ideal(R, [y^2*z - x^3 - 2*x^2*z - z^3])
    @test is_closed_image([x, y - z], X)
  end

  @testset "ConstructibleSetTree - type and iteration" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    L = [z^2, x^2, x*y]
    C = total_image(L; clean = true)
    @test C isa ConstructibleSetTree
    # Tuple-style destructuring must still work.
    (N, E) = C
    @test N === C.nodes
    @test E === C.edges
    @test length(C) == 2
  end

  @testset "is_contained - conic projection" begin
    # Map (x:y:z) -> (y:z) on the conic xy = z^2.
    # Image = P^1 \ {[0:1]}: every (b_0:b_1) with b_0 != 0 is achieved;
    # (0:1) requires y=0, z!=0 on the conic, which gives 0 = z^2 -> z=0, contradiction.
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    X = ideal(R, [x*y - z^2])
    L = [y, z]
    C = total_image(L, X; clean = true)
    @test sort(dim.(C.nodes)) == [1, 2]  # sanity check
    # [1:0] is in the image: source [0:1:0] -> [1:0].
    @test is_contained([1, 0], C)
    # [1:1] is in the image: source [1:1:1] satisfies 1*1=1^2.
    @test is_contained([1, 1], C)
    # [0:1] is NOT in the image.
    @test !is_contained([0, 1], C)
  end

  @testset "is_closed_image - classical embeddings" begin
    # Veronese embedding nu_2: P^1 -> P^2, [s:t] |-> [s^2:st:t^2].
    # This is a closed immersion; image is the smooth conic b_0*b_2 - b_1^2 = 0.
    R1, (s, t) = polynomial_ring(QQ, [:s, :t])
    @test is_closed_image([s^2, s*t, t^2])

    # Twisted cubic embedding nu_3: P^1 -> P^3, [s:t] |-> [s^3:s^2*t:s*t^2:t^3].
    # Also a closed immersion (rational normal curve of degree 3).
    R2, (u, v) = polynomial_ring(QQ, [:u, :v])
    @test is_closed_image([u^3, u^2*v, u*v^2, v^3])
  end

  @testset "is_contained - Veronese and twisted cubic" begin
    # Veronese conic in P^2.
    R1, (s, t) = polynomial_ring(QQ, [:s, :t])
    C2 = total_image([s^2, s*t, t^2]; clean = true)
    # [1:1] maps to [1:1:1] on the conic.
    @test is_contained([1, 1, 1], C2)
    # [1:0:1] is not on b_0*b_2 - b_1^2 = 0.
    @test !is_contained([1, 0, 1], C2)

    # Twisted cubic in P^3.
    R2, (u, v) = polynomial_ring(QQ, [:u, :v])
    C3 = total_image([u^3, u^2*v, u*v^2, v^3]; clean = true)
    # [1:1] maps to [1:1:1:1].
    @test is_contained([1, 1, 1, 1], C3)

    @test !is_contained([1, 0, 1, 0], C3)
  end

end
