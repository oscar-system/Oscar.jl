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

  @testset "partial_image - removable base locus 2" begin
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

end
