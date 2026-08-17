@testset "all tests - ActionMap.jl" verbose=true begin

  @testset "ZZ and QQ" begin
    s_ZZ = action_shift(ZZ)
    d_ZZ = action_derivation(ZZ)
    
    s_QQ = action_shift(QQ)
    d_QQ = action_derivation(QQ)

    @test_throws ArgumentError action_derivation(s_ZZ)
    @test_throws ArgumentError action_shift(d_ZZ)

    @test s_ZZ isa Oscar.TrivialActionShift
    @test s_QQ isa Oscar.TrivialActionShift
    @test d_ZZ isa Oscar.TrivialActionDerivation
    @test d_QQ isa Oscar.TrivialActionDerivation
    
    for am in [s_ZZ, d_ZZ]
      @test domain(am) === ZZ
      @test codomain(am) === ZZ
    end
    
    @test s_ZZ(-7) == -7
    @test s_ZZ(ZZ(5)) == ZZ(5)
    @test d_ZZ(17) == 0
    @test d_ZZ(ZZ(-4)) == ZZ(0)

    @test s_QQ(-4//9) == -4//9
    @test s_QQ(QQ(3//8)) == QQ(3//8)
    @test s_QQ(6) == 6
    @test d_QQ(3//4) == 0
    @test d_QQ(-4) == 0
    @test d_QQ(QQ(4//9)) == 0
  end

  @testset "Univariate polynomial rings" begin
    R, x = polynomial_ring(QQ, :x)
    
    zero_R = action_derivation(map_from_func(R, R, p -> zero(R)))
    @test zero_R isa Oscar.TrivialActionDerivation

    ddx = action_derivation(map_from_func(R, R, derivative))
    @test ddx isa Oscar.NontrivialActionDerivation
    @test domain(ddx) === R
    @test codomain(ddx) === R
    @test ddx(x^3) == 3*x^2
    @test ddx(17*x) == 17

    id_R = action_shift(map_from_func(R, R, p -> p))
    @test id_R isa Oscar.TrivialActionShift
    @test id_R(x^2 + 1) == x^2 + 1

    x_to_zero = action_shift(hom(R, R, zero(R))) # evaluation at zero
    @test x_to_zero isa Oscar.NontrivialActionShift
    @test domain(x_to_zero) === R
    @test codomain(x_to_zero) === R
    @test x_to_zero(x^2 + 1) == 1 
  end

  @testset "Multivariate polynomial rings" begin
    R, (x, y) = polynomial_ring(ZZ, [:x, :y])
    
    swap_hom = hom(R, R, [y, x])
    swap = action_shift(swap_hom)
    
    @test swap isa Oscar.NontrivialActionShift
    @test swap(x^2 * y) == x * y^2
    @test swap(x^2 + y) == x + y^2
    
    mmf_ddx = map_from_func(R, R, p -> derivative(p, x))
    ddx = action_derivation(mmf_ddx)
    
    @test ddx isa Oscar.NontrivialActionDerivation
    @test ddx(x^3 * y^2) == 3*x^2*y^2
    @test ddx(y^2 + 2*y + 1) == 0

    @test_throws ArgumentError action_derivation(swap_hom)
    @test_throws ArgumentError action_shift(mmf_ddx)
  end
  @testset "Rational function fields" begin
    R, x = polynomial_ring(QQ, :x)
    F = fraction_field(R)
    
    ddx_F = action_derivation(map_from_func(F, F, derivative))
    
    @test ddx_F isa Oscar.NontrivialActionDerivation
    @test domain(ddx_F) === F
    @test codomain(ddx_F) === F
    
    @test ddx_F(F(1) // F(x)) == F(-1) // F(x^2) # Quotient rule

    g = F(1) // F(x)^2 
    
    shift_func = r -> F(evaluate(numerator(r), g)) // F(evaluate(denominator(r), g)) # p -> p(x^(-2))
    inv_sq_shift = action_shift(map_from_func(F, F, shift_func))

    @test inv_sq_shift isa Oscar.NontrivialActionShift
    @test domain(inv_sq_shift) === F
    @test codomain(inv_sq_shift) === F
    
    @test inv_sq_shift(F(x)) == g
    @test inv_sq_shift((F(x^2 + 1)) // F(x)) == (F(1 + x^4) // F(x^2))
  end
  @testset "Quadratic number fields" begin
    for d in [-5, -3, -1, 2, 3, 5]
      K, a = quadratic_field(d)
    
      trivial_galois = action_shift(hom(K, K, a))
      nontrivial_galois = action_shift(hom(K, K, -a))
    
      @test trivial_galois isa Oscar.TrivialActionShift
      @test nontrivial_galois isa Oscar.NontrivialActionShift
      @test trivial_galois(a) == a
      @test nontrivial_galois(a) == -a
      @test trivial_galois(a) + nontrivial_galois(a) == trace(a)
      @test trivial_galois(a) * nontrivial_galois(a) == norm(a)
      @test nontrivial_galois(K(3) + 2*a) == K(3) - 2*a
    end
  end
end
