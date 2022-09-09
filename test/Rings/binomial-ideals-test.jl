@testset "Binomial Ideals" begin
  
  @testset "Binomial and unital test" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
    f = x+y
    @test is_binomial(f)
    J = ideal(R, [x^2-y^3, z^2])
    @test is_binomial(J)
    Qxy, (x, y, z, t) = PolynomialRing(FlintQQ, 4)
    I = ideal(elem_type(Qxy)[x*y, z*t^2-t^3, z^2-y^2])
    @test Oscar.is_binomial(I)
    @test Oscar.is_unital(I)

    J = ideal([x*y - z*t^2 + t^3, z*t^2-t^3])
    @test Oscar.is_binomial(J)
    J1 = ideal(Qxy, x^2+y^2+z^2)
    @test !Oscar.is_binomial(J1)
  end
  
  @testset "Cellular decomposition" begin
  
    R, x = PolynomialRing(QQ, "x"=>1:5)
    I = ideal(R, [x[1]^3*x[3]-x[1]^3, x[1]^4, x[1]^2*x[2]*x[4]-x[1]^2*x[2], x[2]^2, x[4]^3-1])
    @test is_cellular(I)[1]
    I = ideal(R, [x[1]*x[4]^2-x[2]*x[5]^2, x[1]^3*x[3]^3-x[2]^4*x[4]^2, x[2]*x[4]^8-x[3]^3*x[5]^6])
    @test !is_cellular(I)[1]
    Qxy, (x, y, z, t) = PolynomialRing(FlintQQ, 4)
    I = ideal(elem_type(Qxy)[x*y, z*t^2-t^3, z^2-y^2])
    @test !Oscar.is_cellular(I)[1]
    lI = Oscar.cellular_decomposition(I)
    lI2 = Oscar.cellular_decomposition_macaulay(I)
    @test length(lI) == length(lI2)
    for x in lI
      @test is_cellular(x)[1]
      @test x in lI2
    end
    R, x = PolynomialRing(QQ, "x"=>1:3)
    I = ideal(R, [x[3]^2*(x[1]^2-x[2]^2), x[3]*(x[1]^4-x[2]^4), x[3]^3])
    ap = cellular_minimal_associated_primes(I)
    @test length(ap) == 1
    RQQAb = base_ring(ap[1])
    @test ap[1] == ideal(RQQAb, [RQQAb[3]])

  end

  @testset "Binomial primary decomposition" begin
    Qxy, (x, y, z, t) = PolynomialRing(FlintQQ, 4)
    I = ideal(elem_type(Qxy)[x*y, z*t^2-t^3, z^2-y^2])
    lP = Oscar.binomial_primary_decomposition(I)
    
    J = lP[1][1]
    for i = 2:length(lP)
      J = intersect(J, lP[i][1])
    end
    RQQAb = base_ring(J)
    QQAb = base_ring(RQQAb)
    x, y, z, t = gens(RQQAb) 
    IQQAb = ideal([x*y, z*t^2-t^3, z^2-y^2])
    @test IQQAb == J
    

    I = Oscar.birth_death_ideal(2, 1)
    lP = Oscar.primary_decomposition(I)
    lP1 = Oscar.binomial_primary_decomposition(I)
    @test length(lP) == length(lP1)
    RQQAb = base_ring(lP1[1][1])
    QQAb = base_ring(RQQAb)
    for x in lP
      y = ideal(RQQAb, [map_coefficients(QQAb, p, parent = RQQAb) for p in gens(x[1])])
      z = ideal(RQQAb, [map_coefficients(QQAb, p, parent = RQQAb) for p in gens(x[2])])
      @test (y, z) in lP1
    end
    
    R, x = PolynomialRing(QQ, "x"=>1:3)
    I = ideal(R, [x[3]^2*(x[1]^2-x[2]^2), x[3]*(x[1]^4-x[2]^4), x[3]^3])
    @test length(Oscar.binomial_primary_decomposition(I)) == 5

  end

end