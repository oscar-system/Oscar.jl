@testset "Binomial Ideals" begin

  @testset "Binomial and unital test" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
    f = x+y
    @test isbinomial(f)
    J = ideal(R, [x^2-y^3, z^2])
    @test isbinomial(J)
    Qxy, (x, y, z, t) = PolynomialRing(FlintQQ, 4)
    I = ideal(elem_type(Qxy)[x*y, z*t^2-t^3, z^2-y^2])
    @test Oscar.isbinomial(I)
    @test Oscar.isunital(I)

    J = ideal([x*y - z*t^2 + t^3, z*t^2-t^3])
    @test Oscar.isbinomial(J)
    J1 = ideal(Qxy, x^2+y^2+z^2)
    @test !Oscar.isbinomial(J1)
  end
  
  @testset "Cellular decomposition" begin
  
    R, x = PolynomialRing(QQ, "x"=>1:5)
    I = ideal(R, [x[1]^3*x[3]-x[1]^3, x[1]^4, x[1]^2*x[2]*x[4]-x[1]^2*x[2], x[2]^2, x[4]^3-1])
    @test iscellular(I)[1]
    I = ideal(R, [x[1]*x[4]^2-x[2]*x[5]^2, x[1]^3*x[3]^3-x[2]^4*x[4]^2, x[2]*x[4]^8-x[3]^3*x[5]^6])
    @test !iscellular(I)[1]
    Qxy, (x, y, z, t) = PolynomialRing(FlintQQ, 4)
    I = ideal(elem_type(Qxy)[x*y, z*t^2-t^3, z^2-y^2])
    @test !Oscar.iscellular(I)[1]
    lI = Oscar.cellular_decomposition(I)
    lI2 = Oscar.cellular_decomposition_macaulay(I)
    @test length(lI) == length(lI2)
    for x in lI
      @test iscellular(x)[1]
      @test x in lI2
    end

  end

  @testset "Binomial primary decomposition" begin
    Qxy, (x, y, z, t) = PolynomialRing(FlintQQ, 4)
    I = ideal(elem_type(Qxy)[x*y, z*t^2-t^3, z^2-y^2])
    lP = Oscar.primary_decomposition(I)
    lP1 = Oscar.binomial_primary_decomposition(I)
    @test length(lP) == length(lP1)
    for x in lP
      @test x in lP1
    end

    I = Oscar.birth_death_ideal(2, 1)
    lP = Oscar.primary_decomposition(I)
    lP1 = Oscar.binomial_primary_decomposition(I)
    @test length(lP) == length(lP1)
    for x in lP
      @test x in lP1
    end
  end

end