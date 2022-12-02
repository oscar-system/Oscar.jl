@testset "singular locus and is_smooth" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  X = Spec(R, I*J, units_of(R))
  Z, _ = singular_locus(X)
  @test is_subset(subscheme(X, [x,y,z]), Z)
  @test is_subset(subscheme(X, [z^2-3, y-2, x-1]), Z)
  @test !is_smooth(X)
  @test is_smooth(Spec(R))
  Y = Spec(R, ideal(R, [x^2 - y^2 + z^2 - 1]))
  @test is_smooth(Y)
  U1 = complement_of_ideal(R, [0, 0, 0])
  X1 = Spec(R, I * J, U1)
  Z1, _ = singular_locus(X1)
  @test is_subset(subscheme(X1,  [x, y, z]), Z1)
  @test is_smooth(Spec(R, U1))
  U2 = complement_of_ideal(R, [1, 1, 0])
  X2 = Spec(R, I * J, U2)
  Z2, _ = singular_locus(X2)
  @test is_empty(Z2)
  @test is_smooth(X2)
  X3 = Spec(R, ideal(R,[x,y,z]), units_of(R))
  is_smooth(X3)
  @test singular_locus(Spec(R))[1] == Spec(R,ideal(R,[one(R)]))
end

@testset "singular_locus_reduced" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [(x^2 + y^2 + z^2)*(x-y)^2])
  X = Spec(R, I, units_of(R))
  Y, _ = singular_locus(X)
  Yr, _ = singular_locus_reduced(X)
  ISL1 = ideal(R, [x-y])
  ISL2 = ideal(R, [x-y,2y^2+z^2])
  SL1 = Spec(R, ISL1,units_of(R))
  SL2= Spec(R,ISL2,units_of(R))
  @test reduced_scheme(Y)[1] == SL1
  @test reduced_scheme(Yr)[1] == SL2
end

@testset "is_equidimensional, is_reduced etc" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [x^2 - y^2 + z^2])
  I2 = ideal(R, [(x^2 - y^2 + z^2)^2])
  J = ideal(R, [x-1, y-2])
  X = Spec(R, I , units_of(R))
  X2 = Spec(R, I2 , units_of(R))
  Y = Spec(R, I*J, units_of(R))
  @test is_equidimensional(X)
  @test is_equidimensional(X)   ## if caching works, this tests a different line of code
  @test is_equidimensional(X2)
  @test !is_equidimensional(Y)
  @test is_reduced(X)
  @test !is_reduced(X2)
  @test reduced_scheme(X2)[1] == X
  U = complement_of_ideal(R,[0,0,0])
  X1 = Spec(R, I ,U)
  X12 = Spec(R, I2, U)
  Y1 = Spec(R, I*J, U)
  @test is_equidimensional(X1)
  @test is_equidimensional(X12)
  @test is_equidimensional(Y1) # pther component does not pass through point
  @test is_reduced(X1)
  @test !is_reduced(X12)
  @test reduced_scheme(X12)[1] == X1
  @test is_equidimensional(Spec(R))
  @test is_equidimensional(Spec(R,U))
  @test is_reduced(Spec(R))
  @test is_reduced(Spec(R,U))
  Z = Spec(R)
  @test reduced_scheme(Z)[1] == Z
end

@testset "derivative in localized ring" begin
  R, (x,y) = QQ["x","y"]
  U = complement_of_ideal(R,[0,0])
  X = Spec(R, U)
  L = OO(X)
  I = ideal(L, [1//(1+x), (1-y)^2//((y-1)*(1-x))])
  @test denominator(derivative(gens(I)[1],1)) == (x^2 + 2*x + 1)
  @test numerator(derivative(gens(I)[1],1)) == -1
  @test denominator(derivative(gens(I)[2],1)) == (x^2 - 2*x + 1)
  @test numerator(derivative(gens(I)[2],1)) == y-1
  @test denominator(derivative(gens(I)[2],2)) == (x - 1)
  @test numerator(derivative(gens(I)[2],2)) == -1
  J=ideal(R, [x^2+y^2])
  X2 = Spec(R, J ,U)
  L2 = OO(X2)
  I2 = ideal(L2, [1//(1+x), (1-y)^2//((y-1)*(1-x))])
  @test denominator(derivative(gens(I2)[1],1)) == (x^2 + 2*x + 1)
  @test numerator(derivative(gens(I2)[1],1)) == -1
  @test denominator(derivative(gens(I2)[2],1)) == (x^2 - 2*x + 1)
  @test numerator(derivative(gens(I2)[2],1)) == y-1
  @test denominator(derivative(gens(I2)[2],2)) == (x - 1)
  @test numerator(derivative(gens(I2)[2],2)) == -1
end


