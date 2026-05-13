@testset "DeterminantalGerm" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  L, _ = localization(complement_of_point_ideal(R, [0,0,0]))
  A = L[x 0 -z; 
        0 y  z]
  X_A = DeterminantalGerm(A, 2)
  @test determinantal_type(X_A) == (2,3,2)
  @test defining_matrix(X_A) === A
  @test dim(X_A) == 1
  @test codim(X_A) == 2
  @test point(X_A) == QQ.([0,0,0])
  @test modulus(OO(X_A)) == ideal(L, [x*y, x*z, y*z])
end

@testset "SymmetricDeterminantalGerm" begin
  R, (v,w,x,y,z) = QQ[:v,:w,:x,:y,:z]
  L, _ = localization(complement_of_point_ideal(R, [0,0,0,0,0]))
  A = L[v w x; 
        w x y; 
        x y z]
  X_A = SymmetricDeterminantalGerm(A, 2)
  @test determinantal_type(X_A) == (3,3,2)
  @test defining_matrix(X_A) === A
  @test dim(X_A) == 2
  @test codim(X_A) == 3
  @test point(X_A) == QQ.([0,0,0,0,0])
end

@testset "SkewSymmetricDeterminantalGerm" begin
  R, (x,y) = QQ[:x,:y]
  L, _ = localization(complement_of_point_ideal(R, [0,0]))
  A = L[0  0 x y^3; 
        0  0 y^2 x; 
       -x -y^2 0 0; 
       -y^3 -x 0 0]
  X_A = SkewSymmetricDeterminantalGerm(A, 2)
  @test determinantal_type(X_A) == (4,4,2)
  @test defining_matrix(X_A) === A
  @test dim(X_A) == 1
  @test codim(X_A) == 1
  @test point(X_A) == QQ.([0,0])
  @test modulus(OO(X_A)) == ideal(L, [x^2 - y^5])
end