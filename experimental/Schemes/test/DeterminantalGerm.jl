@testset "DeterminantalGerm" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  L, _ = localization(complement_of_point_ideal(R, [0,0,0]))
  A = L[x 0 -z; 
        0 y  z]
  X_A = DeterminantalGerm(A, 2)
  @test determinantal_type(X_A) == (2,3,2)
  @test defining_matrix(X_A) === A
  @test Oscar._mat_type(X_A) === Val{:generic}
  @test dim(X_A) == 1
  @test codim(X_A) == 2
  @test point(X_A) == QQ.([0,0,0])
  @test modulus(OO(X_A)) == ideal(L, [x*y, x*z, y*z])
end

@testset "Symmetric DeterminantalGerm" begin
  R, (v,w,x,y,z) = QQ[:v,:w,:x,:y,:z]
  L, _ = localization(complement_of_point_ideal(R, [0,0,0,0,0]))
  A = L[v w x; 
        w x y; 
        x y z]
  X_A = DeterminantalGerm(A, 2, mat_type = :symmetric)
  @test determinantal_type(X_A) == (3,3,2)
  @test defining_matrix(X_A) === A
  @test Oscar._mat_type(X_A) === Val{:symmetric}
  @test dim(X_A) == 2
  @test codim(X_A) == 3
  @test point(X_A) == QQ.([0,0,0,0,0])
  @test modulus(OO(X_A)) == ideal(minors(L[v w x y; w x y z], 2))
end

@testset "Skew-symmetric DeterminantalGerm" begin
  R, (x,y) = QQ[:x,:y]
  L, _ = localization(complement_of_point_ideal(R, [0,0]))
  A = L[0  0 x y^2; 
        0  0 y^3 x; 
       -x -y^3 0 0; 
       -y^2 -x 0 0]
  X_A = DeterminantalGerm(A, 2, mat_type = :skew_symmetric)
  @test determinantal_type(X_A) == (4,4,2)
  @test defining_matrix(X_A) === A
  @test Oscar._mat_type(X_A) === Val{:skew_symmetric}
  @test dim(X_A) == 1
  @test codim(X_A) == 1
  @test point(X_A) == QQ.([0,0])
  @test modulus(OO(X_A)) == ideal(L, [x^2 - y^5])
end

@testset "DeterminantalGerm constructor errors" begin
  R, (x, y) = QQ[:x,:y]
  L, _ = localization(complement_of_point_ideal(R, [0,0]))
  A = L[0 0; 
        0 0]
  B = L[0 x;
        y 0]
  @test_throws ArgumentError DeterminantalGerm(A, 2, mat_type = :wrong_symbol)
  @test_throws ArgumentError DeterminantalGerm(A, 1000000)

  # matrix does not describe a singularity of expected codimension 
  @test_throws ErrorException DeterminantalGerm(A, 2)
  @test_throws ErrorException DeterminantalGerm(A, 2, mat_type = :symmetric)
  @test_throws ErrorException DeterminantalGerm(A, 1, mat_type = :skew_symmetric)
  
  # matrix does not have specified symmetry
  @test_throws ArgumentError DeterminantalGerm(B, 2, mat_type = :symmetric)
  @test_throws DomainError DeterminantalGerm(B, 1, mat_type = :skew_symmetric)
end

@testset "DeterminantalGerm MPoly-matrix constructor" begin
  R, (x,y) = QQ[:x,:y]
  A = R[0 x y;
        x y 0;
        y 0 x^2]
  X_A = DeterminantalGerm(A, 3, [0,0])
  X_A_sym = DeterminantalGerm(A, 3, [0,0], mat_type = :symmetric)
  @test X_A != X_A_sym
  @test defining_matrix(X_A) == A == defining_matrix(X_A_sym)
  @test determinantal_type(X_A) == (3,3,3) == determinantal_type(X_A_sym)
  @test Oscar.underlying_scheme(X_A) == Oscar.underlying_scheme(X_A_sym)

  B = R[0    0  x y^3;
        0    0  y x^2;
       -x   -y  0  0;
      -y^3 -x^2 0  0]
  X_B = DeterminantalGerm(B, 4, [0,0])
  X_B_skew = DeterminantalGerm(B, 2, [0,0], mat_type = :skew_symmetric)
  @test X_B != X_B_skew
  @test defining_matrix(X_B) == B == defining_matrix(X_B_skew)
  @test determinantal_type(X_B) == (4,4,4)
  @test determinantal_type(X_B_skew) == (4,4,2)
end



