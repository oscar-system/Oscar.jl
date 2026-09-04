@testset "DeterminantalGerm" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  L, _ = localization(complement_of_point_ideal(R, [2,1,0]))
  A = L[x-2 0 -z; 
        0 y-1  z]
  X_A = DeterminantalGerm(A, 2)
  @test determinantal_type(X_A) == (2,3,2)
  @test determinantal_ideal(X_A) == ideal(L, [(x-2)*(y-1), (x-2)*z, (y-1)*z])
  @test defining_matrix(X_A) === A
  @test Oscar._mat_type(X_A) === Val{:generic}
  @test dim(X_A) == 1
  @test codim(X_A) == 2
  @test point(X_A) == QQ.([2,1,0])
  @test modulus(OO(X_A)) == ideal(L, [(x-2)*(y-1), (x-2)*z, (y-1)*z])
end

@testset "Symmetric DeterminantalGerm (Pinkham's example)" begin
  R, (v,w,x,y,z) = QQ[:v,:w,:x,:y,:z]
  L, _ = localization(complement_of_point_ideal(R, [23,0,0,0,0]))
  A = L[v-23 w  x; 
         w   x  y; 
         x   y  z]
  X_A = DeterminantalGerm(A, 2, mat_type = :symmetric)
  @test determinantal_type(X_A) == (3,3,2)
  @test determinantal_ideal(X_A) == ideal(minors(L[v-23 w x y; w x y z], 2))
  @test defining_matrix(X_A) === A
  @test Oscar._mat_type(X_A) === Val{:symmetric}
  @test dim(X_A) == 2
  @test codim(X_A) == 3
  @test point(X_A) == QQ.([23,0,0,0,0])
  @test modulus(OO(X_A)) == ideal(minors(L[v-23 w x y; w x y z], 2))
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
  @test determinantal_ideal(X_A) == ideal(L, [x^2 - y^5])
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

@testset "Endowing a determinantal structure onto a scheme" begin
  R, (x,y,z) = QQ[:x,:y,:z]
  I = ideal([x*y, x*z, y*z])
  U = complement_of_point_ideal(R, [0,0,0])
  X = spec(R, I, U)
  A = localized_ring(OO(X))[x 0 z;
                            0 y z]
  X_A = DeterminantalGerm(X, A, 2)
  @test defining_ideal(X_A) == localized_ring(OO(X))(I)
  @test determinantal_ideal(X_A) == ideal(localized_ring(OO(X)), [x*y, x*z, -y*z])
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

@testset "== for DeterminantalGerm" begin
  R, (x,y,z) = QQ[:x,:y,:z]
  A = R[x y; 
        y z]
  X = DeterminantalGerm(A, 2, [0,0,0])
  # is identical
  X2 = X
  @test X == X2
  # different determinantal structure (with and without symmetry)
  X_sym = DeterminantalGerm(A, 2, [0,0,0], mat_type = :symmetric)
  @test X != X_sym
  # different ambient_coordinate_ring
  P, (a,b,c) = QQ[:a,:b,:c]
  B = P[a b;
        b c]
  Y = DeterminantalGerm(B, 2, [0,0,0])
  @test X != Y
  # same defining_matrix
  X3 = DeterminantalGerm(A, 2, [0,0,0])
  @test X == X3
  # different constructors for the same germ
  L, _ = localization(R, complement_of_point_ideal(R, [0,0,0]))
  X4 = DeterminantalGerm(L.(A), 2)
  @test X == X4
  # same determinantal_ideal, different defining matrix
  C = L[z y;
        y x]
  X5 = DeterminantalGerm(C, 2)
  @test determinantal_ideal(X4) == determinantal_ideal(X5)
  @test X4 != X5
end

@testset "T1_GL module" begin
  R, (x, y) = QQ[:x,:y]
  A = R[x  y^2;
       y^2 x^2]
  X_A = DeterminantalGerm(A, 2, [0,0])
  X_A_sym = DeterminantalGerm(A, 2, [0,0], mat_type = :symmetric)
  @test tjurina_GL_number(X_A) == 8
  @test tjurina_GL_number(X_A_sym) == 6
  B = R[x 0  0;
        0 y  x;
        0 x y^2]
  @test vector_space_dim(T1_GL_module(B, mat_type = :generic)[1]) == 11
  @test vector_space_dim(T1_GL_module(B, mat_type = :symmetric)[1]) == 7
  C = R[0     0    x    0;
        0     0    0 x^2+y^3;
       -x     0    0    0;
        0 -x^2-y^3 0    0]
  X_C = DeterminantalGerm(C, 4, [0,0])
  X_C_skew = DeterminantalGerm(C, 2, [0,0], mat_type = :skew_symmetric)
  M1, intr_1 = T1_GL_module(X_C_skew)
  L = localized_ring(OO(X_C_skew))
  M2, intr_2 = T1_GL_module(L.(C), mat_type = :skew_symmetric)
  @test intr_1.(repres.(vector_space_basis(M1))) == intr_2.(repres.(vector_space_basis(M2)))
  @test tjurina_GL_number(X_C) == PosInf()
  @test tjurina_GL_number(X_C_skew) == 16


  R, (v,w,x,y,z) = QQ[:v,:w,:x,:y,:z]
  A = R[x y z;
        v w 0]
  T1_GL_A = T1_GL_module(A)[1]
  @test krull_dim(T1_GL_A) == 2
  @test vector_space_dim(T1_GL_A) == PosInf()
  X_A = DeterminantalGerm(A, 2, [0,0,0,0,0])
  @test !is_EIDS(X_A)
  B = R[x y z;
        v w x]
  @test vector_space_dim(T1_GL_module(B)[1]) == 1
  @test_throws ArgumentError Oscar._T1_mat_module(B, Val{:wrong_symbol}, Val{:GL})
  @test_throws ArgumentError T1_GL_module(B, mat_type = :wrong_symbol)
end

@testset "T1_SL module" begin
  R, (x, y) = QQ[:x,:y]
  A = R[x^2  y;
         y  x^3]
  X_A = DeterminantalGerm(A, 2, [0,0])
  X_A_sym = DeterminantalGerm(A, 2, [0,0], mat_type = :symmetric)
  @test tjurina_SL_number(X_A) == 6
  @test tjurina_SL_number(X_A_sym) == 4
  B = R[x  0 y^2;
        0  y  x;
       y^2 x  0]
  X_B = DeterminantalGerm(B, 3, [0,0])
  X_B_sym = DeterminantalGerm(B, 3, [0,0], mat_type = :symmetric)
  @test vector_space_dim(T1_SL_module(B, mat_type = :generic)[1]) == 12
  @test vector_space_dim(T1_SL_module(B, mat_type = :symmetric)[1]) == 8
  # for ambient dimension 2 is tau_{GL}^{sym} equal to the milnor_number of a symmetric determinantal hypersurface germ
  @test milnor_number(HypersurfaceGerm(representative(X_B_sym), point(X_B_sym))) == tjurina_SL_number(X_B_sym)
  C = R[0    0   x  y^2;
        0    0  y^2 x^2;
       -x  -y^2  0   0;
      -y^2 -x^2  0   0]
  X_C = DeterminantalGerm(C, 4, [0,0])
  X_C_skew = DeterminantalGerm(C, 2, [0,0], mat_type = :skew_symmetric)
  M1, intr_1 = T1_SL_module(X_C_skew)
  L = localized_ring(OO(X_C_skew))
  M2, intr_2 = T1_SL_module(L.(C), mat_type = :skew_symmetric)
  @test intr_1.(repres.(vector_space_basis(M1))) == intr_2.(repres.(vector_space_basis(M2)))
  @test tjurina_SL_number(X_C) == PosInf()
  @test tjurina_SL_number(X_C_skew) == 12


  # for ambient dimension 3 is tau_{GL}^{sq} equal to the milnor_number of a generic determinantal hypersurface germ
  R, (x,y,z) = QQ[:x,:y,:z]
  A = R[  x   −y;
        y+z^7  x]
  X_A = DeterminantalGerm(A, 2, [0,0,0])
  @test tjurina_SL_number(X_A) == 13
  @test milnor_number(HypersurfaceGerm(representative(X_A), point(X_A))) == tjurina_SL_number(X_A)


  R, (v,w,x,y,z) = QQ[:v,:w,:x,:y,:z]
  A = R[x y z;
        v 0 0]
  T1_SL_A = T1_SL_module(A)[1]
  @test krull_dim(T1_SL_A) == 4
  @test vector_space_dim(T1_SL_A) == PosInf()
  B = R[x y z;
        v w x^2+y^2]
  @test vector_space_dim(T1_SL_module(B)[1]) == 3
  @test_throws ArgumentError Oscar._T1_mat_module(B, Val{:wrong_symbol}, Val{:SL})
  @test_throws ArgumentError T1_SL_module(B, mat_type=:wrong_symbol)
end
