@testset "ExteriorAlgebra_naive" begin

### 2023-02-10  Tests for experimental/ExteriorAlgebra/ExteriorAlgebra_naive.jl

######## CONSTRUCTOR TESTs
@test_throws  ArgumentError  exterior_algebra_naive(QQ, 0);
exterior_algebra_naive(QQ, 1);  # --> special case (is commutative)
exterior_algebra_naive(QQ, 2);  # --> first general case
exterior_algebra_naive(QQ, 99);  # -->  with 999 indets takes a long time [~19s]


exterior_algebra_naive(GF(2), 2);  # BUG??  not recognized as commutative!!
exterior_algebra_naive(GF(3), 4);
##  exterior_algebra_naive(GF(2), 1500);   ## limit 1500 on my 32Gbyte machine (!NB printing requires a lot of space!)

### exterior_algebra_naive(GF(1180591620717411303449), 2);  #  --> ERROR prime too big (for GF)


exterior_algebra_naive(ResidueField(ZZ,2), 2);
exterior_algebra_naive(ResidueField(ZZ,3), 4);
exterior_algebra_naive(ResidueField(ZZ,1180591620717411303449), 2);

exterior_algebra_naive(ResidueRing(ZZ,4), 3)  # Coeffs not integral domain


@test_throws  ArgumentError  exterior_algebra_naive(QQ, 2; var_prefix="@");
@test_throws  ArgumentError  exterior_algebra_naive(QQ, ["x", "y", "x"]);


## (reduced) COMPUTATIONAL SPEED TEST

ExtAlg, (x1,x2,x3,x4,x5,x6) = exterior_algebra_naive(QQ, 6; var_prefix="x");
fac1 = x1 + 2*x1*x2 + 3*x1*x2*x3 + 4*x1*x2*x3*x4 + 5*x1*x2*x3*x4*x5 + 6*x1*x2*x3*x4*x5*x6;
fac2 = x6 + 2*x5*x6 + 3*x4*x5*x6 + 4*x3*x4*x5*x6 + 5*x2*x3*x4*x5*x6 + 6*x1*x2*x3*x4*x5*x6;
prod12 = fac1*fac2;
prod21 = fac2*fac1;
expected12 = 35*x1*x2*x3*x4*x5*x6 +4*x1*x2*x3*x4*x6 +6*x1*x2*x3*x5*x6 +6*x1*x2*x4*x5*x6 +4*x1*x3*x4*x5*x6 +3*x1*x2*x3*x6 +4*x1*x2*x5*x6 +3*x1*x4*x5*x6 +2*x1*x2*x6 +2*x1*x5*x6 +x1*x6;
expected21 = - 3*x1*x2*x3*x4*x5*x6 + 4*x1*x2*x3*x4*x6 + 6*x1*x2*x3*x5*x6 + 6*x1*x2*x4*x5*x6 + 4*x1*x3*x4*x5*x6 - 3*x1*x2*x3*x6 + 4*x1*x2*x5*x6 - 3*x1*x4*x5*x6 + 2*x1*x2*x6 + 2*x1*x5*x6 - x1*x6;
@test is_zero(prod12 - expected12);
@test is_zero(prod21 - expected21);

@test is_zero(fac1*prod12);
@test is_zero(prod12*fac1);
@test is_zero(fac2*prod12);
@test is_zero(prod12*fac2);

@test is_zero(fac1*prod21);
@test is_zero(prod21*fac1);
@test is_zero(fac2*prod21);
@test is_zero(prod21*fac2);

end



# Duplicate of test above, but using the singular implementation:

@testset "ExteriorAlgebra_singular" begin

### 2023-02-10  Tests for experimental/ExteriorAlgebra/ExteriorAlgebra_singular.jl

######## CONSTRUCTOR TESTs
@test_throws  ArgumentError  exterior_algebra_singular(QQ, 0);
exterior_algebra_singular(QQ, 1);  # --> special case (is commutative)
exterior_algebra_singular(QQ, 2);  # --> first general case
exterior_algebra_singular(QQ, 99);  # -->  with 999 indets takes a long time [~19s]


exterior_algebra_singular(GF(2), 2);  # BUG??  not recognized as commutative!!
exterior_algebra_singular(GF(3), 4);
##  exterior_algebra_singular(GF(2), 1500);   ## limit 1500 on my 32Gbyte machine (!NB printing requires a lot of space!)

### exterior_algebra_singular(GF(1180591620717411303449), 2);  #  --> ERROR prime too big (for GF)


exterior_algebra_singular(ResidueField(ZZ,2), 2);
exterior_algebra_singular(ResidueField(ZZ,3), 4);
exterior_algebra_singular(ResidueField(ZZ,1180591620717411303449), 2);

exterior_algebra_singular(ResidueRing(ZZ,4), 3)  # Coeffs not integral domain


@test_throws  ArgumentError  exterior_algebra_singular(QQ, 2; var_prefix="@");
@test_throws  ArgumentError  exterior_algebra_singular(QQ, ["x", "y", "x"]);


## (reduced) COMPUTATIONAL SPEED TEST

ExtAlg, (x1,x2,x3,x4,x5,x6) = exterior_algebra_singular(QQ, 6; var_prefix="x");
fac1 = x1 + 2*x1*x2 + 3*x1*x2*x3 + 4*x1*x2*x3*x4 + 5*x1*x2*x3*x4*x5 + 6*x1*x2*x3*x4*x5*x6;
fac2 = x6 + 2*x5*x6 + 3*x4*x5*x6 + 4*x3*x4*x5*x6 + 5*x2*x3*x4*x5*x6 + 6*x1*x2*x3*x4*x5*x6;
prod12 = fac1*fac2;
prod21 = fac2*fac1;
expected12 = 35*x1*x2*x3*x4*x5*x6 +4*x1*x2*x3*x4*x6 +6*x1*x2*x3*x5*x6 +6*x1*x2*x4*x5*x6 +4*x1*x3*x4*x5*x6 +3*x1*x2*x3*x6 +4*x1*x2*x5*x6 +3*x1*x4*x5*x6 +2*x1*x2*x6 +2*x1*x5*x6 +x1*x6;
expected21 = - 3*x1*x2*x3*x4*x5*x6 + 4*x1*x2*x3*x4*x6 + 6*x1*x2*x3*x5*x6 + 6*x1*x2*x4*x5*x6 + 4*x1*x3*x4*x5*x6 - 3*x1*x2*x3*x6 + 4*x1*x2*x5*x6 - 3*x1*x4*x5*x6 + 2*x1*x2*x6 + 2*x1*x5*x6 - x1*x6;
@test is_zero(prod12 - expected12);
@test is_zero(prod21 - expected21);

@test is_zero(fac1*prod12);
@test is_zero(prod12*fac1);
@test is_zero(fac2*prod12);
@test is_zero(prod12*fac2);

@test is_zero(fac1*prod21);
@test is_zero(prod21*fac1);
@test is_zero(fac2*prod21);
@test is_zero(prod21*fac2);

end
