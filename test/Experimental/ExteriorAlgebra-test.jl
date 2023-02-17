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


@test_throws  ArgumentError  exterior_algebra_naive(QQ, ["x", "y", "x"]); # duplicate name


## (reduced) COMPUTATIONAL SPEED TEST

ExtAlg, (e1,e2,e3,e4,e5,e6) = exterior_algebra_naive(QQ, 6);
fac1 = e1 + 2*e1*e2 + 3*e1*e2*e3 + 4*e1*e2*e3*e4 + 5*e1*e2*e3*e4*e5 + 6*e1*e2*e3*e4*e5*e6;
fac2 = e6 + 2*e5*e6 + 3*e4*e5*e6 + 4*e3*e4*e5*e6 + 5*e2*e3*e4*e5*e6 + 6*e1*e2*e3*e4*e5*e6;
prod12 = fac1*fac2;
prod21 = fac2*fac1;
expected12 = 35*e1*e2*e3*e4*e5*e6 +4*e1*e2*e3*e4*e6 +6*e1*e2*e3*e5*e6 +6*e1*e2*e4*e5*e6 +4*e1*e3*e4*e5*e6 +3*e1*e2*e3*e6 +4*e1*e2*e5*e6 +3*e1*e4*e5*e6 +2*e1*e2*e6 +2*e1*e5*e6 +e1*e6;
expected21 = - 3*e1*e2*e3*e4*e5*e6 + 4*e1*e2*e3*e4*e6 + 6*e1*e2*e3*e5*e6 + 6*e1*e2*e4*e5*e6 + 4*e1*e3*e4*e5*e6 - 3*e1*e2*e3*e6 + 4*e1*e2*e5*e6 - 3*e1*e4*e5*e6 + 2*e1*e2*e6 + 2*e1*e5*e6 - e1*e6;
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


@test_throws MethodError exterior_algebra_singular(ZZ, 3)  # Coeffs not field
@test_throws MethodError exterior_algebra_singular(ResidueRing(ZZ,4), 3)  # Coeffs not field
@test_throws  ArgumentError  exterior_algebra_singular(QQ, ["x", "y", "x"]); # duplicate name


## (reduced) COMPUTATIONAL SPEED TEST

ExtAlg, (e1,e2,e3,e4,e5,e6) = exterior_algebra_singular(QQ, 6);
fac1 = e1 + 2*e1*e2 + 3*e1*e2*e3 + 4*e1*e2*e3*e4 + 5*e1*e2*e3*e4*e5 + 6*e1*e2*e3*e4*e5*e6;
fac2 = e6 + 2*e5*e6 + 3*e4*e5*e6 + 4*e3*e4*e5*e6 + 5*e2*e3*e4*e5*e6 + 6*e1*e2*e3*e4*e5*e6;
prod12 = fac1*fac2;
prod21 = fac2*fac1;
expected12 = 35*e1*e2*e3*e4*e5*e6 +4*e1*e2*e3*e4*e6 +6*e1*e2*e3*e5*e6 +6*e1*e2*e4*e5*e6 +4*e1*e3*e4*e5*e6 +3*e1*e2*e3*e6 +4*e1*e2*e5*e6 +3*e1*e4*e5*e6 +2*e1*e2*e6 +2*e1*e5*e6 +e1*e6;
expected21 = - 3*e1*e2*e3*e4*e5*e6 + 4*e1*e2*e3*e4*e6 + 6*e1*e2*e3*e5*e6 + 6*e1*e2*e4*e5*e6 + 4*e1*e3*e4*e5*e6 - 3*e1*e2*e3*e6 + 4*e1*e2*e5*e6 - 3*e1*e4*e5*e6 + 2*e1*e2*e6 + 2*e1*e5*e6 - e1*e6;
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
