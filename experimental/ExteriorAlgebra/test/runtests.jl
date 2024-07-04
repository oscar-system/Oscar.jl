@testset "ExteriorAlgebra" begin

  ### 2023-02-10  Tests for experimental/ExteriorAlgebra/ExteriorAlgebra.jl

  ######## CONSTRUCTOR TESTs
  @test_throws ArgumentError exterior_algebra(QQ, 0)
  @testset "exterior_algebra constructor" for (R,n) in [
      (QQ, 1)   # --> special case (is commutative)
      (QQ, 2)   # --> first general case
      (QQ, 99)  # --> Also tried with 999 indets, but takes a long time [~19s]

      (GF(2), 2)  # BUG??  not recognized as commutative!!
      (GF(3), 4)
      (residue_field(ZZ, 2)[1], 2)
      (residue_field(ZZ, 3)[1], 4)
      #(GF(1180591620717411303449), 2)
      #(residue_field(ZZ, 1180591620717411303449)[1], 2)
      ## (GF(2), 1500);   ## limit 1500 on my 32Gbyte machine (!NB printing requires a lot of space!)
  ]
    A, g = exterior_algebra(R, n)
    @test A isa PBWAlgQuo
    @test g isa Vector{elem_type(A)}
    @test length(g) == n
    @test ngens(A) == n
    @test gens(A) == g
  end

  # Duplicate names are allowed
  exterior_algebra(QQ, ["x", "y", "x"])

  @test_throws MethodError exterior_algebra(ZZ, 3)                  # Coeffs not field
  @test_throws MethodError exterior_algebra(residue_ring(ZZ, 4)[1], 3)  # Coeffs not field

  @test_throws ArgumentError exterior_algebra(QQ, String[])        # empty name list

  ## (reduced) COMPUTATIONAL SPEED TEST

  ExtAlg, (e1, e2, e3, e4, e5, e6) = exterior_algebra(QQ, 6)
  fac1 =
    e1 +
    2 * e1 * e2 +
    3 * e1 * e2 * e3 +
    4 * e1 * e2 * e3 * e4 +
    5 * e1 * e2 * e3 * e4 * e5 +
    6 * e1 * e2 * e3 * e4 * e5 * e6
  fac2 =
    e6 +
    2 * e5 * e6 +
    3 * e4 * e5 * e6 +
    4 * e3 * e4 * e5 * e6 +
    5 * e2 * e3 * e4 * e5 * e6 +
    6 * e1 * e2 * e3 * e4 * e5 * e6
  prod12 = fac1 * fac2
  prod21 = fac2 * fac1
  expected12 =
    35 * e1 * e2 * e3 * e4 * e5 * e6 +
    4 * e1 * e2 * e3 * e4 * e6 +
    6 * e1 * e2 * e3 * e5 * e6 +
    6 * e1 * e2 * e4 * e5 * e6 +
    4 * e1 * e3 * e4 * e5 * e6 +
    3 * e1 * e2 * e3 * e6 +
    4 * e1 * e2 * e5 * e6 +
    3 * e1 * e4 * e5 * e6 +
    2 * e1 * e2 * e6 +
    2 * e1 * e5 * e6 +
    e1 * e6
  expected21 =
    -3 * e1 * e2 * e3 * e4 * e5 * e6 +
    4 * e1 * e2 * e3 * e4 * e6 +
    6 * e1 * e2 * e3 * e5 * e6 +
    6 * e1 * e2 * e4 * e5 * e6 +
    4 * e1 * e3 * e4 * e5 * e6 - 3 * e1 * e2 * e3 * e6 + 4 * e1 * e2 * e5 * e6 -
    3 * e1 * e4 * e5 * e6 +
    2 * e1 * e2 * e6 +
    2 * e1 * e5 * e6 - e1 * e6
  @test is_zero(prod12 - expected12)
  @test is_zero(prod21 - expected21)

  @test is_zero(fac1 * prod12)
  @test is_zero(prod12 * fac1)
  @test is_zero(fac2 * prod12)
  @test is_zero(prod12 * fac2)

  @test is_zero(fac1 * prod21)
  @test is_zero(prod21 * fac1)
  @test is_zero(fac2 * prod21)
  @test is_zero(prod21 * fac2)
end

# # -------------------------------------------------------
# # Duplicate of test above, but using the PBWAlgQuo implementation:

# @testset "ExteriorAlgebra_PBWAlgQuo" begin

# ### 2023-02-10  Tests for experimental/ExteriorAlgebra/ExteriorAlgebra_PBWAlgQuo.jl

# ######## CONSTRUCTOR TESTs
# @test_throws  ArgumentError  exterior_algebra_PBWAlgQuo(QQ, 0);
# exterior_algebra_PBWAlgQuo(QQ, 1);  # --> special case (is commutative)
# exterior_algebra_PBWAlgQuo(QQ, 2);  # --> first general case
# exterior_algebra_PBWAlgQuo(QQ, 99);  # -->  with 999 indets takes a long time [~19s]

# exterior_algebra_PBWAlgQuo(GF(2), 2);  # BUG??  not recognized as commutative!!
# exterior_algebra_PBWAlgQuo(GF(3), 4);
# ##  exterior_algebra_PBWAlgQuo(GF(2), 1500);   ## limit 1500 on my 32Gbyte machine (!NB printing requires a lot of space!)

# ### exterior_algebra_PBWAlgQuo(GF(1180591620717411303449), 2);  #  --> ERROR prime too big (for GF)

# exterior_algebra_PBWAlgQuo(residue_field(ZZ, 2)[1], 2);
# exterior_algebra_PBWAlgQuo(residue_field(ZZ, 3)[1], 4);
# exterior_algebra_PBWAlgQuo(residue_field(ZZ, 1180591620717411303449)[1], 2);

# exterior_algebra_PBWAlgQuo(residue_ring(ZZ,4)[1], 3)  # Coeffs not integral domain

# @test_throws  ArgumentError  exterior_algebra_PBWAlgQuo(QQ, String[]); # empty name list
# exterior_algebra_PBWAlgQuo(QQ, ["x", "y", "x"]); #  !!duplicate names are allowed!!

# ## (reduced) COMPUTATIONAL SPEED TEST

# ExtAlg, (e1,e2,e3,e4,e5,e6) = exterior_algebra_PBWAlgQuo(QQ, 6);
# fac1 = e1 + 2*e1*e2 + 3*e1*e2*e3 + 4*e1*e2*e3*e4 + 5*e1*e2*e3*e4*e5 + 6*e1*e2*e3*e4*e5*e6;
# fac2 = e6 + 2*e5*e6 + 3*e4*e5*e6 + 4*e3*e4*e5*e6 + 5*e2*e3*e4*e5*e6 + 6*e1*e2*e3*e4*e5*e6;
# prod12 = fac1*fac2;
# prod21 = fac2*fac1;
# expected12 = 35*e1*e2*e3*e4*e5*e6 +4*e1*e2*e3*e4*e6 +6*e1*e2*e3*e5*e6 +6*e1*e2*e4*e5*e6 +4*e1*e3*e4*e5*e6 +3*e1*e2*e3*e6 +4*e1*e2*e5*e6 +3*e1*e4*e5*e6 +2*e1*e2*e6 +2*e1*e5*e6 +e1*e6;
# expected21 = - 3*e1*e2*e3*e4*e5*e6 + 4*e1*e2*e3*e4*e6 + 6*e1*e2*e3*e5*e6 + 6*e1*e2*e4*e5*e6 + 4*e1*e3*e4*e5*e6 - 3*e1*e2*e3*e6 + 4*e1*e2*e5*e6 - 3*e1*e4*e5*e6 + 2*e1*e2*e6 + 2*e1*e5*e6 - e1*e6;
# @test is_zero(prod12 - expected12);
# @test is_zero(prod21 - expected21);

# @test is_zero(fac1*prod12);
# @test is_zero(prod12*fac1);
# @test is_zero(fac2*prod12);
# @test is_zero(prod12*fac2);

# @test is_zero(fac1*prod21);
# @test is_zero(prod21*fac1);
# @test is_zero(fac2*prod21);
# @test is_zero(prod21*fac2);

# end
