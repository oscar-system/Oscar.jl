@testset "Hilbert series" begin
  P, x = graded_polynomial_ring(QQ, 5, :x, [1 1 1 1 1; 1 -2 3 -4 5]);
  G =
  [
   x[1]^30*x[2]^16*x[3]^7*x[4]^72*x[5]^31,
   x[1]^11*x[2]^20*x[3]^29*x[4]^93*x[5]^5,
   x[1]^5*x[2]^9*x[3]^43*x[4]^7*x[5]^97,
   x[1]^57*x[2]^6*x[3]^7*x[4]^56*x[5]^37,
   x[1]^22*x[2]^40*x[3]^27*x[4]^67*x[5]^18,
   x[1]^5*x[2]^5*x[3]^7*x[4]^78*x[5]^81,
   x[1]^85*x[2]^15*x[3]^12*x[4]^41*x[5]^37,
   x[1]^13*x[2]^35*x[3]^2*x[4]^64*x[5]^84,
   x[1]^31*x[2]^57*x[3]^12*x[4]^93*x[5]^6,
   x[1]^37*x[2]^59*x[3]^88*x[4]^16*x[5]^7,
   x[1]^47*x[2]^63*x[3]^32*x[4]^28*x[5]^41,
   x[1]^18*x[2]^19*x[3]^95*x[4]^7*x[5]^73,
   x[1]^37*x[2]^84*x[3]^11*x[4]^48*x[5]^32,
   x[1]^85*x[2]*x[3]^73*x[4]^24*x[5]^31,
   x[1]^59*x[2]^56*x[3]^41*x[4]^12*x[5]^50,
   x[1]^58*x[2]^36*x[3]*x[4]^35*x[5]^89,
   x[1]^71*x[2]^12*x[3]^36*x[4]^76*x[5]^25,
   x[1]^58*x[2]^32*x[3]^85*x[4]^44*x[5],
   x[1]^61*x[2]^81*x[3]^15*x[4]^59*x[5]^8,
   x[1]^84*x[2]^49*x[3]^32*x[4]^52*x[5]^7,
   x[1]^40*x[2]^3*x[3]^54*x[4]^39*x[5]^89,
   x[1]^12*x[2]^52*x[3]^100*x[4]^49*x[5]^14,
   x[1]^75*x[2]^78*x[3]^34*x[4]*x[5]^41,
   x[1]^46*x[2]^22*x[3]^99*x[4]^49*x[5]^14,
   x[1]^95*x[2]^11*x[3]^63*x[4]^52*x[5]^10,
   x[1]^32*x[2]^47*x[3]^97*x[4]^32*x[5]^26,
   x[1]^52*x[2]^64*x[3]^62*x[4]^13*x[5]^47,
   x[1]^70*x[2]^46*x[3]^90*x[4]^13*x[5]^21,
   x[1]^3*x[2]^67*x[3]^90*x[4]^45*x[5]^52,
   x[1]^17*x[2]^100*x[3]^58*x[4]^62*x[5]^21,
   x[1]^45*x[2]^54*x[3]^65*x[4]^64*x[5]^32,
   x[1]^24*x[2]^85*x[3]^27*x[4]^49*x[5]^76,
   x[1]^28*x[2]^83*x[3]^9*x[4]^97*x[5]^45,
   x[1]^40*x[2]^52*x[3]^99*x[4]^27*x[5]^49,
   x[1]^88*x[2]^36*x[3]^26*x[4]^30*x[5]^90,
   x[1]^36*x[2]^68*x[3]^76*x[4]^81*x[5]^9,
   x[1]^48*x[2]^25*x[3]^80*x[4]^40*x[5]^83,
   x[1]^93*x[2]^14*x[3]^80*x[4]^89*x[5]^3,
   x[1]^5*x[2]^61*x[3]^65*x[4]^65*x[5]^94,
   x[1]^48*x[2]^91*x[3]^70*x[4]^66*x[5]^16,
   x[1]^39*x[2]^79*x[3]^98*x[4]^2*x[5]^75,
   x[1]^4*x[2]^60*x[3]^74*x[4]^56*x[5]^100,
   x[1]^59*x[2]^100*x[3]^74*x[4]^51*x[5]^12,
   x[1]^90*x[2]^61*x[3]^85*x[4]^43*x[5]^19,
   x[1]^44*x[2]^97*x[3]^39*x[4]^27*x[5]^97,
   x[1]^91*x[2]^37*x[3]^2*x[4]^97*x[5]^80,
   x[1]^67*x[2]^44*x[3]^99*x[4]^5*x[5]^93,
   x[1]^95*x[2]^93*x[3]^88*x[4]^30*x[5]^2,
   x[1]^91*x[3]^84*x[4]^59*x[5]^76,
   x[1]^62*x[2]^3*x[3]^96*x[4]^72*x[5]^84
  ];

  I = ideal(P,G);
  PmodI, _ = quo(P,I);
  (num1,_), (_,_) = multi_hilbert_series(PmodI);
  S = parent(num1)
  (num2,_), (_,_) = multi_hilbert_series(PmodI; parent=S, backend=:Zach);
  @test num2 == num1

  P2, x = graded_polynomial_ring(QQ, 5, :x, [1 1 1 1 1])
  G = P2.(G);

  I = ideal(P2,G);
  PmodI, _ = quo(P2,I);
  num1, _ = hilbert_series(PmodI);
  L, q = laurent_polynomial_ring(ZZ, :q);
  num2, _ = hilbert_series(PmodI; parent=L);
  @test num2 == evaluate(num1, q)
end

@testset "second round of Hilbert series" begin
  # Standard graded
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);
  J = ideal(R, [y-x^2, z-x^3]);
  H = homogenizer(R, "w")
  I = H(J)
  A, _ = quo(base_ring(I), I);
  numer1, denom1 = hilbert_series(A);
  S, t = laurent_polynomial_ring(ZZ, :t);
  numer2, denom2 = hilbert_series(A; parent=S);
  Smult, (T,) = polynomial_ring(ZZ, [:t]);
  numer3, denom3 = hilbert_series(A; parent=Smult);
  @test numer3 == evaluate(numer1, T)
  @test numer2 == evaluate(numer1, t)

  # Negative grading
  RR, (X, Y) = graded_polynomial_ring(QQ, [:X, :Y], [-1, -1])
  JJ = ideal(RR, X^2 - Y^2);
  A, _ = quo(base_ring(JJ), JJ);
  (numer1, denom1), _ = multi_hilbert_series(A);
  S, t = laurent_polynomial_ring(ZZ, :t);
  (numer2, denom2), _ = multi_hilbert_series(A; parent=S);
  Smult, (T,) = polynomial_ring(ZZ, [:t]);
  @test_throws  DomainError  multi_hilbert_series(A; parent=Smult)

  # Graded by commutative group
  G = free_abelian_group(2);
  G, _ = quo(G, [G[1]-3*G[2]]);
  RR, (X, Y) = graded_polynomial_ring(QQ, [:X, :Y], [G[1], G[2]]);
  JJ = ideal(RR, X^2 - Y^6);
  A, _ = quo(base_ring(JJ), JJ);
  (numer1, denom1), (H, iso) = multi_hilbert_series(A);
  @test is_free(H) && isone(torsion_free_rank(H))
  S, t = laurent_polynomial_ring(ZZ, :t);
  (numer2, denom2), (H, iso) = multi_hilbert_series(A; parent=S);
  @test is_free(H) && isone(torsion_free_rank(H))
  Smult, (T,) = polynomial_ring(ZZ, [:t]);
  (numer3, denom3), (H, iso) = multi_hilbert_series(A; parent=Smult);
  @test is_free(H) && isone(torsion_free_rank(H))
  @test numer1 == evaluate(numer2, T)
  @test numer3 == evaluate(numer2, first(gens(parent(numer3))))
end

@testset "Hilbert series part 3" begin
  R, (x, y) = graded_polynomial_ring(QQ, [:x, :y], [1, -1]);
  I = ideal(R, [x]);
  RmodI, _ = quo(R, I);
  @test_throws ArgumentError  hilbert_series(RmodI);  # weights must be non-neg
  @test_throws ArgumentError  multi_hilbert_series(RmodI); # possible infinite dimension
end

@testset "Hilbert series part 4" begin
  # This test verifies that an intermediate overflow that was reported in
  # https://github.com/oscar-system/Oscar.jl/issues/2411 is fixed.
  A, x = graded_polynomial_ring(QQ, :x => 1:37)
  I = ideal([2*x[11] - 2*x[17] - 2*x[24] + 2*x[32] - 111916*x[37], 2*x[4] - 2*x[8] - 2*x[26] + 2*x[34] - 41216*x[37], 2*x[2] - 2*x[9] - 2*x[20] + 2*x[35] + 37974*x[37], x[28] - x[36], x[21] - x[36], x[27] - x[28] + x[33] + x[36], x[26] - x[27] - x[33] + x[34], x[20] - x[21] + x[35] + x[36], x[15] - x[21] - x[28] + x[36], x[10] - x[36], x[25] - x[28] + x[31] + x[36], x[24] - x[25] - x[26] + x[27] - x[31] + x[32] + x[33] - x[34], -x[14] + x[15] + x[18] - x[21] + x[25] - x[28] + x[31] + x[36], x[13] - x[14] + x[18] - x[19] - 2*x[20] + 2*x[21] - x[26] + x[27] + x[33] - x[34] - 2*x[35] - 2*x[36], x[9] - x[10] + x[35] + x[36], x[6] - x[10] - x[28] + x[36], x[19] - x[21] + x[30] + x[36], -x[18] + x[19] + x[23] - x[25] - x[27] + x[28] + x[30] - x[31] - x[33] - x[36], x[17] - x[19] - x[30] + x[32], x[12] - x[14] - x[17] + x[18] - x[27] + x[28] + x[31] - x[32] - x[33] - x[36], x[8] - x[10] + x[34] + x[36], x[5] - x[6] - x[8] + x[10] - x[27] + x[28] - x[34] - x[36], x[3] - x[10] - x[21] + x[36], -x[18] + x[19] + x[20] - x[21] + x[29] + x[30] + x[35] + x[36], x[22] + x[23] + x[24] - x[25] - x[29] - x[30] - x[31] + x[32], x[16] + x[17] + x[18] - x[19] - x[22] - x[23] - x[24] + x[25], x[11] + x[12] + x[13] - x[14] - x[16] - x[17] - x[18] + x[19] + x[22] + x[23] + x[24] - x[25] + x[29] + x[30] + x[31] - x[32], x[7] + x[8] + x[9] - x[10] - x[33] + x[34] + x[35] + x[36], x[4] + x[5] + x[9] - x[10] + x[26] - x[27] + x[35] + x[36], x[2] + x[3] + x[9] - x[10] + x[20] - x[21] + x[35] + x[36], x[1] - x[3] - x[6] + x[10] - x[15] + x[21] + x[28] - x[36], -x[27]*x[36] + x[34]*x[35], -x[25]*x[36] + x[32]*x[35], x[14]*x[36] + x[19]*x[35] + x[25]*x[36] + x[27]*x[36] - x[32]*x[35] - x[34]*x[35], -x[19]*x[36] - x[25]*x[36] + x[32]*x[34] + x[32]*x[35], -x[19]*x[35] - x[19]*x[36] + x[25]*x[34] - x[25]*x[36] + x[32]*x[34] + x[32]*x[35], x[14]*x[36] - x[19]*x[35] + x[25]*x[34] + x[27]*x[32], x[14]*x[35] - x[14]*x[36] + x[19]*x[35] - x[19]*x[36] + x[25]*x[27] - x[25]*x[34] - x[27]*x[32] + x[32]*x[34], x[14]*x[34] + x[19]*x[27] - 2*x[19]*x[35] + 2*x[25]*x[34] - x[25]*x[36] + x[32]*x[35], x[14]*x[32] - 2*x[14]*x[36] + x[19]*x[25] - 2*x[19]*x[35] - x[27]*x[36] + x[34]*x[35]])
  Q, _ = quo(A, I)
  h = hilbert_polynomial(Q)
  @test coefficients(h)[1] == 6
  @test coefficients(h)[2] == 0
  @test degree(Q) == 6
end

@testset "Hilbert series part 5"  begin
  # Test a trivial case -- just to check that the 1 power product prints correctly
  P,(x,y) = graded_polynomial_ring(QQ, ["x","y"]);
  @test  repr(Oscar.PP(degrees(one(P)))) == "1"
  I = ideal(one(P));
  PmodI,_ = quo(P,I); # Q is the zero ring
  HSRing1,_ = polynomial_ring(ZZ, "t");
  @test is_zero(Oscar.HSNum_abbott(PmodI, HSRing1))
end
