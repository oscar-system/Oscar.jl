@testset "tor for modules over localizations at prime ideals" begin
  R, (x0, x1, x2, x3) = polynomial_ring(QQ, [:x0, :x1, :x2, :x3])
  p = ideal(R, [x0-x1,x2,x3])
  Rp, _ = localization(R, complement_of_prime_ideal(p))
  I = ideal(Rp, [x0-x1,x2-x3])
  Ym = ideal(Rp, [x0^3+(x2^2-x3^2)*x0-x1^3+x2^3+x3^3])
  Rp1 = FreeMod(Rp, 1)
  IRp1, inc_I = I*Rp1
  YmRp1, inc_Ym = Ym*Rp1
  M1 = cokernel(inc_I)
  M2 = cokernel(inc_Ym)
  H0, _ = tor(M1, M2, 0)
  @test length(H0) == 3

  R, (x_0, x_1, x_2, x_3) = polynomial_ring(QQ, [:x0, :x1, :x2, :x3]);
  p = ideal(R, [x_0-x_1,x_2,x_3]);
  Rp, _ = localization(R, complement_of_prime_ideal(p););
  I = ideal(Rp, [x_0-x_1,x_2-x_3]);
  Ym = ideal(Rp, [x_0^3+(x_2^2-x_3^2)*x_0-x_1^3+x_2^3+x_3^3]);
  Rp1 = FreeMod(Rp, 1);
  IRp1, inc_I = I*Rp1;
  YmRp1, inc_Ym = Ym*Rp1;
  M1 = cokernel(inc_I);
  M2 = cokernel(inc_Ym);
  H0, _ = tor(M1, M2, 0);
  @test length(H0) == 3

  i = ideal(Rp, [x_0-x_1,x_2^2+x_3^2]);
  Yp = ideal(Rp, [x_0^3+(x_2^2+x_3^2)*x_0-x_1^3+x_2^3+x_3^3]);
  iRp1, inc_i = i*Rp1;
  YpRp1, inc_Yp = Yp*Rp1;
  N1 = cokernel(inc_i);
  N2 = cokernel(inc_Yp);
  H1, _ = tor(N1, N2, 0);
  @test length(H1) == 6

  q = ideal(R, [x_0^2+x_0*x_1+x_1^2,x_2,x_3]);
  Rq, _ = localization(R, complement_of_prime_ideal(q););
  II = ideal(Rq, [x_0^2+x_0*x_1+x_1^2,x_2-x_3]);
  YYm = ideal(Rq, [x_0^3+(x_2^2-x_3^2)*x_0-x_1^3+x_2^3+x_3^3]);
  Rq1 = FreeMod(Rq, 1);
  IRq1, inc_II = II*Rq1;
  YYmRq1, inc_YYm = YYm*Rq1;
  MM1 = cokernel(inc_II);
  MM2 = cokernel(inc_YYm);
  HH0, _ = tor(MM1, MM2, 0);
  @test length(HH0) == 3

  ii = ideal(Rq, [x_0^2+x_0*x_1+x_1^2,x_2^2+x_3^2]);
  YYp = ideal(Rq, [x_0^3+(x_2^2+x_3^2)*x_0-x_1^3+x_2^3+x_3^3]);
  iiRq1, inc_ii = ii*Rq1;
  YYpRq1, inc_YYp = YYp*Rq1;
  NN1 = cokernel(inc_ii);
  NN2 = cokernel(inc_YYp);
  HH1, _ = tor(NN1, NN2, 0);
  @test length(HH1) == 6
end

