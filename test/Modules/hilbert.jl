@testset "Hilbert series and free resolution" begin
  Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]; weights = [3,1,1]);
  F = graded_free_module(Rg, 1);
  A = Rg[x; y];
  B = Rg[x^2+y^6; y^7; z^4];
  M = SubquoModule(F, A, B);
  fr = free_resolution(M)

  # 2023-08-18 There is no length(fr), so instead we do the following
  indexes = range(fr.C);
  fr_len = first(indexes)

  num,_ = hilbert_series(fr[fr_len])
  for i in fr_len:-1:1
    phi = map(fr,i)
    N = cokernel(phi)
    numer,denom = hilbert_series(N)
    numer_next,_ = hilbert_series(fr[i-1])
    @test numer ==  numer_next- num
    num = numer
  end

  num,_ = hilbert_series(fr[fr_len]; parent=laurent_polynomial_ring(QQ, :T; cached=false)[1])
  for i in fr_len:-1:1
    phi = map(fr,i)
    N = cokernel(phi)
    numer,denom = hilbert_series(N; parent=parent(num))
    numer_next,_ = hilbert_series(fr[i-1]; parent=parent(num))
    dummy_num, _ = hilbert_series(N)
    @test_throws ErrorException dummy_num - num
    @test numer ==  numer_next- num
    num = numer
  end
end

