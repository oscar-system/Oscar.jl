@testset "mpoly-localizations" begin
  R, var = ZZ["x", "y"]
  x = var[1]
  y = var[2] 
  f = x^2 + y^2 -1
  m = ideal(R, [x, y])
  I = ideal(R, f)
  S = MPolyComplementOfPrimeIdeal(I)
  V = localize_at(S)
  T = MPolyComplementOfKPointIdeal(R, [ZZ(1), ZZ(0)])
  W = localize_at(T)
  
  k = QQ
  R, var = k["x", "y"]
  x = var[1]
  y = var[2] 
  p = 123
  q = -49
  f = x^2 + y^2 - p^2 - q^2
  m = ideal(R, [x, y])
  I = ideal(R, f)
  S = MPolyComplementOfPrimeIdeal(I)
  @test ambient_ring(S) == R
  @test !( f in S )
  @test x in S
  V = localize_at(S)
  @test original_ring(V) == R
  @test inverted_set(V) == S

  T = MPolyComplementOfKPointIdeal(R, [k(p), k(q)])
  @test ambient_ring(T) == R
  @test typeof(point_coordinates(T)) == Vector{elem_type(k)}
  @test x in T
  @test !(x-p in T)
  W = localize_at(T)

  a = W(R(1))
  b = W(2)
  c = W(y//x)
  @test a*(b+c) == a*b + a*c
  @test a*(b+c)^2 != a*(b^2) + a*c^2
  b = W((x-19)//(y-2))
  W(1)//b
  c = W((f+3)//(y-9))
  @test c//b == 1//(b//c)

  I_loc = ideal(W, gens(I))
  @test gens(I_loc) == [W(f)]
  I_loc = ideal(W, gens(I)[1])
  @test gens(I_loc) == [W(f)]
  #@test W(f) in I_loc
  
  J = ideal(W, [W(x-3), W(y-4)])
  lbpa = groebner_basis(J)
  @test ordering(lbpa) == :negdegrevlex
  @test oscar_gens(lbpa)[1] == W(1)

  J = ideal(W, [f, y-q])
  lbpa = groebner_basis(J)
  @test ordering(lbpa) == :negdegrevlex
  @test oscar_gens(lbpa) == W.([x-p, y-q])

  @test I_loc + J isa MPolyLocalizedIdeal
  @test I_loc * J isa MPolyLocalizedIdeal
  @test J^5 isa MPolyLocalizedIdeal

  @test reduce(W(x)//W(y-q+1), lbpa) == W(p)//W(y-q+1)

  K = ideal(W, f)
  @test f*(x-p+4) in K
  @test !(f+2 in K)
  @test f*(x-p+9) in I_loc
end
