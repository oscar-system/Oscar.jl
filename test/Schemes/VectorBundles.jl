@testset "trivializations of coherent sheaves" begin
  IP2 = projective_space(QQ, 2)
  S = ambient_coordinate_ring(IP2)
  (x,y,z) = gens(S)
  f = x^2 + y*z
  IPC = subscheme(IP2, f)
  C = covered_scheme(IPC)
  set_name!(C, "C")
  WC = cotangent_sheaf(C)
  @test is_locally_free(WC)
  V = oscar.trivializing_covering(WC)
  U = affine_charts(C)
  @test WC(U[1], V[1]) isa ModuleFPHom
  W = PrincipalOpenSubset(V[2], first(gens(OO(V[2]))))
  res = WC(V[1], W)
  v = res(first(gens(WC(V[1]))))
  @test coordinates(v)[1] == -inv(first(gens(OO(W))))^2

  TC = tangent_sheaf(C)
  @test_broken oscar.trivializing_covering(TC) isa Covering

  IP4 = projective_space(QQ, 4)
  S = ambient_coordinate_ring(IP4)
  (x,y,z,u,v) = gens(S)
  A = S[x y z; u v x]
  I = ideal(S, minors(A, 2))
  IPX = subscheme(IP4, I)
  X = covered_scheme(IPX)
  set_name!(X, "X")
  WX = cotangent_sheaf(X)
  C = oscar.trivializing_covering(WX)
  @test WX(C[1]) isa FreeMod
  @test !any(x->x===ambient_scheme(C[1]), affine_charts(X)) # codimension 2 means recursion depth >= 2.
end
