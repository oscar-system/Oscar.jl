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
  @test oscar.trivializing_covering(TC) isa Covering

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
  TX = tangent_sheaf(X)
  CC = oscar.trivializing_covering(TX)

  # Testing transitions across charts while going down in the trees.
  @test all(x->TX(x) isa FreeMod, patches(CC))
  A = PrincipalOpenSubset(CC[1], gens(OO(CC[1]))[1])
  UU = simplify(A)
  @test TX(CC[1], UU) isa ModuleFPHom
  UUU = PrincipalOpenSubset(UU, one(OO(UU)))
  @test TX(CC[1], UUU) isa ModuleFPHom
  B = PrincipalOpenSubset(CC[2], one(OO(CC[2])))
  @test oscar.is_open_func(TX)(A, B)
  @test WX(B, A) isa ModuleFPHom
  @test TX(B, A) isa ModuleFPHom
  @test TX(B, UU) == compose(TX(B, A), TX(A, UU))
  @test TX(B, UUU) == compose(TX(B, A), TX(A, UUU))
  @test TX(B, UUU) == compose(TX(B, UU), TX(UU, UUU))
  BU = simplify(B)
  BUU = PrincipalOpenSubset(BU, one(OO(BU)))
  @test WX(BU, A) isa ModuleFPHom
  @test TX(BU, A) isa ModuleFPHom
  @test TX(BU, UU) == compose(TX(BU, A), TX(A, UU))
  @test TX(BU, UUU) == compose(TX(BU, A), TX(A, UUU))
  @test TX(BU, UUU) == compose(TX(BU, UU), TX(UU, UUU))
  @test WX(BUU, A) isa ModuleFPHom
  @test TX(BUU, A) isa ModuleFPHom
  @test TX(BUU, UU) == compose(TX(BUU, A), TX(A, UU))
  @test TX(BUU, UUU) == compose(TX(BUU, A), TX(A, UUU))
  @test TX(BUU, UUU) == compose(TX(BUU, UU), TX(UU, UUU))
  @test !iszero(TX(BUU, UUU))
end
