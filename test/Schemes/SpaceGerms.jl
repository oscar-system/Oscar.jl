@testset "Space Germ constructors 1" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  X = Spec(R, I*J, units_of(R))
  Y = Spec(R, I, units_of(R))
  X0 = SpaceGerm(X,[0,0,0])
  X1 = SpaceGerm(X,[1,2,1])
  Y0 = SpaceGerm(Y,[0,0,0])
  Y1 = SpaceGerm(Y,[1,2,1])
  U0 = MPolyComplementOfKPointIdeal(R,[0,0,0])
  U1 = MPolyComplementOfKPointIdeal(R,[1,2,1])
  @test U0 == inverted_set(OO(X0))
  @test U1 == inverted_set(OO(X1))
  @test X0 == Y0
  @test X1 !=Y1
  @test isempty(Y1)
end

@testset "Space Germ constructors Spec-Ideal" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  J = ideal(R, [x-y])
  U0 = MPolyComplementOfKPointIdeal(R,[0,0,0])
  U1 = MPolyComplementOfKPointIdeal(R,[1,2,2])
  Xg = Spec(R)
  Xgq = Spec(R,J)
  Xl1 = Spec(R, U0)
  Xl2 = Spec(R, U1)
  Xlq = Spec(R, J ,U0)
  Yg = SpaceGerm(Xg,ideal(OO(Xg),[x,y,z]))
  Ygq = SpaceGerm(Xgq, ideal(OO(Xgq),[x,y,z]))
  Yl1 = SpaceGerm(Xl1, ideal(OO(Xl1),[x,y,z]))
  @test_throws ErrorException("rings are incompatible") Yl2 = SpaceGerm(Xl2, ideal(OO(Xl1),[x,y,z]))
  Ylq = SpaceGerm(Xlq, ideal(OO(Xlq),[x,y,z]))
end

@testset "Space Germ constructors 2" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  S, (u,v,w) = QQ["u", "v", "w"]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  X = Spec(R, I*J, units_of(R))
  Y = Spec(R, I, units_of(R))
  X0 = SpaceGerm(X,ideal(R,[x,y,z]))
  X1 = SpaceGerm(X,ideal(R,[x-1,y-2,z-2]))
  Y0 = SpaceGerm(Y,ideal(R,[x,y,z]))
  Y1 = SpaceGerm(Y,ideal(R,[x-1,y-2,z-2]))
  U0 = MPolyComplementOfKPointIdeal(R,[0,0,0])
  U1 = MPolyComplementOfKPointIdeal(R,[1,2,2])
  @test U0 == inverted_set(OO(X0))
  @test U1 == inverted_set(OO(X1))
  @test X0 == Y0
  @test X1 !=Y1
  @test isempty(Y1)
end

@testset "Space Germ constructors 3" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  Q,_= quo(R,ideal(R,[z-y]))
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  Z = Spec(R, ideal(R,[z-y]))
  X = subscheme(Z, I*J)
  Y = subscheme(Z, I)
  @test_throws ErrorException("rings are incompatible") Z0 = SpaceGerm(Z,ideal(Q,[x,y,z]))
  X0 = SpaceGerm(X,ideal(OO(X),[x,y,z]))
  @test_throws ErrorException("Ideal does not describe finite set of points") X1 = SpaceGerm(X,ideal(OO(X),[x-1,y-2,z-1]))
  @test_throws ErrorException("Ideal does not describe a single K-point") X1 = SpaceGerm(Z, ideal(OO(Z),[x^2-1,y-2,z-2]))
  X1 = SpaceGerm(X,ideal(OO(X),[x-1,y-2,z-2]))
  Y0 = SpaceGerm(Y,[0,0,0])
  Y1 = SpaceGerm(Y,[1,2,2])
end

@testset "germ_at_point constructors 1" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  S, (u,v,w) = QQ["u", "v", "w"]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  X = Spec(R, I*J, units_of(R))
  @test_throws AssertionError("base_ring(I) === R") Spec(R,ideal(S,[u,v,w]))
  @test_throws ErrorException("rings are not compatible") SpaceGerm(X,ideal(S,[u,v,w]))
  @test_throws ErrorException("the number of variables in the ring does not coincide with the number of coordinates") SpaceGerm(X,[1,1,1,1])
  X0,phi0 = germ_at_point(X,ideal(R,[x,y,z]))
  X1 = SpaceGerm(X,ideal(R,[x,y,z]))
  X2,phi2 = germ_at_point(X,[0,0,0])
  U0 = MPolyComplementOfKPointIdeal(R,[0,0,0])
  @test U0 == inverted_set(OO(X0))
  @test U0 == inverted_set(OO(X2))
  @test X0 == X1
  @test X0 == X2
  @test phi0 == phi2
end

@testset "germ_at_point constructors 2" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  K = ideal(R, [x,y])
  Q,_ = quo(R, I*J)
  Q1,_ = quo(R,I)
  Q2,_ = quo(R,J)
  Q3,_ = quo(R,K)
  X0,phi0 = germ_at_point(R,[0,0,0])
  X1,phi1 = germ_at_point(R,ideal(R,[x,y,z]))
  X2,phi2 = germ_at_point(Q1,[0,0,0])
  @test_throws ErrorException("rings are incompatible") germ_at_point(Q2,(ideal(Q,[x,y,z])))
  @test_throws ErrorException("Ideal does not describe finite set of points") germ_at_point(Q2,(ideal(Q2,[x,y,z])))
  X3,phi3 = germ_at_point(Q3,(ideal(Q3,[x,y,z])))
  X4,phi4 = germ_at_point(OO(X0))
  X5,phi5 = germ_at_point(OO(X3))
  X6 = union(X2,X3)
end

@testset "space germ utilities" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x, y])
  Q1,_ = quo(R,I)
  Q2,_ = quo(R,J)
  X0,phi0 = germ_at_point(R,[0,0,0])
  X1,phi1 = germ_at_point(Q1,[0,0,0])
  X2,phi2 = germ_at_point(Q2,[0,0,0])
  @test point(X0)==[0,0,0]
  Y = representative(X0)
  Y1 = representative(X1)
  @test ideal(X0) != ideal(R,zero(R))
  @test ideal(X0) == ideal(OO(X0),[0])
  @test ideal(X2) == ideal(localized_ring(OO(X2)),[x,y])
  @test inverted_set(OO(ambient_germ(X2))) == MPolyComplementOfKPointIdeal(R,[0,0,0])
  @test ambient_germ(X0) == X0
end

@testset "SpaceGerm from Spec" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  K = ideal(R, [x,y])
  L = ideal(R, [x,y,z^2])
  X = Spec(R, I*J, units_of(R))
  Y = Spec(R, I, units_of(R))
  Z = Spec(R, J ,units_of(R))
  W = Spec(R, K ,units_of(R))
  V = Spec(R, L ,units_of(R))
  SY = Spec(R, ideal(R,[x,y,z]))
  X0 = SpaceGerm(X,[0,0,0])
  Y0 = SpaceGerm(Y,[0,0,0])
  Z0 = SpaceGerm(Z,[0,0,0])
  W0 = SpaceGerm(W,[0,0,0])
  SY0 = SpaceGerm(SY,[0,0,0])
  @test isempty(Z0)
  @test issubset(Z0,X0)
  @test issubset(Y0,X0)
  @test issubset(Z0,Y0)
  @test issubset(X0,Y0)
  @test !issubset(W0,X0)
  @test !issubset(X0,Z0)
  @test Y0 == intersect(X0,Y0)
  V0=intersect(W0,X0)
  @test V0 == SpaceGerm(V,[0,0,0])
  @test SY0 == singular_locus(Y0)
end
