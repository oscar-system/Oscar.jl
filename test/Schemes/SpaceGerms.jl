@testset "Space Germ constructors 1"
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
  @test_broken X0 == Y0
  @test_broken  X1 !=Y1
  @test isempty(Y1)
end

@testset "Space Germ constructors 2"
  R, (x,y,z) = QQ["x", "y", "z"]
  S, (u,v,w) = QQ["u", "v", "w"]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  X = Spec(R, I*J, units_of(R))
  Y = Spec(R, I, units_of(R))
  @test_throws Spec(R,ideal(S,[u,v,w]))
  X0 = SpaceGerm(X,ideal(R,[x,y,z]))
  X1 = SpaceGerm(X,ideal(R,[x-1,y-2,z-2]))
  Y0 = SpaceGerm(Y,ideal(R,[x,y,z]))
  Y1 = SpaceGerm(Y,ideal(R,[x-1,y-2,z-2]))
  U0 = MPolyComplementOfKPointIdeal(R,[0,0,0])
  U1 = MPolyComplementOfKPointIdeal(R,[1,2,2])
  @test U0 == inverted_set(OO(X0))
  @test U1 == inverted_set(OO(X1))
  @test_broken X0 == Y0
  @test_broken  X1 !=Y1
  @test isempty(Y1)
end

#@testset "Space Germ constructors 3" 
  R, (x,y,z) = QQ["x", "y", "z"]
  Q,_= quo(R,ideal(R,[z-y]))
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  Z = Spec(R, ideal(R,[z-y]))
  X = subscheme(Z, I*J)
  Y = subscheme(Z, I)
  @test_throws  X0 = SpaceGerm(X,ideal(Q,[x,y,z]))
  X0 = SpaceGerm(X,ideal(OO(X),[x,y,z]))
  @test_throws  X1 = SpaceGerm(X,ideal(OO(X),[x-1,y-2,z-1]))
  X1 = SpaceGerm(X,ideal(OO(X),[x-1,y-2,z-2]))
  Y0 = SpaceGerm(Y,[0,0,0])
  Y1 = SpaceGerm(Y,[1,2,2])
end

@testset "germ_at_point constructors 1"
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  X = Spec(R, I*J, units_of(R))
  @test_throws Spec(R,ideal(S,[u,v,w]))
  @test_throws SpaceGerm(X,ideal(S,[u,v,w]))
  @test_throws SpaceGerm(X,ideal(R)[x,y])
  @test_throws SpaceGerm(X,[1,1,1,1])
  X0,phi0 = germ_at_point(X,ideal(R,[x,y,z]))
  X1 = SpaceGerm(X,ideal(R,[x,y,z]))
  X2,phi2 = germ_at_point(X,[0,0,0])
  U0 = MPolyComplementOfKPointIdeal(R,[0,0,0])
  @test U0 == inverted_set(OO(X0))
  @test U0 == inverted_set(OO(X2))
  @test_broken X0 == X1
  @test_broken X0 == X2
  @test_broken phi0 == phi2
end

@testset "germ_at_point constructors 2";
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
  @test_throws germ_at_point(Q2,(ideal(Q,[x,y,z])))
  @test_throws germ_at_point(Q2,(ideal(Q2,[x,y,z])))
  X3,phi3 = germ_at_point(Q3,(ideal(Q3,[x,y,z])))
  X4,phi4 = germ_at_point(OO(X0))
  X5,phi5 = germ_at_point(OO(X3))
  X6 = union(X2,X3)
end

@testset "space germ utilities"
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
  @test_broken ambient_germ(X0) == X0
end