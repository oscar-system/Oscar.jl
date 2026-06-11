@testset "Space Germ constructors 1" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  X = spec(R, I*J, units_of(R))
  Y = spec(R, I, units_of(R))
  X0 = SpaceGerm(X,[0,0,0])
  X1 = SpaceGerm(X,[1,2,1])
  Y0 = SpaceGerm(Y,[0,0,0])
  Y1 = SpaceGerm(Y,[1,2,1])
  Z0 = HypersurfaceGerm(X,[0,0,0])
  @test defining_ring_element(Z0) isa elem_type(Oscar.localized_ring_type(ring_type(Z0)))
  Z1 = CompleteIntersectionGerm(X,[1,2,1])
  @test defining_ring_elements(Z1) isa Vector{elem_type(Oscar.localized_ring_type(ring_type(Z0)))}
  U0 = Oscar.MPolyComplementOfKPointIdeal(R,[0,0,0])
  U1 = Oscar.MPolyComplementOfKPointIdeal(R,[1,2,1])
  @test U0 == inverted_set(OO(X0))
  @test U1 == inverted_set(OO(X1))
  @test X0 == Y0
  @test X1 !=Y1
  @test isempty(Y1)
  @test milnor_number(Z0) == 1
  @test milnor_number(Z1) == 0
  XG = spec(R,I)
  @test milnor_number(XG) == 1
  K = ideal(R,[x^2+y^2-z^2,x*y])
  YG = spec(R,K)
  @test milnor_number(YG) == 5
  K = ideal(R,[x*y,x^2+y^2-z^2])
  ZG = spec(R,K)
  @test milnor_number(YG) == 5
end

@testset "SpaceGerms at geometric points" begin
  X = affine_space(QQ,3)
  OOX = coordinate_ring(X)
  (x,y,z) = gens(OOX)
  I = ideal(OOX,[x+y+z^2+1])
  J = ideal(OOX,[x+y,z^2+1])
  QI,_ = quo(OOX,I)
  QJ,_ = quo(OOX,J)
  U = complement_of_prime_ideal(ideal(OOX,[x,y,z^2+1]))
  RI,_ = localization(QI,U)
  RJ,_ = localization(QJ,U)
  Y = SpaceGerm(RI)
  @test OO(Y) == RI
#  HypersurfaceGerm(RI)      does not work without standard bases
#  CompleteIntersectionGerm(TI)     does not work without standard bases
  V = complement_of_prime_ideal(ideal(OOX,[x,z^2+1]))
  SI,_ = localization(QI,V)
  Z = SpaceGerm(SI)
  @test_broken dim(Z) == 1
end

@testset "Space Germ constructors AffineScheme-Ideal" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  J = ideal(R, [x-y])
  U0 = Oscar.MPolyComplementOfKPointIdeal(R,[0,0,0])
  U1 = Oscar.MPolyComplementOfKPointIdeal(R,[1,2,2])
  Xg = spec(R)
  Xgq = spec(R,J)
  Xl1 = spec(R, U0)
  Xl2 = spec(R, U1)
  Xlq = spec(R, J ,U0)
  Yg = SpaceGerm(Xg,ideal(OO(Xg),[x,y,z]))
  Ygq = SpaceGerm(Xgq, ideal(OO(Xgq),[x,y,z]))
  Yl1 = SpaceGerm(Xl1, ideal(OO(Xl1),[x,y,z]))
  @test_throws ErrorException("rings are incompatible") Yl2 = SpaceGerm(Xl2, ideal(OO(Xl1),[x,y,z]))
  Ylq = SpaceGerm(Xlq, ideal(OO(Xlq),[x,y,z]))
end

@testset "Space Germ constructors 2" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  S, (u,v,w) = QQ[:u, :v, :w]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  X = spec(R, I*J, units_of(R))
  Y = spec(R, I, units_of(R))
  X0 = SpaceGerm(X,ideal(R,[x,y,z]))
  X1 = SpaceGerm(X,ideal(R,[x-1,y-2,z-2]))
  Y0 = SpaceGerm(Y,ideal(R,[x,y,z]))
  Y1 = SpaceGerm(Y,ideal(R,[x-1,y-2,z-2]))
  U0 = Oscar.MPolyComplementOfKPointIdeal(R,[0,0,0])
  U1 = Oscar.MPolyComplementOfKPointIdeal(R,[1,2,2])
  @test U0 == inverted_set(OO(X0))
  @test U1 == inverted_set(OO(X1))
  @test X0 == Y0
  @test X1 !=Y1
  @test isempty(Y1)
end

@testset "Space Germ constructors 3" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  Q,_= quo(R,ideal(R,[z-y]))
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  Z = spec(R, ideal(R,[z-y]))
  X = subscheme(Z, I*J)
  Y = subscheme(Z, I)
  @test_throws ErrorException("rings are incompatible") Z0 = SpaceGerm(Z,ideal(Q,[x,y,z]))
  X0 = SpaceGerm(X,ideal(OO(X),[x,y,z]))
  @test_throws ArgumentError("Ideal does not describe finite set of points") X1 = SpaceGerm(X,ideal(OO(X),[x-1,y-2,z-1]))
  @test_throws ArgumentError("Ideal does not describe a single K-point") X1 = SpaceGerm(Z, ideal(OO(Z),[x^2-1,y-2,z-2]))
  X1 = SpaceGerm(X,ideal(OO(X),[x-1,y-2,z-2]))
  Y0 = SpaceGerm(Y,[0,0,0])
  Y1 = SpaceGerm(Y,[1,2,2])
end

@testset "germ_at_point constructors 1" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  S, (u,v,w) = QQ[:u, :v, :w]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  X = spec(R, I*J, units_of(R))
  @test_throws ArgumentError("Base rings must be the same.") spec(R,ideal(S,[u,v,w]))
  @test_throws ErrorException("rings are not compatible") SpaceGerm(X,ideal(S,[u,v,w]))
  @test_throws ErrorException("the number of variables in the ring does not coincide with the number of coordinates") SpaceGerm(X,[1,1,1,1])
  X0,phi0 = germ_at_point(X,ideal(R,[x,y,z]))
  X1 = SpaceGerm(X,ideal(R,[x,y,z]))
  X2,phi2 = germ_at_point(X,[0,0,0])
  U0 = Oscar.MPolyComplementOfKPointIdeal(R,[0,0,0])
  @test U0 == inverted_set(OO(X0))
  @test U0 == inverted_set(OO(X2))
  @test X0 == X1
  @test X0 == X2
  @test phi0 == phi2
end

@testset "germ_at_point constructors 2" begin
  R, (x,y,z) = QQ[:x, :y, :z]
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
  @test_throws ArgumentError("Ideal does not describe finite set of points") germ_at_point(Q2,(ideal(Q2,[x,y,z])))
  X3,phi3 = germ_at_point(Q3,(ideal(Q3,[x,y,z])))
  X4,phi4 = germ_at_point(OO(X0))
  X5,phi5 = germ_at_point(OO(X3))
  X6 = union(X2,X3)
end

@testset "space germ utilities" begin
  R, (x,y,z) = QQ[:x, :y, :z]
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
  @test defining_ideal(X0) != ideal(R,zero(R))
  @test defining_ideal(X0) == ideal(OO(X0),[0])
  @test defining_ideal(X2) == ideal(localized_ring(OO(X2)),[x,y])
  @test inverted_set(OO(ambient_germ(X2))) == Oscar.MPolyComplementOfKPointIdeal(R,[0,0,0])
  @test ambient_germ(X0) == X0
end

@testset "SpaceGerm from AffineScheme" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  K = ideal(R, [x,y])
  L = ideal(R, [x,y,z^2])
  X = spec(R, I*J, units_of(R))
  Y = spec(R, I, units_of(R))
  Z = spec(R, J ,units_of(R))
  W = spec(R, K ,units_of(R))
  V = spec(R, L ,units_of(R))
  SY = spec(R, ideal(R,[x,y,z]))
  X0 = SpaceGerm(X,[0,0,0])
  Y0 = SpaceGerm(Y,[0,0,0])
  Z0 = SpaceGerm(Z,[0,0,0])
  W0 = SpaceGerm(W,[0,0,0])
  SY0 = SpaceGerm(SY,[0,0,0])
  @test is_empty(Z0)
  @test is_subset(Z0,X0)
  @test is_subset(Y0,X0)
  @test is_subset(Z0,Y0)
  @test is_subset(X0,Y0)
  @test !is_subset(W0,X0)
  @test !is_subset(X0,Z0)
  @test Y0 == intersect(X0,Y0)
  V0=intersect(W0,X0)
  @test V0 == SpaceGerm(V,[0,0,0])
  @test SY0 == singular_locus(Y0)[1]
end
