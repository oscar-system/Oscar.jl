@testset "projective_schemes_1" begin

# test for relative projective space over a polynomial ring
R, (x,y) = QQ["x", "y"]
R_ext, _ = PolynomialRing(R, ["u", "v"])
S, (u,v) = grade(R_ext, [1,1])

I = ideal(S, [x*v - y*u])
X = ProjectiveScheme(S, I)
CX = affine_cone(X)
p = covered_projection_to_base(X)
@test OO(CX).(homogeneous_coordinates(X)) == [homog_to_frac(X)(g) for g in gens(S)]
hc = homogeneous_coordinates(X)
frac_to_homog_pair(X)(hc[1]*hc[2])

phi = ProjectiveSchemeMor(X, X, [-u, -v])

g = map_on_affine_cones(phi)
@test is_well_defined(phi)

# test for projective space over a field
R_ext, _ = PolynomialRing(QQ, ["u", "v"])
S, (u,v) = grade(R_ext, [1,1])

I = ideal(S, [u])
X = ProjectiveScheme(S, I)
CX = affine_cone(X)
@test OO(CX).(homogeneous_coordinates(X)) == [homog_to_frac(X)(g) for g in gens(S)]
hc = homogeneous_coordinates(X)
frac_to_homog_pair(X)(hc[1]*hc[2])

phi = ProjectiveSchemeMor(X, X, [u^2, v^2])

g = map_on_affine_cones(phi)

@test is_well_defined(phi)

# test for relative projective space over MPolyQuoLocalizedRings
Y = Spec(R)
Q = OO(Y)
R_ext, _ = PolynomialRing(Q, ["u", "v"])
S, (u,v) = grade(R_ext, [1,1])
X = ProjectiveScheme(S)

phi = ProjectiveSchemeMor(X, X, [u^2, v^2])

g = map_on_affine_cones(phi)

@test is_well_defined(phi)
end

@testset "projective_schemes_2" begin
  R, (x, y, z) = QQ["x", "y", "z"]
  I = ideal(R, [x^2-y*z])
  X = Spec(R, I)
  U = SpecOpen(X, [x, y])
  P = projective_space(OO(U), 1)
  S = ambient_ring(P)
  Y = subscheme(P, [x*S[1]- y*S[2], z*S[1] - x*S[2]])
  C = affine_cone(Y)
  phi = homog_to_frac(Y)
  s1 = phi(S[2])
  s0 = phi(S[1])
  a = OO(U)([x//y, z//x])
  @test pullback(projection_to_base(Y))(a)*s0 == s1
end


@testset "projective schemes as covered schemes" begin
  P3 = projective_space(QQ,3)
  S = ambient_ring(P3)
  F = subscheme(P3,ideal(S,S[1]^4+S[2]^4+S[3]^4+S[4]^4))
  Fc = as_covered_scheme(F)
  U = patches(Fc)[1]
  V = patches(Fc)[2]
  oscar.intersect_in_covering(U,V,Fc[1])
end

@testset "Fermat lines" begin
  K,a = cyclotomic_field(8)
  P3 = projective_space(K,3)
  S = ambient_ring(P3)
  @test Oscar.homogeneous_coordinate(P3,1) == S[1]
  F = subscheme(P3, ideal(S, S[1]^4 + S[2]^4 + S[3]^4 + S[4]^4))
  Fc = as_covered_scheme(F)
  U = patches(Fc)[1]
  V = patches(Fc)[2]
  oscar.intersect_in_covering(U,V,Fc[1]);
  line = subscheme(F, ideal(S, [S[1]+a*S[2],S[3]+a*S[4]]))
  groebner_basis(defining_ideal(line))

end

@testset "Issue #1580" begin 
  R,(x,) = PolynomialRing(GF(3),["x"])
  Rx,i = localization(R, x)
  x = Rx(x)
  P2 = projective_space(Rx, 2)
  affine_cone(P2) 
  @test_broken as_covered_scheme(P2)
end

@testset "affine cone" begin
  R,(x,) = PolynomialRing(GF(3),["x"])
  Rx,i = localization(R, x)
  x = Rx(x)
  Rq,j = quo(Rx,ideal(Rx,x))
  P2 = projective_space(Rx, 2)
  affine_cone(P2)
  @test_broken as_covered_scheme(P2)  # too exotic
end

@testset "morphisms of projective schemes I" begin
  R, (x,y) = QQ["x", "y"]

  IP2_X = projective_space(R, 2, var_name="u")
  projective_space(R, ["u", "v", "w"])
  X = base_scheme(IP2_X)
  projective_space(X, ["u", "v", "w"])

  inc = ClosedEmbedding(X, ideal(OO(X), [x^2 - y^2]))
  Y = domain(inc)

  IP2_Y = projective_space(Y, 2, var_name="v")
  projective_space(Y, ["u", "v", "w"])
  projective_space(OO(Y), ["u", "v", "w"])

  U = PrincipalOpenSubset(X, x)

  IP2_U = projective_space(OO(U), 2, var_name="u")
  projective_space(U, ["u", "v", "w"])
  projective_space(OO(U), ["u", "v", "w"])
  set_base_scheme!(IP2_U, U)
  @test base_scheme(IP2_U) === U

  UY = intersect(U, Y)

  IP2_UY = projective_space(UY, 2, var_name="v")
  projective_space(UY, ["u", "v", "w"])
  projective_space(OO(UY), ["u", "v", "w"])

  @test projective_scheme_type(OO(X)) == typeof(IP2_X)
  @test projective_scheme_type(OO(Y)) == typeof(IP2_Y)
  @test projective_scheme_type(OO(U)) == typeof(IP2_U)
  @test projective_scheme_type(OO(UY)) == typeof(IP2_UY)

  @test projective_scheme_type(X) == typeof(IP2_X)
  @test projective_scheme_type(Y) == typeof(IP2_Y)
  @test projective_scheme_type(U) == typeof(IP2_U)
  @test projective_scheme_type(UY) == typeof(IP2_UY)

  @test base_ring_type(IP2_X) == typeof(OO(X))
  @test base_ring_type(IP2_Y) == typeof(OO(Y))
  @test base_ring_type(IP2_U) == typeof(OO(U))
  @test base_ring_type(IP2_UY) == typeof(OO(UY))

  @test ring_type(IP2_X) == typeof(ambient_ring(IP2_X))
  @test ring_type(IP2_Y) == typeof(ambient_ring(IP2_Y))
  @test ring_type(IP2_U) == typeof(ambient_ring(IP2_U))
  @test ring_type(IP2_UY) == typeof(ambient_ring(IP2_UY))

  CX = affine_cone(IP2_X)
  CY = affine_cone(IP2_Y)
  CU = affine_cone(IP2_U)
  CUY = affine_cone(IP2_UY)

  pCX = projection_to_base(IP2_X)
  pCY = projection_to_base(IP2_Y)
  pCU = projection_to_base(IP2_U)
  pCUY = projection_to_base(IP2_UY)

  @test domain(pCX) == CX
  @test domain(pCY) == CY
  @test domain(pCU) == CU
  @test domain(pCUY) == CUY
  @test codomain(pCX) == X
  @test codomain(pCY) == Y
  @test codomain(pCU) == U
  @test codomain(pCUY) == UY

  @test sprint(show, IP2_X) isa String
  @test sprint(show, IP2_Y) isa String
  @test sprint(show, IP2_U) isa String
  @test sprint(show, IP2_UY) isa String


  W = SpecOpen(UY, [x-1, y-1])
  IP2_W = projective_space(W, 2, var_name="w")
  CW = affine_cone(IP2_W)
  pCW = projection_to_base(IP2_W)

  WY = SpecOpen(Y, [x-y, x+y-2])
  WtoWY = inclusion_morphism(W, WY)
  IP2_WY = projective_space(WY, 2, var_name="w")
  _, map = fiber_product(pullback(WtoWY), IP2_WY)

  map_on_affine_cones(map)

  IP2_Xh = subscheme(IP2_X, gens(ambient_ring(IP2_X))[1])
  ProjectiveSchemeMor(IP2_Xh, IP2_X, gens(ambient_ring(IP2_Xh)))
  IP2_Yh = subscheme(IP2_Y, gens(ambient_ring(IP2_Y))[1])
  ProjectiveSchemeMor(IP2_Yh, IP2_Y, gens(ambient_ring(IP2_Yh)))
  IP2_Uh = subscheme(IP2_U, gens(ambient_ring(IP2_U))[1])
  ProjectiveSchemeMor(IP2_Uh, IP2_U, gens(ambient_ring(IP2_Uh)))
  IP2_UYh = subscheme(IP2_UY, gens(ambient_ring(IP2_UY))[1])
  ProjectiveSchemeMor(IP2_UYh, IP2_UY, gens(ambient_ring(IP2_UYh)))
  IP2_Wh = subscheme(IP2_W, gens(ambient_ring(IP2_W))[1])
  ProjectiveSchemeMor(IP2_Wh, IP2_W, gens(ambient_ring(IP2_Wh)))

  incYtoX = inclusion_morphism(Y, X)
  h = hom(ambient_ring(IP2_X), ambient_ring(IP2_Y), pullback(incYtoX), gens(ambient_ring(IP2_Y)))
  YtoX = ProjectiveSchemeMor(IP2_Y, IP2_X, 
                             h, 
                             incYtoX
                            );
  @test YtoX == inclusion_morphism(IP2_Y, IP2_X)
  IP2_Y2, map = fiber_product(incYtoX, IP2_X)
  h2 = hom(ambient_ring(IP2_Y2), ambient_ring(IP2_Y), gens(ambient_ring(IP2_Y)))
  YtoY2 = ProjectiveSchemeMor(IP2_Y, IP2_Y2, h2)
  @test compose(YtoY2, map) == YtoX

  incUtoX = inclusion_morphism(U, X)
  h = hom(ambient_ring(IP2_X), ambient_ring(IP2_U), pullback(incUtoX), gens(ambient_ring(IP2_U)))
  UtoX = ProjectiveSchemeMor(IP2_U, IP2_X, 
                             h, 
                             incUtoX
                            );
  @test UtoX == inclusion_morphism(IP2_U, IP2_X)
  IP2_U2, map = fiber_product(incUtoX, IP2_X)
  h2 = hom(ambient_ring(IP2_U2), ambient_ring(IP2_U), gens(ambient_ring(IP2_U)))
  UtoU2 = ProjectiveSchemeMor(IP2_U, IP2_U2, h2)
  @test compose(UtoU2, map) == UtoX

  incUYtoY = inclusion_morphism(UY, Y)
  h = hom(ambient_ring(IP2_Y), ambient_ring(IP2_UY), pullback(incUYtoY), gens(ambient_ring(IP2_UY)))
  UYtoY = ProjectiveSchemeMor(IP2_UY, IP2_Y, 
                              h, 
                              incUYtoY
                             );
  @test UYtoY == inclusion_morphism(IP2_UY, IP2_Y)
  IP2_UY2, map = fiber_product(incUYtoY, IP2_Y)
  h2 = hom(ambient_ring(IP2_UY2), ambient_ring(IP2_UY), gens(ambient_ring(IP2_UY)))
  UYtoUY2 = ProjectiveSchemeMor(IP2_UY, IP2_UY2, h2)
  @test compose(UYtoUY2, map) == UYtoY

  incUYtoX = inclusion_morphism(UY, X)
  h = hom(ambient_ring(IP2_X), ambient_ring(IP2_UY), pullback(incUYtoX), gens(ambient_ring(IP2_UY)))
  UYtoX = ProjectiveSchemeMor(IP2_UY, IP2_X, 
                              h, 
                              incUYtoX
                             );
  IP2_UY2, map = fiber_product(incUYtoX, IP2_X)
  h2 = hom(ambient_ring(IP2_UY2), ambient_ring(IP2_UY), gens(ambient_ring(IP2_UY)))
  UYtoUY2 = ProjectiveSchemeMor(IP2_UY, IP2_UY2, h2)
  @test compose(UYtoUY2, map) == UYtoX

  @test UYtoX == inclusion_morphism(IP2_UY, IP2_X)
  @test compose(UYtoY, YtoX) == UYtoX

  WW = hypersurface_complement(ambient(W), [x-y])
  phi = MapFromFunc(restriction_map(W, WW), OO(W), OO(WW))
  IP2_WW, map = fiber_product(phi, IP2_W)
  @test base_scheme(IP2_WW) == WW
  @test !(base_scheme(IP2_WW) === WW)
end
