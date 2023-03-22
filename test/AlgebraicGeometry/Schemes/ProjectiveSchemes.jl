@testset "projective_schemes_1" begin

# test for relative projective space over a polynomial ring
R, (x,y) = QQ["x", "y"]
R_ext, _ = polynomial_ring(R, ["u", "v"])
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
R_ext, _ = polynomial_ring(QQ, ["u", "v"])
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
R_ext, _ = polynomial_ring(Q, ["u", "v"])
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
  S = ambient_coordinate_ring(P)
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
  S = ambient_coordinate_ring(P3)
  F = subscheme(P3,ideal(S,S[1]^4+S[2]^4+S[3]^4+S[4]^4))
  Fc = covered_scheme(F)
  U = patches(Fc)[1]
  V = patches(Fc)[2]
  oscar.intersect_in_covering(U,V,Fc[1])
end

@testset "Fermat lines" begin
  K,a = cyclotomic_field(8)
  P3 = projective_space(K,3)
  S = ambient_coordinate_ring(P3)
  @test Oscar.homogeneous_coordinate(P3,1) == S[1]
  F = subscheme(P3, ideal(S, S[1]^4 + S[2]^4 + S[3]^4 + S[4]^4))
  Fc = covered_scheme(F)
  U = patches(Fc)[1]
  V = patches(Fc)[2]
  oscar.intersect_in_covering(U,V,Fc[1]);
  line = subscheme(F, ideal(S, [S[1]+a*S[2],S[3]+a*S[4]]))
  groebner_basis(defining_ideal(line))

end

@testset "Issue #1580" begin 
  R,(x,) = polynomial_ring(GF(3),["x"])
  Rx,i = localization(R, x)
  x = Rx(x)
  P2 = projective_space(Rx, 2)
  affine_cone(P2) 
  @test covered_scheme(P2) isa CoveredScheme
end

@testset "affine cone" begin
  R,(x,) = polynomial_ring(GF(3),["x"])
  Rx,i = localization(R, x)
  x = Rx(x)
  Rq,j = quo(Rx,ideal(Rx,x))
  P2 = projective_space(Rx, 2)
  affine_cone(P2)
  @test covered_scheme(P2)  isa CoveredScheme
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

  @test ring_type(IP2_X) == typeof(ambient_coordinate_ring(IP2_X))
  @test ring_type(IP2_Y) == typeof(ambient_coordinate_ring(IP2_Y))
  @test ring_type(IP2_U) == typeof(ambient_coordinate_ring(IP2_U))
  @test ring_type(IP2_UY) == typeof(ambient_coordinate_ring(IP2_UY))

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

  IP2_Xh = subscheme(IP2_X, gens(ambient_coordinate_ring(IP2_X))[1])
  ProjectiveSchemeMor(IP2_Xh, IP2_X, gens(ambient_coordinate_ring(IP2_Xh)))
  IP2_Yh = subscheme(IP2_Y, gens(ambient_coordinate_ring(IP2_Y))[1])
  ProjectiveSchemeMor(IP2_Yh, IP2_Y, gens(ambient_coordinate_ring(IP2_Yh)))
  IP2_Uh = subscheme(IP2_U, gens(ambient_coordinate_ring(IP2_U))[1])
  ProjectiveSchemeMor(IP2_Uh, IP2_U, gens(ambient_coordinate_ring(IP2_Uh)))
  IP2_UYh = subscheme(IP2_UY, gens(ambient_coordinate_ring(IP2_UY))[1])
  ProjectiveSchemeMor(IP2_UYh, IP2_UY, gens(ambient_coordinate_ring(IP2_UYh)))
  IP2_Wh = subscheme(IP2_W, gens(ambient_coordinate_ring(IP2_W))[1])
  ProjectiveSchemeMor(IP2_Wh, IP2_W, gens(ambient_coordinate_ring(IP2_Wh)))

  incYtoX = inclusion_morphism(Y, X)
  h = hom(ambient_coordinate_ring(IP2_X), ambient_coordinate_ring(IP2_Y), pullback(incYtoX), gens(ambient_coordinate_ring(IP2_Y)))
  YtoX = ProjectiveSchemeMor(IP2_Y, IP2_X, 
                             h, 
                             incYtoX
                            );
  @test YtoX == inclusion_morphism(IP2_Y, IP2_X)
  IP2_Y2, map = fiber_product(incYtoX, IP2_X)
  h2 = hom(ambient_coordinate_ring(IP2_Y2), ambient_coordinate_ring(IP2_Y), gens(ambient_coordinate_ring(IP2_Y)))
  YtoY2 = ProjectiveSchemeMor(IP2_Y, IP2_Y2, h2)
  @test compose(YtoY2, map) == YtoX

  incUtoX = inclusion_morphism(U, X)
  h = hom(ambient_coordinate_ring(IP2_X), ambient_coordinate_ring(IP2_U), pullback(incUtoX), gens(ambient_coordinate_ring(IP2_U)))
  UtoX = ProjectiveSchemeMor(IP2_U, IP2_X, 
                             h, 
                             incUtoX
                            );
  @test UtoX == inclusion_morphism(IP2_U, IP2_X)
  IP2_U2, map = fiber_product(incUtoX, IP2_X)
  h2 = hom(ambient_coordinate_ring(IP2_U2), ambient_coordinate_ring(IP2_U), gens(ambient_coordinate_ring(IP2_U)))
  UtoU2 = ProjectiveSchemeMor(IP2_U, IP2_U2, h2)
  @test compose(UtoU2, map) == UtoX

  incUYtoY = inclusion_morphism(UY, Y)
  h = hom(ambient_coordinate_ring(IP2_Y), ambient_coordinate_ring(IP2_UY), pullback(incUYtoY), gens(ambient_coordinate_ring(IP2_UY)))
  UYtoY = ProjectiveSchemeMor(IP2_UY, IP2_Y, 
                              h, 
                              incUYtoY
                             );
  @test UYtoY == inclusion_morphism(IP2_UY, IP2_Y)
  IP2_UY2, map = fiber_product(incUYtoY, IP2_Y)
  h2 = hom(ambient_coordinate_ring(IP2_UY2), ambient_coordinate_ring(IP2_UY), gens(ambient_coordinate_ring(IP2_UY)))
  UYtoUY2 = ProjectiveSchemeMor(IP2_UY, IP2_UY2, h2)
  @test compose(UYtoUY2, map) == UYtoY

  incUYtoX = inclusion_morphism(UY, X)
  h = hom(ambient_coordinate_ring(IP2_X), ambient_coordinate_ring(IP2_UY), pullback(incUYtoX), gens(ambient_coordinate_ring(IP2_UY)))
  UYtoX = ProjectiveSchemeMor(IP2_UY, IP2_X, 
                              h, 
                              incUYtoX
                             );
  IP2_UY2, map = fiber_product(incUYtoX, IP2_X)
  h2 = hom(ambient_coordinate_ring(IP2_UY2), ambient_coordinate_ring(IP2_UY), gens(ambient_coordinate_ring(IP2_UY)))
  UYtoUY2 = ProjectiveSchemeMor(IP2_UY, IP2_UY2, h2)
  @test compose(UYtoUY2, map) == UYtoX

  @test UYtoX == inclusion_morphism(IP2_UY, IP2_X)
  @test compose(UYtoY, YtoX) == UYtoX

  WW = hypersurface_complement(ambient_scheme(W), [x-y])
  phi = MapFromFunc(restriction_map(W, WW), OO(W), OO(WW))
  IP2_WW, map = fiber_product(phi, IP2_W)
  @test base_scheme(IP2_WW) == WW
  @test !(base_scheme(IP2_WW) === WW)

  IP2_QQ = projective_space(QQ, 2)
  id = identity_map(IP2_QQ)
  @test id == identity_map(IP2_QQ)
end

@testset "warham_preparations" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  S, _ = grade(R)
  X = ProjectiveScheme(S)
  I = ideal(S, x^2 - y*z)
  Q, _ = quo(S, I)
  C = ProjectiveScheme(Q)
  @test affine_algebra(C) === Q
  @test dim(C) == 1
  @test degree(C) == 2
  @test is_smooth(C)
  @test arithmetic_genus(C) == 0

  R, (x,y,z,w) = QQ["x", "y", "z", "w"]
  S, _ = grade(R)
  I = ideal(S, [x^4 + y^4 + z^4 + w^4])
  Q, _ = quo(S, I)
  Y = ProjectiveScheme(Q)
  @test affine_algebra(Y) === Q
  @test dim(Y) == 2
  @test degree(Y) == 4
  @test is_smooth(Y)
  @test arithmetic_genus(Y) == 1


  Base.:/(x::RingElem, y::RingElement) = divexact(x, y; check=true)

  function abelian_d10_pi6(;ring::MPolyRing=graded_polynomial_ring(GF(31991), ["x", "y", "z", "u", "v"])[1])
    (x,y,z,u,v) = gens(ring)
    F = coefficient_ring(ring)
    I = ideal(ring,
    [-x^3*y^2*z-7318*x^4*z^2-x^2*z^3*u-x*y^3*u^2-14636*x^2*y*z*u^2-y*z^2*u^3-7318*y^2*u^4+y^5*v-8856*x^2*y*z^2*v+z^5*v-8856*y^2*z*u^2*v-7318*x^3*y*v^2-F(99)/52*x*z^3*v^2-F(99)/52*y^3*u*v^2-4*x*y*z*u*v^2-7318*z*u^3*v^2-8856*x*y^2*v^3-8856*z^2*u*v^3-F(99)/52*y*z*v^4+7318*x*u*v^4+v^6,-x*y^2*z*u^2-7318*x^2*z^2*u^2-z^3*u^3-7318*y*z*u^4+y^4*z*v+7318*x*y^2*z^2*v-7318*z^4*u*v-8856*y*z^2*u^2*v+x^2*y^2*v^2-F(99)/52*y^2*z*u*v^2-x*z^2*u*v^2+7318*x*y*u^2*v^2+7318*y^3*v^3+y*u*v^4,-x*y^2*u^3-7318*x^2*z*u^3-z^2*u^4-7318*y*u^5+y^4*u*v+7318*x*y^2*z*u*v-7318*z^3*u^2*v-8856*y*z*u^3*v-7318*x^3*u*v^2-F(99)/52*y^2*u^2*v^2-2*x*z*u^2*v^2-y^2*z*v^3-7318*x*z^2*v^3-8856*x*y*u*v^3-x^2*v^4-7318*y*v^5,y^5*z+7318*x*y^3*z^2-7318*y*z^4*u+x^3*y*u^2-8856*y^2*z^2*u^2-F(99)/52*x*y*z*u^3+7318*x^2*u^4+z*u^5+x^2*y^3*v-F(99)/52*y^3*z*u*v-x*y*z^2*u*v+14636*x*y^2*u^2*v+7318*z^2*u^3*v+7318*y^4*v^2+x*u^3*v^2+y^2*u*v^3,y^4*z^2+7318*x*y^2*z^3-7318*z^5*u+x^3*z*u^2-8856*y*z^3*u^2-F(99)/52*x*z^2*u^3-x*y*u^4-7318*u^6+2*x^2*y^2*z*v+7318*x^3*z^2*v-F(99)/52*y^2*z^2*u*v+y^3*u^2*v-2719*x*y*z*u^2*v-8856*z*u^4*v+x^4*v^2+7318*y^3*z*v^2-F(99)/52*x^2*z*u*v^2-F(99)/52*y*u^3*v^2+7318*x^2*y*v^3-8856*x*u^2*v^3-7318*u*v^5,-x*y^2*z^3-7318*x^2*z^4-z^5*u-7318*y*z^3*u^2-x^3*z^2*v+F(99)/52*x*z^3*u*v+x*y*z*u^2*v+7318*z*u^4*v-y^3*z*v^2-14636*x*y*z^2*v^2+8856*z^2*u^2*v^2-x^2*y*v^3+F(99)/52*y*z*u*v^3-7318*x*u^2*v^3-7318*y^2*v^4-u*v^5,y^3*z^3+7318*x*y*z^4-x^4*y*u+7318*y^4*z*u+8856*x*y^2*z^2*u+F(99)/52*x^2*y*z*u^2-7318*x^3*u^3-x*z*u^4+x^2*y*z^2*v-7318*x^2*y^2*u*v+y^2*z*u^2*v-7318*x*z^2*u^2*v+7318*y^2*z^2*v^2-x^2*u^2*v^2,-y^2*z^4-7318*x*z^5+x^4*z*u-7318*y^3*z^2*u-8856*x*y*z^3*u-F(99)/52*x^2*z^2*u^2-x^2*y*u^3-7318*x*u^5-x^2*z^3*v+7318*x^2*y*z*u*v-2*y*z^2*u^2*v-7318*y^2*u^3*v-8856*x*z*u^3*v-7318*y*z^3*v^2-u^4*v^2-7318*z*u^2*v^3,y^6+7318*x*y^4*z-7318*y^2*z^3*u-7318*x^4*u^2-8856*y^3*z*u^2-F(99)/52*x*y^2*u^3-x^2*z*u^3+y*u^5-7318*x^3*y^2*v-F(99)/52*y^4*u*v-4*x*y^2*z*u*v-14636*x^2*z^2*u*v-8856*x^2*y*u^2*v-z^3*u^2*v-8856*x*y^3*v^2-7318*z^4*v^2-x^3*u*v^2-8856*y*z^2*u*v^2-F(99)/52*y^2*z*v^3-x*z^2*v^3+y*v^5,x^4*y^2-7318*y^5*z-8856*x*y^3*z^2-7318*z^6-F(99)/52*x^2*y^2*z*u+x^3*z^2*u-8856*y*z^4*u+7318*x^3*y*u^2-F(99)/52*x*z^3*u^2-7318*z*u^5+7318*x^2*y^3*v-F(99)/52*y^2*z^3*v-x*z^4*v-2719*x*y*z^2*u*v-8856*z^2*u^3*v+2*x^2*y*u*v^2-F(99)/52*y*z*u^2*v^2+7318*x*u^3*v^2+y*z^2*v^3+7318*y^2*u*v^3+u^2*v^4,-x^3*y^3-7318*x^4*y*z-x^2*y*z^2*u-7318*x^2*y^2*u^2-7318*x*y^4*v-8856*x^2*y^2*z*v+y*z^4*v+7318*y^2*z^2*u*v-F(99)/52*x*y*z^2*v^2-x*y^2*u*v^2+7318*x^2*z*u*v^2+z^2*u^2*v^2+7318*z^3*v^3+x*z*v^4,-x^2*y^4-7318*x^3*y^2*z-2*x*y^2*z^2*u-7318*x^2*z^3*u-7318*x*y^3*u^2-z^4*u^2-7318*y*z^2*u^3-7318*y^5*v-8856*x*y^3*z*v-7318*z^5*v-8856*y*z^3*u*v-F(99)/52*y^2*z^2*v^2-x*z^3*v^2-y^3*u*v^2+7318*x*y*z*u*v^2+y*z*v^4,-x*y^5-7318*x^2*y^3*z-y^3*z^2*u-7318*y^4*u^2+7318*x^4*y*v+F(99)/52*x*y^3*u*v+x^2*y*z*u*v-y^2*u^3*v+8856*x^2*y^2*v^2-y*z^3*v^2-14636*y^2*z*u*v^2+F(99)/52*x*y*z*v^3-7318*x^2*u*v^3-z*u^2*v^3-7318*z^2*v^4-x*v^5,-x*y^4*u-7318*x^2*y^2*z*u-y^2*z^2*u^2-7318*y^3*u^3+7318*x^4*u*v+F(99)/52*x*y^2*u^2*v+x^2*z*u^2*v-y*u^4*v+x*y^2*z*v^2+7318*x^2*z^2*v^2+8856*x^2*y*u*v^2-7318*y*z*u^2*v^2+x^3*v^3+7318*x*y*v^4,-x^2*y^2*z^2-7318*x^3*z^3-x*z^4*u-7318*x*y*z^2*u^2-x^4*z*v+F(99)/52*x^2*z^2*u*v+x^2*y*u^2*v+7318*x*u^4*v-7318*x^2*y*z*v^2+y*z^2*u*v^2+7318*y^2*u^2*v^2+8856*x*z*u^2*v^2+u^3*v^3+7318*z*u*v^4,-x*y^3*z-7318*x^2*y*z^2-y*z^3*u-7318*y^2*z*u^2-x^3*y*v+F(99)/52*x*y*z*u*v-7318*x^2*u^2*v-z*u^3*v-7318*x*y^2*v^2-7318*z^2*u*v^2-x*u*v^3,7318*x^5+7318*y^5+8856*x*y^3*z+7318*z^5+F(99)/52*x^2*y^2*u+8856*y*z^3*u+F(99)/52*x*z^2*u^2+7318*u^5+8856*x^3*y*v+F(99)/52*y^2*z^2*v-4599*x*y*z*u*v+8856*z*u^3*v+F(99)/52*x^2*z*v^2+F(99)/52*y*u^2*v^2+8856*x*u*v^3+7318*v^5,-x^2*y^2*u-7318*x^3*z*u-x*z^2*u^2-7318*x*y*u^3-y^2*z^2*v-7318*x*z^3*v-7318*y^3*u*v-8856*x*y*z*u*v-x^2*z*v^2-y*u^2*v^2-7318*y*z*v^3,-x^5-y^5+8856*x^2*y*z^2-z^5+F(99)/52*x^3*z*u+8856*y^2*z*u^2+F(99)/52*x*y*u^3-u^5+F(99)/52*x*z^3*v+F(99)/52*y^3*u*v+5*x*y*z*u*v+8856*x^2*u^2*v+8856*x*y^2*v^2+8856*z^2*u*v^2+F(99)/52*y*z*v^3-v^5]);
    return I
  end

  I = abelian_d10_pi6()
  Q, _ = quo(base_ring(I), I)
  Z = ProjectiveScheme(Q)

  @test dim(Y) == 2
  @test degree(Y) == 4
  # @test is_smooth(Y) # too expensive!
  @test arithmetic_genus(Y) == 1
end
