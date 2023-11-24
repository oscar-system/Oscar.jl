@testset "projective_schemes_1" begin
  # test for relative projective space over a polynomial ring
  R, (x,y) = QQ["x", "y"]
  S, (u,v) = graded_polynomial_ring(R, ["u", "v"])

  I = ideal(S, [x*v - y*u])
  X = ProjectiveScheme(S, I)
  CX, id = affine_cone(X)
  p = covered_projection_to_base(X)
  @test OO(CX).(Oscar.homogeneous_coordinates_on_affine_cone(X)) == [id(g) for g in gens(S)]
  hc = Oscar.homogeneous_coordinates_on_affine_cone(X)

  phi = ProjectiveSchemeMor(X, X, [-u, -v])

  g = map_on_affine_cones(phi)
  #@test is_well_defined(phi) # deprecated

  # test for projective space over a field
  S, (u,v) = graded_polynomial_ring(QQ, ["u", "v"])

  I = ideal(S, [u])
  X = ProjectiveScheme(S, I)
  CX, id = affine_cone(X)
  @test OO(CX).(Oscar.homogeneous_coordinates_on_affine_cone(X)) == [id(g) for g in gens(S)]
  hc = Oscar.homogeneous_coordinates_on_affine_cone(X)

  phi = ProjectiveSchemeMor(X, X, [u^2, v^2])

  g = map_on_affine_cones(phi)

  #@test is_well_defined(phi) # deprecated

  # test for relative projective space over MPolyQuoLocalizedRings
  Y = Spec(R)
  Q = OO(Y)
  S, (u,v) = graded_polynomial_ring(Q, ["u", "v"])
  X = ProjectiveScheme(S)

  phi = ProjectiveSchemeMor(X, X, [u^2, v^2])

  g = map_on_affine_cones(phi)

  #@test is_well_defined(phi) # deprecated
end

@testset "projective_schemes_2" begin
  R, (x, y, z) = QQ["x", "y", "z"]
  I = ideal(R, [x^2-y*z])
  X = Spec(R, I)
  U = SpecOpen(X, [x, y])
  P = projective_space(OO(U), 1)
  S = homogeneous_coordinate_ring(P)
  Y = subscheme(P, [OO(U)(x)*S[1]- OO(U)(y)*S[2], OO(U)(z)*S[1] - OO(U)(x)*S[2]]) # Coercion needs to be carried out manually.
  C, phi = affine_cone(Y)
  s1 = phi(S[2])
  s0 = phi(S[1])
  a = OO(U)([x//y, z//x])
  pullback(projection_to_base(Y))(a)*s0 == s1
  @test pullback(projection_to_base(Y))(a)*s0 == s1
end


@testset "projective schemes as covered schemes" begin
  P3 = projective_space(QQ,3)
  S = homogeneous_coordinate_ring(P3)
  F = subscheme(P3,ideal(S,S[1]^4+S[2]^4+S[3]^4+S[4]^4))
  Fc = covered_scheme(F)
  U = patches(Fc)[1]
  V = patches(Fc)[2]
  Oscar.intersect_in_covering(U,V,Fc[1])
end

@testset "Fermat lines" begin
  K,a = cyclotomic_field(8)
  P3 = projective_space(K,3)
  S = homogeneous_coordinate_ring(P3)
  F = subscheme(P3, ideal(S, S[1]^4 + S[2]^4 + S[3]^4 + S[4]^4))
  Fc = covered_scheme(F)
  U = patches(Fc)[1]
  V = patches(Fc)[2]
  Oscar.intersect_in_covering(U,V,Fc[1]);
  SF = homogeneous_coordinate_ring(F)
  line = subscheme(F, ideal(SF, [SF[1]+a*SF[2],SF[3]+a*SF[4]]))
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

#  @test projective_scheme_type(OO(X)) == typeof(IP2_X)
#  @test projective_scheme_type(OO(Y)) == typeof(IP2_Y)
#  @test projective_scheme_type(OO(U)) == typeof(IP2_U)
#  @test projective_scheme_type(OO(UY)) == typeof(IP2_UY)
# 
#  @test projective_scheme_type(X) == typeof(IP2_X)
#  @test projective_scheme_type(Y) == typeof(IP2_Y)
#  @test projective_scheme_type(U) == typeof(IP2_U)
#  @test projective_scheme_type(UY) == typeof(IP2_UY)
# 
#  @test base_ring_type(IP2_X) == typeof(OO(X))
#  @test base_ring_type(IP2_Y) == typeof(OO(Y))
#  @test base_ring_type(IP2_U) == typeof(OO(U))
#  @test base_ring_type(IP2_UY) == typeof(OO(UY))
# 
#  @test ring_type(IP2_X) == typeof(homogeneous_coordinate_ring(IP2_X))
#  @test ring_type(IP2_Y) == typeof(homogeneous_coordinate_ring(IP2_Y))
#  @test ring_type(IP2_U) == typeof(homogeneous_coordinate_ring(IP2_U))
#  @test ring_type(IP2_UY) == typeof(homogeneous_coordinate_ring(IP2_UY))

  CX, _ = affine_cone(IP2_X)
  CY, _ = affine_cone(IP2_Y)
  CU, _ = affine_cone(IP2_U)
  CUY, _ = affine_cone(IP2_UY)

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

  IP2_Xh = subscheme(IP2_X, gens(homogeneous_coordinate_ring(IP2_X))[1])
  ProjectiveSchemeMor(IP2_Xh, IP2_X, gens(homogeneous_coordinate_ring(IP2_Xh)))
  IP2_Yh = subscheme(IP2_Y, gens(homogeneous_coordinate_ring(IP2_Y))[1])
  ProjectiveSchemeMor(IP2_Yh, IP2_Y, gens(homogeneous_coordinate_ring(IP2_Yh)))
  IP2_Uh = subscheme(IP2_U, gens(homogeneous_coordinate_ring(IP2_U))[1])
  ProjectiveSchemeMor(IP2_Uh, IP2_U, gens(homogeneous_coordinate_ring(IP2_Uh)))
  IP2_UYh = subscheme(IP2_UY, gens(homogeneous_coordinate_ring(IP2_UY))[1])
  ProjectiveSchemeMor(IP2_UYh, IP2_UY, gens(homogeneous_coordinate_ring(IP2_UYh)))
  IP2_Wh = subscheme(IP2_W, gens(homogeneous_coordinate_ring(IP2_W))[1])
  ProjectiveSchemeMor(IP2_Wh, IP2_W, gens(homogeneous_coordinate_ring(IP2_Wh)))

  incYtoX = inclusion_morphism(Y, X)
  h = hom(homogeneous_coordinate_ring(IP2_X), homogeneous_coordinate_ring(IP2_Y), pullback(incYtoX), gens(homogeneous_coordinate_ring(IP2_Y)))
  YtoX = ProjectiveSchemeMor(IP2_Y, IP2_X, 
                             h, 
                             incYtoX
                            );
  @test YtoX == inclusion_morphism(IP2_Y, IP2_X)
  IP2_Y2, map = fiber_product(incYtoX, IP2_X)
  h2 = hom(homogeneous_coordinate_ring(IP2_Y2), homogeneous_coordinate_ring(IP2_Y), gens(homogeneous_coordinate_ring(IP2_Y)))
  YtoY2 = ProjectiveSchemeMor(IP2_Y, IP2_Y2, h2)
  @test compose(YtoY2, map) == YtoX

  incUtoX = inclusion_morphism(U, X)
  h = hom(homogeneous_coordinate_ring(IP2_X), homogeneous_coordinate_ring(IP2_U), pullback(incUtoX), gens(homogeneous_coordinate_ring(IP2_U)))
  UtoX = ProjectiveSchemeMor(IP2_U, IP2_X, 
                             h, 
                             incUtoX
                            );
  @test UtoX == inclusion_morphism(IP2_U, IP2_X)
  IP2_U2, map = fiber_product(incUtoX, IP2_X)
  h2 = hom(homogeneous_coordinate_ring(IP2_U2), homogeneous_coordinate_ring(IP2_U), gens(homogeneous_coordinate_ring(IP2_U)))
  UtoU2 = ProjectiveSchemeMor(IP2_U, IP2_U2, h2)
  @test compose(UtoU2, map) == UtoX

  incUYtoY = inclusion_morphism(UY, Y)
  h = hom(homogeneous_coordinate_ring(IP2_Y), homogeneous_coordinate_ring(IP2_UY), pullback(incUYtoY), gens(homogeneous_coordinate_ring(IP2_UY)))
  UYtoY = ProjectiveSchemeMor(IP2_UY, IP2_Y, 
                              h, 
                              incUYtoY
                             );
  @test UYtoY == inclusion_morphism(IP2_UY, IP2_Y)
  IP2_UY2, map = fiber_product(incUYtoY, IP2_Y)
  h2 = hom(homogeneous_coordinate_ring(IP2_UY2), homogeneous_coordinate_ring(IP2_UY), gens(homogeneous_coordinate_ring(IP2_UY)))
  UYtoUY2 = ProjectiveSchemeMor(IP2_UY, IP2_UY2, h2)
  @test compose(UYtoUY2, map) == UYtoY

  incUYtoX = inclusion_morphism(UY, X)
  h = hom(homogeneous_coordinate_ring(IP2_X), homogeneous_coordinate_ring(IP2_UY), pullback(incUYtoX), gens(homogeneous_coordinate_ring(IP2_UY)))
  UYtoX = ProjectiveSchemeMor(IP2_UY, IP2_X, 
                              h, 
                              incUYtoX
                             );
  IP2_UY2, map = fiber_product(incUYtoX, IP2_X)
  h2 = hom(homogeneous_coordinate_ring(IP2_UY2), homogeneous_coordinate_ring(IP2_UY), gens(homogeneous_coordinate_ring(IP2_UY)))
  UYtoUY2 = ProjectiveSchemeMor(IP2_UY, IP2_UY2, h2)
  @test compose(UYtoUY2, map) == UYtoX

  @test UYtoX == inclusion_morphism(IP2_UY, IP2_X)
  @test compose(UYtoY, YtoX) == UYtoX

  WW = hypersurface_complement(ambient_scheme(W), [x-y])
  phi = MapFromFunc(OO(W), OO(WW), restriction_map(W, WW))
  IP2_WW, map = fiber_product(phi, IP2_W)
  @test base_scheme(IP2_WW) == WW
  @test !(base_scheme(IP2_WW) === WW)

  IP2_QQ = projective_space(QQ, 2)
  id = identity_map(IP2_QQ)
  @test id == identity_map(IP2_QQ)
end

@testset "properties of projective schemes" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  S, _ = grade(R)
  X = ProjectiveScheme(S)
  I = ideal(S, x^2 - y*z)
  Q, _ = quo(S, I)
  C = ProjectiveScheme(Q)
  @test homogeneous_coordinate_ring(C) === Q
  @test dim(C) == 1
  @test degree(C) == 2
  @test is_smooth(C)
  @test arithmetic_genus(C) == 0

  R, (x,y,z,w) = QQ["x", "y", "z", "w"]
  S, _ = grade(R)
  I = ideal(S, [x^4 + y^4 + z^4 + w^4])
  Q, _ = quo(S, I)
  Y = ProjectiveScheme(Q)
  @test homogeneous_coordinate_ring(Y) === Q
  @test dim(Y) == 2
  @test degree(Y) == 4
  @test is_smooth(Y)
  @test arithmetic_genus(Y) == 1
  @test is_reduced(Y)
  @test is_integral(Y)
  @test is_geometrically_integral(Y)
  @test !is_empty(Y)
  @test is_geometrically_reduced(Y)
end

@testset "simple projective spaces" begin
  # Make sure we don't find complicated qoutient rings if not necessary
  P = projective_space(QQ, 1)
  X = covered_scheme(P)
  @test all(U->(OO(U) isa MPolyRing), affine_charts(X))
end

@testset "equality of projective schemes" begin
  P2 = projective_space(GF(17), 2)
  (x0,x1,x2) = homogeneous_coordinates(P2)
  R = homogeneous_coordinate_ring(P2)
  X1 = subscheme(P2, ideal(R, [x0*x1]))
  X2 = subscheme(P2, ideal(R, [x0^2*x1,x0*x1^2,x0*x1*x2]))
  X3 = subscheme(P2, ideal(R,[x0]))
  @test X1 == X2
  @test issubset(X1, P2)
  @test issubset(X1, X2)
end

@testset "closed embeddings" begin
  IP2 = projective_space(QQ, 2)
  S = homogeneous_coordinate_ring(IP2)
  (x, y, z) = gens(S)
  I = ideal(S, [x^2 + y^2 + z^2])
  X, inc = sub(IP2, I)
  X, inc = sub(IP2, gens(I))

  J = ideal(S, [x+y+z])
  X2, inc2 = sub(IP2, J)
  inc_cov = covered_scheme_morphism(inc)
  inc2_cov = covered_scheme_morphism(inc2)
  j1, j2 = fiber_product(inc_cov, inc2_cov)
  @test pushforward(inc_cov)(image_ideal(j2)) == pushforward(inc2_cov)(image_ideal(j1))
  
  @test X === domain(inc)
  @test IP2 === codomain(inc)
  T = homogeneous_coordinate_ring(X)
  Y, inc_Y = sub(X, T[1]*T[2] - T[3]^2)
  @test domain(inc_Y) === Y
  @test codomain(inc_Y) === X
  @test image_ideal(inc_Y) == ideal(T, T[1]*T[2] - T[3]^2)
  map_on_affine_cones(inc_Y)
  inc_comp = compose(inc_Y, inc)
  @test inc_comp isa Oscar.ProjectiveClosedEmbedding
  phi = hom(homogeneous_coordinate_ring(codomain(inc_comp)), 
            homogeneous_coordinate_ring(domain(inc_comp)),
            pullback(inc_comp).(gens(homogeneous_coordinate_ring(codomain(inc_comp))))
           )
  K = kernel(phi)
  @test K == I + ideal(S, S[1]*S[2] - S[3]^2)
end

@testset "empty charts" begin
  P = projective_space(QQ, 2)
  S = homogeneous_coordinate_ring(P)
  X = subscheme(P, S[2])
  Y = covered_scheme(X)
  @test length(affine_charts(Y)) == 2
end
