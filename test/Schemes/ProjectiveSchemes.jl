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

