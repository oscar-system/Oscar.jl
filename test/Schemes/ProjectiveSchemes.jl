@testset "projective_schemes_1" begin

# test for relative projective space over a polynomial ring
R, (x,y) = QQ["x", "y"]
R_ext, _ = PolynomialRing(R, ["u", "v"])
S, (u,v) = grade(R_ext, [1,1])

I = ideal(S, [x*v - y*u])
X = ProjectiveScheme(S, I)
CX = affine_cone(X)
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
