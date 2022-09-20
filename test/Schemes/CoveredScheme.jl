@testset "Covered schemes 1" begin
  R, (x,y) = PolynomialRing(QQ, ["x", "y"])
  X = subscheme(Spec(R), [x^2+y^2])
  P = projective_space(X, 3)
  S = homogeneous_poly_ring(P)
  (u, v) = gens(S)[1], gens(S)[2]
  h = u^3 
  h = u^3 + u^2
  h = u^3 + (u^2)*v 
  h = u^3 + u^2*v - OO(X)(x)*v^3
  Z = subscheme(P, h)
  C = standard_covering(Z)
  f = dehomogenize(Z, 1)
  @test f(u) == gens(OO(C[2]))[1]
end

@testset "Covered schemes 2" begin
  P = projective_space(QQ, ["x", "y", "z", "w"])
  S = homogeneous_poly_ring(P)
  (x,y,z,w) = gens(S)
  X = subscheme(P, [x*w - y*z])
  C = standard_covering(X)
  D, i, j = simplify(C) # not functional for the moment
  @test all( x->(ngens(ambient_ring(x)) == 3), collect(D)) # should be replaced by == 2 when fixed
end
