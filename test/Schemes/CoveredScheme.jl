@testset "Covered schemes 1" begin
  R, (x,y) = PolynomialRing(QQ, ["x", "y"])
  X = subscheme(Spec(R), [x^2+y^2])
  P = projective_space(X, 3)
  S = ambient_ring(P)
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
  Pc = as_covered_scheme(P)
  S = ambient_ring(P)
  (x,y,z,w) = gens(S)
  X = subscheme(P, [x*w - y*z])
  @test dim(Pc)==3
  @test dim(as_covered_scheme(X))==2
  Y = subscheme(P, [x*z,y*z])
  @test dim(as_covered_scheme(Y)) == 2
  C = standard_covering(X)
  D, i, j = simplify(C) # not functional for the moment
  @test all( x->(ngens(ambient_ring(x)) == 3), collect(D)) # should be replaced by == 2 when fixed
  @test_broken transition_graph(Pc[1])
  @test_broken transition_graph(C)
end

@testset "standard_covering" begin
  R, t = PolynomialRing(QQ,["t"])
  T = Oscar.standard_spec(subscheme(Spec(R),t))
  Pt= projective_space(T, 2)
  X = as_covered_scheme(Pt)
  @test dim(X) == 2
end
