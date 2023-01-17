@testset "pullbacks of Cartier divisors" begin 
  P = projective_space(QQ, 2)
  SP = ambient_coordinate_ring(P)
  (x, y, z) = gens(SP)

  A = [1 2 4; 1 3 9; 1 5 25]
  v = A*[x, y, z]
  psi = hom(SP, SP, v)
  g = ProjectiveSchemeMor(P, P, psi)
  g_cov = covered_scheme_morphism(g)

  X = covered_scheme(P)

  D = IdDict{AbsSpec, RingElem}()
  H = x^2 + y^2 + z^2
  for U in affine_charts(X)
    D[U] = dehomogenize(P, U)(H)
  end

  C = CartierDivisor(X, D)
  D = pullback(g_cov)(C)
  for U in patches(trivializing_covering(D))
    @test length(D(U)) == 1
  end
end
