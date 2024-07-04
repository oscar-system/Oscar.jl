@testset "common refinements" begin
  IP2 = projective_space(QQ, [:x, :y, :z])
  S = homogeneous_coordinate_ring(IP2)
  X = covered_scheme(IP2)
  (x, y, z) = gens(S)
  Phi=hom(S,S,[x+y+z,3*x-5*y+2*z,7*x+y-3*z])
  phi = ProjectiveSchemeMor(IP2, IP2, Phi)
  f = covered_scheme_morphism(phi)
  cov1 = domain(covering_morphism(f))

  Psi=hom(S,S,[7*x+3*y+z,x-5*y+2*z,-2x+3*y-3*z])
  psi = ProjectiveSchemeMor(IP2, IP2, Psi)
  g = covered_scheme_morphism(psi)
  cov2 = domain(covering_morphism(g))

  common_refinement([cov1, cov2], default_covering(X))

  h1 = compose(f, f)
  h2 = compose(f, g)
  h3 = compose(g, f)
  h4 = compose(g, g)

  common_refinement(domain_covering.([h1, h2, f]), default_covering(X))
end
