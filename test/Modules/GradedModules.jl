@testset "Graded MOdules 1" begin
  Qx, g = PolynomialRing(QQ, 3)
  R = grade(Qx, [1,2,3])[1]
  F = FreeModule(R, [i^2*R.D[1] for i=1:4])
  s = sub(F, [R(g[1])^2*gen(F, 1), R(g[2])*gen(F, 1)])
  q = quo(s, [R(g[1])^4*gen(F, 1), R(g[2])^2*gen(F, 1)]) 
  @test length(presentation(q)) == 3

  h = hom(F,F, [R(g[1]), R(g[2]), R(g[2]), R(g[3])] .* basis(F)) 

  @test length([h(g) for g = basis(F)]) == dim(F)

  @test !Oscar.ishomogeneous(h)

  c = homogeneous_components(h)
  @test length(c) == 3

  for g = keys(c)
    @test degree(c[g]) == g
  end

  H, phi = hom(F, F)
  
  phi(rand(gens(H)))

  H, mH = hom(q,q )

  free_resolution(H)

  homogeneous_component(H, grading_group(F)[0])
end

@testset "Graded Modules 2" begin

  D = free_abelian_group(2)
  Qx, x, y = PolynomialRing(QQ, :x=>1:2, :y=>1:2)
  R = grade(Qx, [D[1], D[1], D[2], D[2]])[1]

  f = x[1]^3-x[2]^3
  g = y[1]^5+y[2]^5
  h = f+g

  F = FreeModule(R, 1)

  s = sub(F, [derivative(h, i)*F[1] for i=1:4])

  free_resolution(s)

  H, mH = hom(s, quo(F, s))

  free_resolution(H)

  tensor_product(s, s, task = :none)

end
