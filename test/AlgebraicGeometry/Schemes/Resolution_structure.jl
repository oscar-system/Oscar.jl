@testset "embedded desingularization of curves" begin
  R,(x,y) = polynomial_ring(QQ,2)
  I=ideal(R,[x^2-y^5])
  W = CoveredScheme(Spec(R))
  IS = IdealSheaf(W,affine_charts(W)[1],gens(I))
  inc_X = Oscar.CoveredClosedEmbedding(W,IS)
  phi = Oscar.embedded_desingularization(inc_X)
  @test length(phi.ex_div) == 4
  @test is_empty(singular_locus(domain(phi.embeddings[end]))[1])
  @test is_one(ideal_sheaf(phi.ex_div[1]) + ideal_sheaf(phi.ex_div[2]) + image_ideal(phi.embeddings[end]))
  @test is_one(ideal_sheaf(phi.ex_div[1]) + ideal_sheaf(phi.ex_div[3]) + image_ideal(phi.embeddings[end]))
  @test is_one(ideal_sheaf(phi.ex_div[1]) + ideal_sheaf(phi.ex_div[4]) + image_ideal(phi.embeddings[end]))
  @test is_one(ideal_sheaf(phi.ex_div[2]) + ideal_sheaf(phi.ex_div[3]) + image_ideal(phi.embeddings[end]))
  @test is_one(ideal_sheaf(phi.ex_div[2]) + ideal_sheaf(phi.ex_div[4]) + image_ideal(phi.embeddings[end]))
  @test is_one(ideal_sheaf(phi.ex_div[3]) + ideal_sheaf(phi.ex_div[4]) + image_ideal(phi.embeddings[end]))
  @test !is_empty(singular_locus(domain(phi.embeddings[2]))[1])
  @test is_empty(singular_locus(domain(phi.embeddings[3]))[1])
end

@testset "non-embedded desingularization of curves" begin
  R,(x,y) = polynomial_ring(QQ,2)
  I=ideal(R,[x^2-y^5])
  W = CoveredScheme(Spec(R))
  IS = IdealSheaf(W,affine_charts(W)[1],gens(I))
  X = subscheme(IS)
  phi = Oscar.desingularization(X)
  @test length(phi.ex_div) == 2
  aff_charts = affine_charts(domain(phi))
  sl1,_ = singular_locus(aff_charts[1])
  sl2,_ = singular_locus(aff_charts[2])
  sl3,_ = singular_locus(aff_charts[3])
  @test is_one(modulus(OO(sl1)))
  @test is_one(modulus(OO(sl2)))
  @test is_one(modulus(OO(sl3)))
end

@testset "order of an ideal" begin
  R,(x,y,z) = polynomial_ring(QQ,3)
  W = CoveredScheme(spec(R))
  I = ideal(R,[x*y,x^3+y^3+z^3-x*y*z])
  IS = IdealSheaf(W,affine_charts(W)[1], gens(I))
  Y = Oscar.max_order_locus(IS)
end

