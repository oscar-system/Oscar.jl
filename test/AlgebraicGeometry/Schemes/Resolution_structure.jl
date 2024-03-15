@testset "desingularization of curves" begin
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