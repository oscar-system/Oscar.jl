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
  I = ideal(R,[(x+y)^3+z^4])
  IS = IdealSheaf(W,affine_charts(W)[1], gens(I))
  IY = Oscar.max_order_locus(IS)
  decomp = Oscar.maximal_associated_points(IY)
  @test length(decomp) == 1
  @test decomp[1] == IdealSheaf(W,affine_charts(W)[1],[z,x+y])
  JS = IdealSheaf(W, affine_charts(W)[1],[x^2*y^2-z^5,x^3])
  li = Oscar._delta_list(JS)
  @test length(li) == 3
  @test li[1] == JS
  @test Oscar.radical(li[2]) == IdealSheaf(W,affine_charts(W)[1],[x,z])
  @test Oscar.radical(li[2]) == Oscar.radical(Oscar.locus_of_order_geq_b(JS,2))
  @test Oscar.radical(li[3]) == IdealSheaf(W,affine_charts(W)[1],[x,y,z])
end

