@testset "embedded desingularization of curves" begin
  R,(x,y) = polynomial_ring(QQ,2)
  I=ideal(R,[x^2-y^5])
  W = AffineScheme(R)
  IS = IdealSheaf(W,I)
  WC = scheme(IS)
  inc_X = Oscar.CoveredClosedEmbedding(WC,IS)
  phi = embedded_desingularization(inc_X)
  H = IdealSheaf(W, ideal(R, x + 3*y), covered_scheme=WC)
  @test is_subset(total_transform(phi, H), strict_transform(phi, H))
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
  @test length(components(exceptional_divisor(phi))) == 4
  @test length(components(exceptional_locus(phi))) == 4
end

@testset "non-embedded desingularization of curves" begin
  R,(x,y) = polynomial_ring(QQ,2)
  I=ideal(R,[x^2-y^5])
  W = AffineScheme(R)
  IS = IdealSheaf(W,I)
  X = subscheme(IS)
  U = first(affine_charts(X))
  phi = desingularization_only_blowups(X)
  H = IdealSheaf(U, ideal(OO(U), x + 3*y), covered_scheme=X)
  @test is_subset(total_transform(phi, H), strict_transform(phi, H))
  @test dim(center(Oscar.last_map(phi))) == 0
  @test length(phi.ex_div) == 2
  aff_charts = affine_charts(domain(phi))
  sl1,_ = singular_locus(aff_charts[1])
  sl2,_ = singular_locus(aff_charts[2])
  sl3,_ = singular_locus(aff_charts[3])
  @test is_one(modulus(OO(sl1)))
  @test is_one(modulus(OO(sl2)))
  @test is_one(modulus(OO(sl3)))
  @test length(components(exceptional_divisor(phi))) == 1
  @test length(components(exceptional_locus(phi))) == 1
  psi = desingularization(X)
  @test length(morphisms(psi)) == 1
  @test morphism(psi,1) isa NormalizationMorphism
end

@testset "non-embedded desingularization Lipman (dim=2)" begin
  R,(x,y,z) = polynomial_ring(QQ,3)
  I=ideal(R,[(x^2-y^5)*(x^2+y^2+z^4)])
  W = AffineScheme(R)
  IS = IdealSheaf(W,I)
  X = subscheme(IS)
  U = first(affine_charts(X))
  phi = desingularization(X)
  @test morphism(phi,1) isa NormalizationMorphism
  @test morphism(phi,2) isa BlowupMorphism
  @test morphism(phi,3) isa BlowupMorphism
  aff_charts = affine_charts(domain(phi))
  sl1,_ = singular_locus(aff_charts[1])
  @test is_one(modulus(OO(sl1)))
  sl5,_ = singular_locus(aff_charts[5])
  @test is_one(modulus(OO(sl5)))
  sl6,_ = singular_locus(aff_charts[6])
  @test is_one(modulus(OO(sl6)))
  @test length(components(exceptional_divisor(phi))) == 2
  @test_broken length(components(exceptional_locus(phi))) == 2
end

@testset "order of an ideal" begin
  R,(x,y,z) = polynomial_ring(QQ,3)
  W = AffineScheme(R)
  I = ideal(R,[(x+y)^3+z^4])
  IS = IdealSheaf(W,I)
  WC = scheme(IS)
  IY = locus_of_maximal_order(IS)
  decomp = Oscar.maximal_associated_points(IY)
  @test length(decomp) == 1
  @test decomp[1] == IdealSheaf(W,ideal(R,[z,x+y]), covered_scheme=WC)
  JS = IdealSheaf(W, ideal(R,[x^2*y^2-z^5,x^3]), covered_scheme=WC)
  li = Oscar._delta_list(JS)
  @test length(li) == 3
  @test li[1] == JS
  @test radical(li[2]) == IdealSheaf(W,ideal(R,[x,z]), covered_scheme=WC)
  @test radical(li[2]) == radical(Oscar.locus_of_order_geq_b(JS,2))
  @test radical(li[3]) == IdealSheaf(W,ideal(R,[x,y,z]),covered_scheme=WC)
end

