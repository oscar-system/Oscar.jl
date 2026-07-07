@testset "singular fibers of families of hypersurfaces" begin
  R, t = QQ[:t]
  P, (x, y) = R[:x, :y]
  f = (x-1)^2 + (y-1)^2 - t

  PP, (x, y, t) = QQ[:x, :y, :t]
  fl = hom(P, PP, hom(R, PP, t), [x, y])
  ff = fl(f)
  I = ideal(PP, push!([derivative(ff, i) for i in 1:ngens(P)], ff))
  d0 = only(gens(eliminate(I, [x, y])))
  d0 = evaluate(d0, [zero(R), zero(R), gen(R)])
  d1 = Oscar.discriminant(f)
  @test is_associated(d0, d1)
  
  supps = Oscar._support_sets([f, derivative(f, 1), derivative(f, 2)])
  @test !is_zero(Oscar.a_resultant(supps))
end

