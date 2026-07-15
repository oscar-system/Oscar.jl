@testset "single cyclotomics" begin
  # from `QQAbFieldElem`
  F, z = abelian_closure(QQ)
  e5 = GAP.Globals.E(5)
  
  @test GAP.Obj(z(5) + z(5)^4) == e5 + e5^4
end


@testset "GapGroup and GapGroupElem" begin
  # `GapGroup` to GAP group, Perm
  G = symmetric_group(5)
  val = GAP.Globals.SymmetricGroup(5)
  @test GAP.Obj(G) == val

  # `GapGroupElem` to GAP group element, Perm
  g = perm(G, [2,3,1,5,4])
  val = GAP.evalstr("(1,2,3)(4,5)")
  @test GAP.Obj(g) == val
end
