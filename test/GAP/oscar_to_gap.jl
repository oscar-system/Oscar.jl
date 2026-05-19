@testset "single cyclotomics" begin
  # from `QQAbFieldElem`
  F, z = abelian_closure(QQ)
  @test GAP.Obj(z(5) + z(5)^4) == e5 + e5^4

  # not supported conversions
  F, z = quadratic_field(5)
  @test_throws ArgumentError GAP.Obj(z)

  @test_throws ArgumentError GAP.Obj(matrix(F, 1, 1, [z]))
end
