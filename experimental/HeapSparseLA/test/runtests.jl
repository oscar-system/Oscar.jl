@testset "heap based sparse linear algebra" begin
  n = Oscar.HeapNode(1, ZZ(5))
  n.left = Oscar.HeapNode(2, ZZ(7))
  n.right = Oscar.HeapNode(2, ZZ(-3))
  n.right.right = Oscar.HeapNode(3, ZZ(8))
  v = Oscar.HeapSRow(ZZ, n)
  w = Oscar.HeapSRow(ZZ, [1 => ZZ(5), 2 => ZZ(4), 3 => ZZ(11), 3=>ZZ(-3)])
  @test !isdefined(v, :index_cache)
  @test v[3] == 8
  @test isdefined(v, :index_cache) && v.index_cache[3] == 8
  @test Oscar.tree_size(v.top) == 4
  @test v[2] == 4
  @test v.index_cache[2] == 4
  @test Oscar.tree_size(v.top) == 3

  u = v + w
  @test u == 2*v
end

