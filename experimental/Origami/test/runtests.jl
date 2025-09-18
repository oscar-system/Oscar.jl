@testset "Basic origami operations" begin
  h = @perm (1,2)
  v = @perm (2, 3)
  o = origami(h, v)

  @test degree(o) == max(degree(h), degree(v))

  G = symmetric_group(degree(o))

  @test horizontal_perm(o) == G(h)
  @test vertical_perm(o) == G(v)
  @test stratum(o) == [2]
  @test genus(o) == 2
  @test index_monodromy_group(o) == 1

  h2 = @perm (1,2)(3,4)
  v2 = @perm (2,5,3,6)(4,7)
  o2 = origami(h2, v2)

  @test stratum(o2) == [2, 1, 1]
  @test genus(o2) == 3
  @test index_monodromy_group(o2) == 2

end
