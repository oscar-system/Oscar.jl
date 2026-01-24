@testset "Basic origami operations" begin
  h = @perm (1, 2)
  v = @perm (2, 3)
  o = origami(h, v)

  @test degree(o) == max(degree(h), degree(v))
  @test (@inferred degree(o)) == 3

  G = symmetric_group(degree(o))

  @test horizontal_perm(o) == G(h)
  @test vertical_perm(o) == G(v)
  @test stratum(o) == [2]
  @test genus(o) == 2
  @test index_monodromy_group(o) == 1

  h2 = @perm (1, 2)(3, 4)
  v2 = @perm (2, 5, 3, 6)(4, 7)
  o2 = origami(h2, v2)

  @test (@inferred stratum(o2)) == [2, 1, 1]
  @test (@inferred genus(o2)) == 3
  @test (@inferred index_monodromy_group(o2)) == 2
  @test (@inferred sum_of_lyapunov_exponents(o2)) == 11//6

  h3 = @perm (1, 2, 3, 4, 5, 6)
  v3 = @perm (1, 6, 5, 4, 3, 2)
  o3 = origami(h3, v3)

  @test issetequal(
    translations(o3),
    [one(perm_group(o3)), cperm([1, 2, 3, 4, 5, 6]),
      cperm([1, 3, 5], [2, 4, 6]), cperm([1, 4], [2, 5], [3, 6]),
      cperm([1, 5, 3], [2, 6, 4]),
      cperm([1, 6, 5, 4, 3, 2])],
  )
  @test (@inferred is_hyperelliptic(o3)) == true

  h4 = @perm (1, 2, 3, 4, 5, 6)
  v4 = @perm (1, 4)(2, 6, 5, 3)
  o4 = origami(h4, v4)
  @test (@inferred translations(o4)) ==
        [one(perm_group(o4)), cperm([1, 4], [2, 5], [3, 6])] ||
    (@inferred translations(
    o4
  )) == [
    cperm([1, 4], [2, 5], [3, 6]), one(perm_group(o4))
  ]
  @test (@inferred is_hyperelliptic(o4)) == false
end
