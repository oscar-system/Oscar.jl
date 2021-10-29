@testset "integer-localizations" begin
  S = FmpzPowersOfElement([5,14])
  @test 25 in S
  @test !(33 in S)
  @test 5*5*7 in S
  @test ambient_ring(S) == ZZ
  W = localize_at(S)
  @test original_ring(W) == ambient_ring(S)
  @test inverted_set(W) == S
  a = W(3)
  b = W(7, 5)
  @test a+b == W(3*5+7, 5)
  @test a-b == W(3*5-7, 5)
  @test 1//b == W(5,7)
  @test b^2 == W(49, 25)
  c = W(4, 5)
  @test c//c == one(W)

  U = FmpzComplementOfPrimeIdeal(13)
  @test !(13^5 in U)
  @test !(13*4289729837 in U)
  @test 5783790198374098 in U
  @test ambient_ring(U) == ZZ
  W = localize_at(U)
  @test original_ring(W) == ambient_ring(U)
  @test inverted_set(W) == U
  a = W(4, 17)
  b = W(4*17)
  b = b//W(19)
  @test a//b == W(19//17^2)
  @test a - b == W( 4//17 - 4*17//19 )
  @test a + b == W( 4//17 + 4*17//19 )
  @test a * b == W( 4//17 * 4*17//19 )

  O = FmpzComplementOfZeroIdeal()
  @test 234890 in O
  @test !(0 in O)
  @test ambient_ring(O) == ZZ
  W = localize_at(O)
  @test original_ring(W) == ambient_ring(O)
  @test inverted_set(W) == O
  a = W(4, 17)
  b = W(4*17)
  b = b//W(19)
  @test a//b == W(19//17^2)
  @test a - b == W( 4//17 - 4*17//19 )
  @test a + b == W( 4//17 + 4*17//19 )
  @test a * b == W( 4//17 * 4*17//19 )
end
