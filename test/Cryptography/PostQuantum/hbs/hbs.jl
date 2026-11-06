@testset "hbs hash functions" begin
  # Test domain separation and the encoding of bits for the hash-based schemes
  seed = UInt8[1, 2, 3, 4]
  a1 = Oscar._hbs_adrs(Oscar._HBS_WOTS_HASH, 0, 1)
  a2 = Oscar._hbs_adrs(Oscar._HBS_WOTS_HASH, 1, 0)
  @test length(a1) == 24
  @test a1 != a2
  x = UInt8[9, 9]
  @test length(Oscar._hbs_thash(12, seed, a1, x)) == 12
  @test Oscar._hbs_thash(12, seed, a1, x) == Oscar._hbs_thash(12, seed, a1, x)
  @test Oscar._hbs_thash(12, seed, a1, x) != Oscar._hbs_thash(12, seed, a2, x)
  @test Oscar._hbs_prf(8, seed, a1) != Oscar._hbs_prf(8, seed, a2)
  @test Oscar._hbs_indices(UInt8[0b10110000, 0b01000000], 3, 5) == [22, 1, 0]
  @test_throws ArgumentError Oscar._hbs_indices(UInt8[0xff], 2, 5)
end
