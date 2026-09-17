@testset "keccak shake_256" begin
  # Test against known SHAKE256 test vectors, including messages of
  # length 135, 136 and 272 at the block boundaries of the 136-byte rate
  @test bytes2hex(Oscar.shake_256(UInt8[], 64)) == "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762fd75dc4ddd8c0f200cb05019d67b592f6fc821c49479ab48640292eacb3b7c4be"
  @test bytes2hex(Oscar.shake_256([0x61, 0x62, 0x63], 32)) == "483366601360a8771c6863080cc4114d8db44530f8f1e1ee4f94ea37e78b5739"
  @test bytes2hex(Oscar.shake_256(collect(UInt8, 0:134), 32)) == "c45dae624ad8a2f5aa7bac9d7557737fd91c96eedb70a6be5574d57a844eade0"
  @test bytes2hex(Oscar.shake_256(collect(UInt8, 0:135), 32)) == "b7ff4073b3f5a8eabd6e17705ca7f6761a31058f9df781a6a47e3a3063b9d67a"
  @test bytes2hex(Oscar.shake_256(collect(UInt8, 0:136), 32)) == "01d90952c642a5eb2a8fc9d713f843a45d7ac05132dddcb2efc9bebc27e37bcb"
  @test bytes2hex(Oscar.shake_256(UInt8[(7i + 3) % 256 for i in 0:271], 32)) == "fbb7df100461f5db3224c3b715603b5a52f92bd0f761d9361d4aa0613d47033b"
end
