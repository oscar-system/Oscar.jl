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

@testset "hbs merkle trees" begin
  # Test building trees, authentication paths and their verification
  leaves = [UInt8[i, 2i] for i in 1:8]
  levels = merkle_tree(leaves)
  root = merkle_root(levels)
  @test length.(levels) == [8, 4, 2, 1]
  @test levels[2][1] == shake_256(vcat(leaves[1], leaves[2]), 2)
  for i in 0:7
    path = merkle_auth_path(levels, i)
    @test length(path) == 3
    @test merkle_verify(root, leaves[i + 1], i, path)
    @test !merkle_verify(root, leaves[i + 1], (i + 1) % 8, path)
    @test !merkle_verify(root, UInt8[0, 0], i, path)
  end
  @test length(merkle_tree([UInt8[1]])) == 1
  @test_throws ArgumentError merkle_tree([UInt8[i] for i in 1:3])
  @test_throws ArgumentError merkle_auth_path(levels, 8)

  # a custom node function is used for building and verifying
  node = (h, i, l, r) -> shake_256(vcat(UInt8[h, i], l, r), 4)
  levels = merkle_tree([UInt8[i] for i in 1:4]; node)
  path = merkle_auth_path(levels, 2)
  @test merkle_verify(merkle_root(levels), UInt8[3], 2, path; node)
  @test !merkle_verify(merkle_root(levels), UInt8[3], 2, path)
end

@testset "hbs lamport signatures" begin
  # Test signing and verification with Lamport's one-time scheme
  rng = Random.MersenneTwister(1)
  for n in (1, 4, 16)
    pk, sk = lamport_keygen(n; rng)
    msg = Vector{UInt8}("oscar")
    sig = lamport_sign(n, msg, sk)
    @test length(sig) == 8n * n
    @test lamport_verify(n, sig, msg, pk)
    @test !lamport_verify(n, sig, Vector{UInt8}("oscaR"), pk)
    bad = copy(sig)
    bad[1] ⊻= 0x01
    @test !lamport_verify(n, bad, msg, pk)
  end
  @test_throws ArgumentError lamport_keygen(0)
  @test_throws ArgumentError lamport_verify(4, UInt8[], UInt8[], UInt8[])
end
