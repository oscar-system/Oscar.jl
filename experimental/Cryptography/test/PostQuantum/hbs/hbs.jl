using Random
using Test
using Oscar

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

@testset "hbs winternitz signatures" begin
  # Test the digits, signing and recovering the public key of WOTS+
  @test Oscar._wots_lens(16, 16) == (32, 3)
  @test Oscar._wots_lens(8, 4) == (32, 4)
  @test Oscar._wots_lens(16, 256) == (16, 2)
  @test Oscar._wots_digits(1, 16, UInt8[0xa5]) == [10, 5, 0, 15]
  for (n, w) in ((8, 4), (16, 16), (16, 256), (5, 8))
    sk_seed = rand(UInt8, n)
    seed = rand(UInt8, n)
    loc = (1, 2, 3)
    msg = rand(UInt8, n)
    d = Oscar._wots_digits(n, w, msg)
    @test length(d) == Oscar._wots_len(n, w)
    @test all(0 .<= d .< w)
    pk = Oscar._wots_pk(n, w, sk_seed, seed, loc)
    sig = Oscar._wots_sign(n, w, msg, sk_seed, seed, loc)
    @test Oscar._wots_pk_from_sig(n, w, sig, msg, seed, loc) == pk
    msg2 = copy(msg)
    msg2[1] ⊻= 0x01
    @test Oscar._wots_pk_from_sig(n, w, sig, msg2, seed, loc) != pk
    @test Oscar._wots_pk_from_sig(n, w, sig, msg, seed, (1, 2, 4)) != pk
  end
end

@testset "hbs parameter sets" begin
  # Test the validation of parameters and the sizes of the predefined sets
  @test length(hbs_all) == 7
  @test hbs_mss_toy.scheme == :mss && hbs_fors_toy.scheme == :hypertree
  @test hbs_mss_toy.sig_sz == 8 + 8 + 36 * 8 + 4 * 8
  @test hbs_fors_standard.hp == 3 && hbs_fors_standard.len == 35
  @test hbs_fors_standard.sig_sz == 16 + 33 * 7 * 16 + 20 * 38 * 16
  @test hbs_horst_toy.sig_sz == 8 + (8 * 5 * 8 + 4 * 8) + 2 * (36 * 8 + 3 * 8)
  @test_throws ArgumentError instantiate_hbs(0, 16, 4, "bad")
  @test_throws ArgumentError instantiate_hbs(16, 12, 4, "bad")
  @test_throws ArgumentError instantiate_hbs(16, 16, 4, "bad"; d = 2)
  @test_throws ArgumentError instantiate_hbs(16, 16, 5, "bad"; d = 2, fts = :fors, tau = 4, k = 4)
  @test_throws ArgumentError instantiate_hbs(16, 16, 4, "bad"; fts = :horst, tau = 4, k = 4, x = 5)
  @test_throws ArgumentError instantiate_hbs(16, 16, 4, "bad"; fts = :foo, tau = 4, k = 4)
  @test sprint(show, hbs_mss_toy) == "Hash-based signature parameters hbs-mss-toy (mss, n = 8, w = 4, h = 4)"
end

@testset "hbs fors signatures" begin
  # Test signing and recovering the public key of FORS, also with a wrong digest
  for (n, a, k) in ((8, 4, 8), (16, 6, 5), (8, 3, 7))
    sk_seed = rand(UInt8, n)
    seed = rand(UInt8, n)
    loc = (0, 5, 2)
    md = rand(UInt8, Oscar._fors_md_size(a, k))
    pk = Oscar._fors_pk(n, a, k, sk_seed, seed, loc)
    sig = Oscar._fors_sign(n, a, k, md, sk_seed, seed, loc)
    @test length(sig) == k * (a + 1) * n
    @test Oscar._fors_pk_from_sig(n, a, k, sig, md, seed, loc) == pk
    md2 = copy(md)
    md2[1] ⊻= 0x80
    @test Oscar._fors_pk_from_sig(n, a, k, sig, md2, seed, loc) != pk
    sig2 = copy(sig)
    sig2[end] ⊻= 0x01
    @test Oscar._fors_pk_from_sig(n, a, k, sig2, md, seed, loc) != pk
    @test Oscar._fors_pk_from_sig(n, a, k, sig, md, seed, (0, 5, 3)) != pk
  end
end

@testset "hbs horst signatures" begin
  # Test signing and recovering the public key of HORST, also with modified signatures
  for (n, tau, k, x) in ((8, 6, 8, 2), (16, 5, 6, 0), (8, 4, 5, 4))
    sk_seed = rand(UInt8, n)
    seed = rand(UInt8, n)
    loc = (0, 1, 7)
    md = rand(UInt8, Oscar._horst_md_size(tau, k))
    pk = Oscar._horst_pk(n, tau, sk_seed, seed, loc)
    sig = Oscar._horst_sign(n, tau, k, x, md, sk_seed, seed, loc)
    @test length(sig) == k * (tau - x + 1) * n + 2^x * n
    @test Oscar._horst_pk_from_sig(n, tau, k, x, sig, md, seed, loc) == pk
    md2 = copy(md)
    md2[1] ⊻= 0x80
    @test Oscar._horst_pk_from_sig(n, tau, k, x, sig, md2, seed, loc) != pk
    sig2 = copy(sig)
    sig2[1] ⊻= 0x01
    @test Oscar._horst_pk_from_sig(n, tau, k, x, sig2, md, seed, loc) != pk
    sig3 = copy(sig)
    sig3[end] ⊻= 0x01
    @test Oscar._horst_pk_from_sig(n, tau, k, x, sig3, md, seed, loc) != pk
  end
end

@testset "hbs trees over winternitz keys" begin
  # Test signing with a leaf of a tree and recovering the root from the signature
  n, w, hp = 8, 4, 3
  sk_seed = rand(UInt8, n)
  seed = rand(UInt8, n)
  levels = Oscar._xmss_tree(n, w, hp, sk_seed, seed, 2, 5)
  root = merkle_root(levels)
  m = rand(UInt8, n)
  for leaf in (0, 3, 7)
    sig = Oscar._xmss_sign(n, w, levels, m, sk_seed, seed, 2, 5, leaf)
    @test length(sig) == (Oscar._wots_len(n, w) + hp) * n
    @test Oscar._xmss_root_from_sig(n, w, hp, sig, m, seed, 2, 5, leaf) == root
    @test Oscar._xmss_root_from_sig(n, w, hp, sig, m, seed, 2, 5, (leaf + 1) % 8) != root
    @test Oscar._xmss_root_from_sig(n, w, hp, sig, m, seed, 2, 6, leaf) != root
  end
end

@testset "hbs single Merkle tree scheme" begin
  # Test the stateful scheme: all signatures of a key verify, the counter runs out
  p = hbs_mss_toy
  pk, sk = Oscar._mss_keygen(p, Random.MersenneTwister(3))
  @test length(pk) == p.pk_sz && length(sk) == p.sk_sz
  sigs = Vector{Vector{UInt8}}()
  for i in 1:2^p.h
    msg = UInt8[i, 7]
    sig = Oscar._mss_sign(p, msg, sk)
    @test length(sig) == p.sig_sz
    @test Oscar._mss_verify(p, sig, msg, pk)
    @test !Oscar._mss_verify(p, sig, UInt8[i, 8], pk)
    push!(sigs, sig)
  end
  @test length(unique(sig[1:8] for sig in sigs)) == 2^p.h
  @test_throws ArgumentError Oscar._mss_sign(p, UInt8[1], sk)
  bad = copy(sigs[1])
  bad[end] ⊻= 0x01
  @test !Oscar._mss_verify(p, bad, UInt8[1, 7], pk)
  bad = copy(sigs[1])
  bad[1:8] = Oscar._hbs_adrs(2^p.h)
  @test !Oscar._mss_verify(p, bad, UInt8[1, 7], pk)
end

@testset "hbs hypertree scheme" begin
  # Test signing and verification with HORST and FORS in a hypertree of several trees
  for p in (hbs_horst_toy, hbs_fors_toy)
    pk, sk = Oscar._ht_keygen(p, Random.MersenneTwister(5))
    @test length(pk) == p.pk_sz && length(sk) == p.sk_sz
    for msg in (UInt8[], UInt8[1, 2, 3], rand(UInt8, 100))
      sig = Oscar._ht_sign(p, msg, sk)
      @test length(sig) == p.sig_sz
      @test Oscar._ht_verify(p, sig, msg, pk)
      @test !Oscar._ht_verify(p, sig, vcat(msg, UInt8[0]), pk)
      for pos in (1, p.n + 1, p.sig_sz)
        bad = copy(sig)
        bad[pos] ⊻= 0x01
        @test !Oscar._ht_verify(p, bad, msg, pk)
      end
    end
    # signing is deterministic and does not change the key
    @test Oscar._ht_sign(p, UInt8[1], sk) == Oscar._ht_sign(p, UInt8[1], sk)
  end
end

@testset "hbs interface" begin
  # Test the public interface for all parameter sets that are fast enough
  for p in (hbs_mss_toy, hbs_horst_toy, hbs_fors_toy, hbs_horst_light, hbs_fors_light)
    pk, sk = keygen(p; rng = Random.MersenneTwister(11))
    msg = Vector{UInt8}("hash-based signatures")
    sig = sign(p, msg, sk)
    @test length(sig) == p.sig_sz
    @test verify(p, sig, msg, pk)
    @test !verify(p, sig, Vector{UInt8}("other message"), pk)
    @test open(p, vcat(msg, sig), pk) == msg
    @test open(p, vcat(UInt8[0], sig), pk) === nothing
    @test open(p, sig[1:end-1], pk) === nothing
    @test_throws ArgumentError verify(p, sig[1:end-1], msg, pk)
  end
  # keys from different seeds do not accept each other's signatures
  pk1, sk1 = keygen(hbs_fors_toy)
  pk2, _ = keygen(hbs_fors_toy)
  @test !verify(hbs_fors_toy, sign(hbs_fors_toy, UInt8[1], sk1), UInt8[1], pk2)
  # the scheme with a single tree is stateful, the one with a hypertree is not
  pk, sk = keygen(hbs_mss_toy)
  sk0 = copy(sk)
  sign(hbs_mss_toy, UInt8[1], sk)
  @test sk != sk0
  pk, sk = keygen(hbs_fors_toy)
  sk0 = copy(sk)
  sign(hbs_fors_toy, UInt8[1], sk)
  @test sk == sk0
end

@testset "hbs standard parameters" begin
  # Test one signature with the largest parameter set
  p = hbs_fors_standard
  pk, sk = keygen(p)
  sig = sign(p, UInt8[1, 2, 3], sk)
  @test length(sig) == p.sig_sz
  @test verify(p, sig, UInt8[1, 2, 3], pk)
end
