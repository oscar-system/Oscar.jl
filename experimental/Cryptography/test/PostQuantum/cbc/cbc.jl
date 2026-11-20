using Random
using Test
using Oscar

@testset "cbc goppa codes" begin
  # Test that Patterson's algorithm corrects up to t errors of a Goppa code
  rng = Random.Xoshiro(1)
  for (m, n, t) in ((4, 16, 2), (5, 32, 3), (6, 64, 5))
    C = Oscar._goppa_code(rng, m, n, t)
    @test size(C.H) == (m * t, n)
    @test size(C.G) == (n - m * t, n)
    @test is_zero(C.G * transpose(C.H))
    for w in 0:t
      pos = sort(randperm(rng, n)[1:w])
      e = zeros(Int, n)
      e[pos] .= 1
      s = [is_one(c) ? 1 : 0 for c in vec(collect(C.H * matrix(GF(2), n, 1, e)))]
      @test Oscar._goppa_decode(C, s) == pos
    end
  end
end

@testset "mceliece parameters" begin
  # Test that parameter sets are defined and validated
  @test mceliece_toy.k == 34
  @test mceliece_light.k == 524
  @test mceliece_standard.k == 2720
  @test length(mceliece_all) == 3
  @test_throws ArgumentError instantiate_mceliece(6, 100, 5, "bad")
  @test_throws ArgumentError instantiate_mceliece(6, 64, 1, "bad")
end

@testset "mceliece encryption and decryption" begin
  # Test that decryption inverts encryption for all parameter sets
  for p in mceliece_all
    pk, sk = keygen(p)
    @test size(pk) == (p.k, p.n)
    for msg in (UInt8[], UInt8[0xff], rand(UInt8, 16))
      ct = encrypt(p, pk, msg)
      @test ncols(ct) == p.n
      @test decrypt(p, sk, ct) == msg
    end
  end
end

@testset "mceliece errors and keys" begin
  # Test that ciphertexts carry t errors and that a wrong key fails
  p = mceliece_toy
  rng = Random.Xoshiro(2)
  pk, sk = keygen(p; rng)
  msg = rand(rng, UInt8, 12)
  ct = encrypt(p, pk, msg; rng)
  bits = Oscar._cbc_pad(Oscar._bytes_to_bits(msg), p.k)
  X = matrix(GF(2), length(bits) ÷ p.k, p.k, bits)
  @test all(i -> count(is_one, ct[i, :] - (X * pk)[i, :]) == p.t, 1:nrows(ct))
  pk2, sk2 = keygen(p; rng)
  @test decrypt(p, sk2, ct) != msg
end

@testset "niederreiter parameters" begin
  # Test that parameter sets are defined and validated
  @test niederreiter_toy.l == 22
  @test niederreiter_light.l == ndigits(binomial(big(1024), 50); base = 2) - 1
  @test length(niederreiter_all) == 3
  @test_throws ArgumentError instantiate_niederreiter(6, 64, 1, "bad")
  @test_throws ArgumentError instantiate_niederreiter(17, 64, 5, "bad")
end

@testset "niederreiter constant weight encoding" begin
  # Test that numbers and vectors of weight t correspond bijectively
  n, t = 10, 3
  words = [Oscar._constant_weight_encode(big(N), n, t) for N in 0:binomial(n, t) - 1]
  @test allunique(words)
  @test all(w -> length(w) == t && issorted(w) && 1 <= w[1] && w[end] <= n, words)
  @test all(N -> Oscar._constant_weight_decode(words[N + 1]) == N, 0:binomial(n, t) - 1)
  N = binomial(big(1024), 50) - 1
  @test Oscar._constant_weight_decode(Oscar._constant_weight_encode(N, 1024, 50)) == N
end

@testset "niederreiter encryption and decryption" begin
  # Test that decryption inverts encryption for all parameter sets
  for p in niederreiter_all
    pk, sk = keygen(p)
    @test size(pk) == (p.m * p.t, p.n)
    for msg in (UInt8[], UInt8[0xff], rand(UInt8, 16))
      ct = encrypt(p, pk, msg)
      @test ncols(ct) == p.m * p.t
      @test decrypt(p, sk, ct) == msg
    end
  end
end

@testset "niederreiter keys" begin
  # Test that encryption is deterministic and that a wrong key fails
  p = niederreiter_toy
  rng = Random.Xoshiro(3)
  pk, sk = keygen(p; rng)
  msg = rand(rng, UInt8, 12)
  ct = encrypt(p, pk, msg)
  @test ct == encrypt(p, pk, msg)
  pk2, sk2 = keygen(p; rng)
  @test decrypt(p, sk2, ct) != msg
end

@testset "hqc reed-solomon and reed-muller codes" begin
  # Test that the codes of HQC correct errors up to their capacities
  rng = Random.Xoshiro(4)
  for (n1, k) in ((10, 4), (46, 16))
    msg = rand(rng, UInt8, k)
    c = Oscar._rs_encode(msg, n1)
    @test length(c) == n1
    for w in 0:(n1 - k) ÷ 2
      r = copy(c)
      for i in randperm(rng, n1)[1:w]
        r[i] ⊻= rand(rng, 0x01:0xff)
      end
      @test Oscar._rs_decode(r, k) == msg
    end
    r = copy(c)
    for i in 1:n1 - k
      r[i] ⊻= 0x01
    end
    @test Oscar._rs_decode(r, k) != msg
  end
  for b in (0x00, 0x01, 0xa5, 0xff), mult in (1, 3)
    word = Oscar._rm_encode(b, mult)
    @test length(word) == 128 * mult
    @test Oscar._rm_decode(word, mult) == b
    word[randperm(rng, 128 * mult)[1:31 * mult]] .⊻= 1
    @test Oscar._rm_decode(word, mult) == b
  end
end

@testset "hqc parameters" begin
  # Test that parameter sets are defined and validated
  @test hqc_1.n == 17669
  @test hqc_1.n2 == 384
  @test hqc_5.k == 32
  @test length(hqc_all) == 4
  @test_throws ArgumentError instantiate_hqc(1279, 10, 4, 1, 8, 8, 8, "bad")
  @test_throws ArgumentError instantiate_hqc(1290, 10, 4, 1, 8, 8, 8, "bad")
  @test_throws ArgumentError instantiate_hqc(1283, 10, 10, 1, 8, 8, 8, "bad")
end

@testset "hqc encryption and decryption" begin
  # Test that decryption inverts encryption for all parameter sets
  for p in hqc_all
    pk, sk = keygen(p)
    for msg in (UInt8[], UInt8[0xff], rand(UInt8, 3 * p.k))
      ct = encrypt(p, pk, msg)
      @test length(ct) == cld(8 * length(msg) + 1, 8 * p.k)
      @test all(((u, v),) -> degree(u) < p.n && degree(v) < p.n1 * p.n2, ct)
      @test decrypt(p, sk, ct) == msg
    end
  end
end

@testset "hqc keys" begin
  # Test the relation between the keys and that a wrong key fails
  p = hqc_toy
  rng = Random.Xoshiro(5)
  (h, s), (x, y) = keygen(p; rng)
  _, x_ = polynomial_ring(GF(2), :x)
  @test s == mod(x + h * y, x_^p.n - 1)
  @test all(z -> count(!is_zero, coefficients(z)) == p.w, (x, y))
  msg = rand(rng, UInt8, 64)
  ct = encrypt(p, (h, s), msg; rng)
  @test decrypt(p, (x, y), ct) == msg
  pk2, sk2 = keygen(p; rng)
  @test decrypt(p, sk2, ct) != msg
end
