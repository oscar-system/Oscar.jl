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

