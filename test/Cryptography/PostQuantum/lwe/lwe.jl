@testset "lwe parameters" begin
  # Test that parameter sets are defined and validated
  @test lwe_toy.n == 16
  @test lwe_standard.n == 512
  @test length(lwe_all) == 3
  @test_throws ArgumentError instantiate_lwe(16, 32, 100, 2, "bad")
  @test_throws ArgumentError instantiate_lwe(0, 32, 1031, 2, "bad")
end

@testset "lwe encryption and decryption" begin
  # Test that decryption inverts encryption for all parameter sets
  for p in lwe_all
    pk, sk = lwe_keygen(p)
    for msg in (UInt8[], UInt8[0x00], UInt8[0xff], rand(UInt8, 16))
      ct = lwe_encrypt(p, pk, msg)
      @test lwe_decrypt(p, sk, ct) == msg
    end
  end
end

@testset "lwe key properties" begin
  # Test that the public key has small noise and that a wrong key fails
  p = lwe_toy
  (A, b), s = lwe_keygen(p)
  noise = mod.(b - A * s, p.q)
  @test all(x -> min(x, p.q - x) <= p.eta, noise)
  msg = rand(UInt8, 32)
  pk2, sk2 = lwe_keygen(p)
  ct = lwe_encrypt(p, (A, b), msg)
  @test lwe_decrypt(p, sk2, ct) != msg
end
