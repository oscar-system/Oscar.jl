@testset "uov module structure" begin
  # Test that all parameter sets are defined
  @test isdefined(Oscar, :uov_1p)
  @test isdefined(Oscar, :uov_1s)
  @test isdefined(Oscar, :uov_3)
  @test isdefined(Oscar, :uov_5)

  # Test parameter values
  @test uov_1p.gf == 256
  @test uov_1p.n == 112
  @test uov_1p.m == 44

  @test uov_1s.gf == 16
  @test uov_1s.n == 160
  @test uov_1s.m == 64

  @test uov_3.gf == 256
  @test uov_3.n == 184
  @test uov_3.m == 72

  @test uov_5.gf == 256
  @test uov_5.n == 244
  @test uov_5.m == 96

  # Test all variants exist
  @test length(uov_all) == 12
end

@testset "uov functions exported" begin
  # Test that functions are exported
  @test isdefined(Oscar, :keygen)
  @test isdefined(Oscar, :sign)
  @test isdefined(Oscar, :verify)
  @test isdefined(Oscar, :open)
end

@testset "uov signing and verification" begin
  # Test signing and verification
  msg = rand(UInt8, 32)

  # Test uov_1p (and its pkc/skc variants)
  for uov in (uov_1p, uov_1p_pkc, uov_1p_pkc_skc)
      pk, sk = keygen(uov)
      sig = sign(uov, msg, sk)
      @test length(sig) == uov.sig_sz
      @test verify(uov, sig, msg, pk)
  end

  # Test uov_1s (GF(16) field)
  for uov in (uov_1s, uov_1s_pkc, uov_1s_pkc_skc)
      pk, sk = keygen(uov)
      sig = sign(uov, msg, sk)
      @test length(sig) == uov.sig_sz
      @test verify(uov, sig, msg, pk)
  end

  # Test uov_3
  for uov in (uov_3, uov_3_pkc, uov_3_pkc_skc)
      pk, sk = keygen(uov)
      sig = sign(uov, msg, sk)
      @test length(sig) == uov.sig_sz
      @test verify(uov, sig, msg, pk)
  end

  # Test uov_5 (largest parameter set)
  for uov in (uov_5, uov_5_pkc, uov_5_pkc_skc)
      pk, sk = keygen(uov)
      sig = sign(uov, msg, sk)
      @test length(sig) == uov.sig_sz
      @test verify(uov, sig, msg, pk)
  end

  # Test that verification rejects a tampered signature
  pk, sk = keygen(uov_1p)
  sig = sign(uov_1p, msg, sk)
  sig_bad = copy(sig)
  sig_bad[1] = sig_bad[1] ⊻ 0xFF
  @test !verify(uov_1p, sig_bad, msg, pk)
end

@testset "uov message opening" begin
  # Test opening of a signed message
  uov = uov_1p
  pk, sk = keygen(uov)
  msg = UInt8[1, 2, 3]
  sm = vcat(msg, sign(uov, msg, sk))
  uov_open(uov, sm, pk)
end
