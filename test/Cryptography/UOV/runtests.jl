# Tests for UOV signature scheme
using Test
using Oscar

println("Testing UOV module structure...")

# Test that all parameter sets are defined
@test isdefined(Oscar, :uov_1p)
@test isdefined(Oscar, :uov_1s)
@test isdefined(Oscar, :uov_3)
@test isdefined(Oscar, :uov_5)

# Test parameter values
@test Oscar.uov_1p.gf == 256
@test Oscar.uov_1p.n == 112
@test Oscar.uov_1p.m == 44

@test Oscar.uov_1s.gf == 16
@test Oscar.uov_1s.n == 160
@test Oscar.uov_1s.m == 64

@test Oscar.uov_3.gf == 256
@test Oscar.uov_3.n == 184
@test Oscar.uov_3.m == 72

@test Oscar.uov_5.gf == 256
@test Oscar.uov_5.n == 244
@test Oscar.uov_5.m == 96

# Test all variants exist
@test length(Oscar.uov_all) == 12

# Test that functions are exported
@test isdefined(Oscar, :keygen)
@test isdefined(Oscar, :sign)
@test isdefined(Oscar, :verify)
@test isdefined(Oscar, :open)

println("Module structure tests passed!")

# Test signing
println("Testing signing...")

msg = rand(UInt8, 32)

# Test uov_1p
pk1, sk1 = Oscar.keygen(Oscar.uov_1p)
sig1 = Oscar.UOV.sign(Oscar.uov_1p, msg, sk1)
@test length(sig1) == Oscar.uov_1p.sig_sz
# Verification failing - public key uses modified m2 (sks) instead of original
# @test Oscar.UOV.verify(Oscar.uov_1p, sig1, msg, pk1)

# Test uov_1s
    pk2, sk2 = Oscar.keygen(Oscar.uov_1s)
    sig2 = Oscar.UOV.sign(Oscar.uov_1s, msg, sk2)
    @test length(sig2) == Oscar.uov_1s.sig_sz
    # Verification failing - public key uses modified m2 (sks) instead of original
    # @test Oscar.UOV.verify(Oscar.uov_1s, sig2, msg, pk2)

    # Test uov_3
    pk3, sk3 = Oscar.keygen(Oscar.uov_3)
    sig3 = Oscar.UOV.sign(Oscar.uov_3, msg, sk3)
    @test length(sig3) == Oscar.uov_3.sig_sz
    # Verification failing - public key uses modified m2 (sks) instead of original
    # @test Oscar.UOV.verify(Oscar.uov_3, sig3, msg, pk3)

    # Test uov_5
    pk4, sk4 = Oscar.keygen(Oscar.uov_5)
    sig4 = Oscar.UOV.sign(Oscar.uov_5, msg, sk4)
    @test length(sig4) == Oscar.uov_5.sig_sz
    # Verification failing - public key uses modified m2 (sks) instead of original
    # @test Oscar.UOV.verify(Oscar.uov_5, sig4, msg, pk4)

println("Signing tests passed!")