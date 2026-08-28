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

# Test signing and verification
println("Testing signing and verification...")

msg = rand(UInt8, 32)

# Test uov_1p (and its pkc/skc variants)
for uov in (Oscar.uov_1p, Oscar.uov_1p_pkc, Oscar.uov_1p_pkc_skc)
    pk, sk = Oscar.keygen(uov)
    sig = Oscar.UOV.sign(uov, msg, sk)
    @test length(sig) == uov.sig_sz
    @test Oscar.UOV.verify(uov, sig, msg, pk)
end

# Test uov_1s (GF(16) field)
for uov in (Oscar.uov_1s, Oscar.uov_1s_pkc, Oscar.uov_1s_pkc_skc)
    pk, sk = Oscar.keygen(uov)
    sig = Oscar.UOV.sign(uov, msg, sk)
    @test length(sig) == uov.sig_sz
    @test Oscar.UOV.verify(uov, sig, msg, pk)
end

# Test uov_3
for uov in (Oscar.uov_3, Oscar.uov_3_pkc, Oscar.uov_3_pkc_skc)
    pk, sk = Oscar.keygen(uov)
    sig = Oscar.UOV.sign(uov, msg, sk)
    @test length(sig) == uov.sig_sz
    @test Oscar.UOV.verify(uov, sig, msg, pk)
end

# Test uov_5 (largest parameter set)
for uov in (Oscar.uov_5, Oscar.uov_5_pkc, Oscar.uov_5_pkc_skc)
    pk, sk = Oscar.keygen(uov)
    sig = Oscar.UOV.sign(uov, msg, sk)
    @test length(sig) == uov.sig_sz
    @test Oscar.UOV.verify(uov, sig, msg, pk)
end

# Test that verification rejects a tampered signature
pk, sk = Oscar.keygen(Oscar.uov_1p)
sig = Oscar.UOV.sign(Oscar.uov_1p, msg, sk)
sig_bad = copy(sig)
sig_bad[1] = sig_bad[1] ⊻ 0xFF
@test !Oscar.UOV.verify(Oscar.uov_1p, sig_bad, msg, pk)

println("Signing and verification tests passed!")
