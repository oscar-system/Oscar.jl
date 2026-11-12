@testset "InjectiveModuleHom" begin
    Q = affine_semigroup([1 2; 1 3])
    kQ = monoid_algebra(Q,QQ)
    F = faces(kQ)
    G = Oscar.InjectiveResolutions.ambient_group(Q)
    I1 = Oscar.InjectiveResolutions.IndecInj(F[2], G([1,0])) 
    I2 = Oscar.InjectiveResolutions.IndecInj(F[1], G([0,-2])) 
    I = Oscar.InjectiveResolutions.InjMod(kQ,[I1])
    J = Oscar.InjectiveResolutions.InjMod(kQ,[I2])
    phi = Oscar.InjectiveResolutions.InjectiveModuleHom(I,J,matrix(QQ,1,1,[1]))
    @test !is_zero(phi)
    
    phi_reversed = Oscar.InjectiveResolutions.InjectiveModuleHom(J,I,matrix(QQ,1,1,[1]))
    @test is_zero(phi_reversed)

    KQ = quotient_ring_as_module(ideal(kQ,[]))
    e = Oscar.InjectiveResolutions.InjectiveModuleHom(KQ,I,matrix(QQ,1,1,[1]))
    @test is_zero(e)
    
    I1 = Oscar.InjectiveResolutions.IndecInj(F[2], G([0,0])) 
    I = Oscar.InjectiveResolutions.InjMod(kQ,[I1])
    e = Oscar.InjectiveResolutions.InjectiveModuleHom(KQ,I,matrix(QQ,1,1,[1]))
    @test !is_zero(e)

    phi2 = Oscar.InjectiveResolutions.InjectiveModuleHom(I, J, matrix(QQ, 1, 1, [1]))
    @test !is_zero(Oscar.InjectiveResolutions.compose(e, phi2))
end

@testset "shifted localization membership" begin
    Q = affine_semigroup([1 2; 1 3])
    kQ = monoid_algebra(Q,QQ)
    F = faces(kQ)[2]
    @test Oscar.InjectiveResolutions.in_localization([1,2],Q,F.semigroup)
    @test Oscar.InjectiveResolutions.in_localization([0,1],Q,F.semigroup)
    @test !Oscar.InjectiveResolutions.in_localization([0,-1],Q,F.semigroup)

    G = Oscar.InjectiveResolutions.ambient_group(Q)
    @test Oscar.InjectiveResolutions.in_shifted_localization(G([3,4]),G([3,5]),F,Q)
    @test !Oscar.InjectiveResolutions.in_shifted_localization(G([3,9]),G([3,5]),F,Q)
end