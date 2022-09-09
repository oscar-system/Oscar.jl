using FTheoryTools
using Test
using Oscar

#############################################################
# 1: Construction and properties of global Weierstrass models
#############################################################

w = GenericGlobalWeierstrassModelOverProjectiveSpace(3)

@testset "Test properties of Weierstrass models" begin
    @test parent(poly_f(w)) == cox_ring(toric_base_space(w))
    @test parent(poly_g(w)) == cox_ring(toric_base_space(w))
end


#############################################################
# 2: Test if errors are raised where we want them
#############################################################

P3 = projective_space(NormalToricVariety,3)
f = sum([rand(Int)*b for b in basis_of_global_sections(anticanonical_bundle(P3)^4)]);
g = sum([rand(Int)*b for b in basis_of_global_sections(anticanonical_bundle(P3)^6)]);
F2 = hirzebruch_surface(2)
g2 = sum([rand(Int)*b for b in basis_of_global_sections(anticanonical_bundle(F2)^6)]);

@testset "Test errors to be raised in Weierstrass models" begin
    @test_throws ArgumentError GenericGlobalWeierstrassModelOverProjectiveSpace(-1)
    @test_throws ArgumentError GenericGlobalWeierstrassModelOverProjectiveSpace(0)
    @test_throws ArgumentError toric_base_space(GlobalWeierstrassModel([f,g]))
    @test_throws ArgumentError GlobalWeierstrassModel([f,g2])
end
