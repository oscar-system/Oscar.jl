using FTheoryTools
using Test
using Oscar


#############################################################
# 1: Construction and properties of global Weierstrass models
#############################################################

w = GenericGlobalWeierstrassModelOverProjectiveSpace(3)
Base.show(w)

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
    @test_throws ArgumentError GlobalWeierstrassModel([f,g,g2])
end


#############################################################
# 3: Construction and properties of global Tate models
#############################################################

test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")
t = GenericGlobalTateModel(base)
Base.show(t)

@testset "Test properties of global Tate models" begin
    @test parent(a1(t)) == cox_ring(toric_ambient_space(t))
    @test parent(a2(t)) == cox_ring(toric_ambient_space(t))
    @test parent(a3(t)) == cox_ring(toric_ambient_space(t))
    @test parent(a4(t)) == cox_ring(toric_ambient_space(t))
    @test parent(a6(t)) == cox_ring(toric_ambient_space(t))
    @test parent(pt(t)) == cox_ring(toric_ambient_space(t))
    @test dim(toric_base_space(t)) == 3
    @test dim(toric_ambient_space(t)) == 5
    @test is_smooth(toric_ambient_space(t)) == false
end

@testset "Test errors to be raised in Tate models" begin
    @test_throws ArgumentError GenericGlobalTateModel(hirzebruch_surface(1))
end
