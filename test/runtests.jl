using FTheoryTools
using Test
using Oscar


#############################################################
# 1: Compute fibrations
#############################################################

test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")
t = GenericGlobalTateModel(base)
Base.show(t)
w = GenericGlobalWeierstrassModel(base)
Base.show(w)

#############################################################
# 1: Test global Weierstrass models
#############################################################

@testset "Attributes of global Weierstrass models" begin
    @test parent(weierstrass_section_f(w)) == cox_ring(toric_ambient_space(w))
    @test parent(weierstrass_section_g(w)) == cox_ring(toric_ambient_space(w))
    @test parent(weierstrass_polynomial(w)) == cox_ring(toric_ambient_space(w))
    @test dim(toric_base_space(w)) == 3
    @test dim(toric_ambient_space(w)) == 5
    @test is_smooth(toric_ambient_space(w)) == false
end

@testset "Error messages in global Weierstrass models" begin
    @test_throws ArgumentError GenericGlobalWeierstrassModel(hirzebruch_surface(1))
end


#############################################################
# 2: Test global Tate models
#############################################################

@testset "Attributes of global Tate models" begin
    @test parent(tate_section_a1(t)) == cox_ring(toric_ambient_space(t))
    @test parent(tate_section_a2(t)) == cox_ring(toric_ambient_space(t))
    @test parent(tate_section_a3(t)) == cox_ring(toric_ambient_space(t))
    @test parent(tate_section_a4(t)) == cox_ring(toric_ambient_space(t))
    @test parent(tate_section_a6(t)) == cox_ring(toric_ambient_space(t))
    @test parent(tate_polynomial(t)) == cox_ring(toric_ambient_space(t))
    @test dim(toric_base_space(t)) == 3
    @test dim(toric_ambient_space(t)) == 5
    @test is_smooth(toric_ambient_space(t)) == false
end

@testset "Error messages in global Weierstrass models" begin
    @test_throws ArgumentError GenericGlobalTateModel(hirzebruch_surface(1))
end
