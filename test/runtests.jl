using FTheoryTools
using Test
using Oscar


#############################################################
# 1: Compute base space for explicit fibrations
#############################################################

test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")


#############################################################
# 2: Test global Weierstrass models
#############################################################

w = GenericGlobalWeierstrassModel(base)
Base.show(w)

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
# 3: Test global Tate models
#############################################################

t = GenericGlobalTateModel(base)
Base.show(t)

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

@testset "Error messages in global Tate models" begin
    @test_throws ArgumentError GenericGlobalTateModel(hirzebruch_surface(1))
end


#############################################################
# 4: Test global Tate models over generic base space
#############################################################

auxiliary_base_ring, (a10,a21,a32,a43,a65,w) = QQ["a10", "a21", "a32", "a43", "a65", "w"]
a1 = a10
a2 = a21 * w
a3 = a32 * w^2
a4 = a43 * w^3
a6 = a65 * w^5
ais = [a1, a2, a3, a4, a6]
t2 = GlobalTateModel(ais, auxiliary_base_ring)

@testset "Attributes of global Tate models over generic base spaces" begin
    @test parent(tate_section_a1(t2)) == cox_ring(auxiliary_ambient_space(t2))
    @test parent(tate_section_a2(t2)) == cox_ring(auxiliary_ambient_space(t2))
    @test parent(tate_section_a3(t2)) == cox_ring(auxiliary_ambient_space(t2))
    @test parent(tate_section_a4(t2)) == cox_ring(auxiliary_ambient_space(t2))
    @test parent(tate_section_a6(t2)) == cox_ring(auxiliary_ambient_space(t2))
    @test parent(tate_polynomial(t2)) == cox_ring(auxiliary_ambient_space(t2))
    @test dim(auxiliary_base_space(t2)) == 6
    @test dim(auxiliary_ambient_space(t2)) == 8
    @test is_smooth(auxiliary_ambient_space(t2)) == false
end

@testset "Error messages in global Tate models over generic base space" begin
    @test_throws ArgumentError GlobalTateModel([a1], auxiliary_base_ring)
end
