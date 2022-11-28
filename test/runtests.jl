using FTheoryTools
using Test
using Oscar


#############################################################
# 1: Global Weierstrass models over concrete base space
#############################################################

w = GlobalWeierstrassModel(TestBase())
Base.show(w)

@testset "Attributes of global Weierstrass models" begin
    @test parent(weierstrass_section_f(w)) == cox_ring(toric_ambient_space(w))
    @test parent(weierstrass_section_g(w)) == cox_ring(toric_ambient_space(w))
    @test parent(weierstrass_polynomial(w)) == cox_ring(toric_ambient_space(w))
    @test parent(discriminant(w)) == cox_ring(toric_ambient_space(w))
    @test dim(toric_base_space(w)) == 3
    @test dim(toric_ambient_space(w)) == 5
    @test is_smooth(toric_ambient_space(w)) == false
    @test toric_variety(cy_hypersurface(w)) == toric_ambient_space(w)
end


#############################################################
# 2: Global Weierstrass models over generic base space
#############################################################

auxiliary_base_ring, (f, g) = QQ["f", "g"]
auxiliary_ring2, (x, y) = QQ["x", "y"]
w2 = GlobalWeierstrassModel(f, g, auxiliary_base_ring)

@testset "Attributes of global Weierstrass models over generic base spaces" begin
    @test parent(weierstrass_section_f(w2)) == cox_ring(toric_ambient_space(w2))
    @test parent(weierstrass_section_g(w2)) == cox_ring(toric_ambient_space(w2))
    @test parent(weierstrass_polynomial(w2)) == cox_ring(toric_ambient_space(w2))
    @test parent(discriminant(w2)) == cox_ring(toric_ambient_space(w2))
    @test dim(toric_base_space(w2)) == 2
    @test dim(toric_ambient_space(w2)) == 4
    @test is_smooth(toric_ambient_space(w2)) == false
    @test toric_variety(cy_hypersurface(w2)) == toric_ambient_space(w2)
end

@testset "Error messages in global Weierstrass models over generic base space" begin
    @test_throws ArgumentError GlobalWeierstrassModel(f, x^2, auxiliary_base_ring)
end


#############################################################
# 3: Global Tate models over concrete base space
#############################################################

t = GlobalTateModel(TestBase())
Base.show(t)

@testset "Attributes of global Tate models" begin
    @test parent(tate_section_a1(t)) == cox_ring(toric_ambient_space(t))
    @test parent(tate_section_a2(t)) == cox_ring(toric_ambient_space(t))
    @test parent(tate_section_a3(t)) == cox_ring(toric_ambient_space(t))
    @test parent(tate_section_a4(t)) == cox_ring(toric_ambient_space(t))
    @test parent(tate_section_a6(t)) == cox_ring(toric_ambient_space(t))
    @test parent(tate_polynomial(t)) == cox_ring(toric_ambient_space(t))
    @test parent(discriminant(t)) == cox_ring(toric_ambient_space(t))
    @test dim(toric_base_space(t)) == 3
    @test dim(toric_ambient_space(t)) == 5
    @test base_fully_specified(t) == true
    @test base_fully_specified(t) == base_fully_specified(global_weierstrass_model(t))
    @test is_smooth(toric_ambient_space(t)) == false
    @test toric_variety(cy_hypersurface(t)) == toric_ambient_space(t)
end


#############################################################
# 4: Global Tate models over generic base space
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
    @test parent(tate_section_a1(t2)) == cox_ring(toric_ambient_space(t2))
    @test parent(tate_section_a2(t2)) == cox_ring(toric_ambient_space(t2))
    @test parent(tate_section_a3(t2)) == cox_ring(toric_ambient_space(t2))
    @test parent(tate_section_a4(t2)) == cox_ring(toric_ambient_space(t2))
    @test parent(tate_section_a6(t2)) == cox_ring(toric_ambient_space(t2))
    @test parent(tate_polynomial(t2)) == cox_ring(toric_ambient_space(t2))
    @test parent(discriminant(t2)) == cox_ring(toric_ambient_space(t2))
    @test dim(toric_base_space(t2)) == 6
    @test dim(toric_ambient_space(t2)) == 8
    @test base_fully_specified(t2) == false
    @test base_fully_specified(t2) == base_fully_specified(global_weierstrass_model(t2))
    @test is_smooth(toric_ambient_space(t2)) == false
    @test toric_variety(cy_hypersurface(t2)) == toric_ambient_space(t2)
end

@testset "Error messages in global Tate models over generic base space" begin
    @test_throws ArgumentError GlobalTateModel([a1], auxiliary_base_ring)
end
