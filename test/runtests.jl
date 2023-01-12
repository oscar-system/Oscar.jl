using FTheoryTools
using Test
using Oscar


#############################################################
# 0: Set up parameters for tests
#############################################################

# test base
base = TestBase()

# sections, used to test error messages
sec_f = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(projective_space(NormalToricVariety,3))^4)])
sec_g = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)])
sec_a1 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(projective_space(NormalToricVariety,3)))])
sec_a2 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^2)])
sec_a3 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^3)])
sec_a4 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)])
sec_a6 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)])

# specific Tate models, over arbitrary base
istar0_s_auxiliary_base_ring, (a1pp, a2pp, a3pp, a4pp, a6pp, vp, mp) = QQ["a1pp", "a2pp", "a3pp", "a4pp", "a6pp", "vp", "mp"];
t_istar0_s = GlobalTateModel([a1pp * vp^1, mp * vp^1 + a2pp * vp^2, a3pp * vp^2, mp^2 * vp^2 + a4pp * vp^3, a6pp * vp^4], istar0_s_auxiliary_base_ring, 3);
tate_auxiliary_base_ring, (a1p, a2p, a3p, a4p, a6p, v) = QQ["a1p", "a2p", "a3p", "a4p", "a6p", "v"];
t_i1 = GlobalTateModel([a1p * v^0, a2p * v^0, a3p * v^1, a4p * v^1, a6p * v^1], tate_auxiliary_base_ring, 3);
t_i2_ns = GlobalTateModel([a1p * v^0, a2p * v^0, a3p * v^1, a4p * v^1, a6p * v^2], tate_auxiliary_base_ring, 3);
t_i2_s = GlobalTateModel([a1p * v^0, a2p * v^1, a3p * v^1, a4p * v^1, a6p * v^2], tate_auxiliary_base_ring, 3);
t_i3_ns = GlobalTateModel([a1p * v^0, a2p * v^0, a3p * v^2, a4p * v^2, a6p * v^3], tate_auxiliary_base_ring, 3);
t_i3_s = GlobalTateModel([a1p * v^0, a2p * v^1, a3p * v^1, a4p * v^2, a6p * v^3], tate_auxiliary_base_ring, 3);
t_i4_ns = GlobalTateModel([a1p * v^0, a2p * v^0, a3p * v^2, a4p * v^2, a6p * v^4], tate_auxiliary_base_ring, 3);
t_i4_s = GlobalTateModel([a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^2, a6p * v^4], tate_auxiliary_base_ring, 3);
t_i5_ns = GlobalTateModel([a1p * v^0, a2p * v^0, a3p * v^3, a4p * v^3, a6p * v^5], tate_auxiliary_base_ring, 3);
t_i5_s = GlobalTateModel([a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5], tate_auxiliary_base_ring, 3);
t_i6_ns = GlobalTateModel([a1p * v^0, a2p * v^0, a3p * v^3, a4p * v^3, a6p * v^6], tate_auxiliary_base_ring, 3);
t_i6_s = GlobalTateModel([a1p * v^0, a2p * v^1, a3p * v^3, a4p * v^3, a6p * v^6], tate_auxiliary_base_ring, 3);
t_i7_ns = GlobalTateModel([a1p * v^0, a2p * v^0, a3p * v^4, a4p * v^4, a6p * v^7], tate_auxiliary_base_ring, 3);
t_i7_s = GlobalTateModel([a1p * v^0, a2p * v^1, a3p * v^3, a4p * v^4, a6p * v^7], tate_auxiliary_base_ring, 3);
t_ii = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^1, a4p * v^1, a6p * v^1], tate_auxiliary_base_ring, 3);
t_iii = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^1, a4p * v^1, a6p * v^2], tate_auxiliary_base_ring, 3);
t_iv_ns = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^1, a4p * v^2, a6p * v^2], tate_auxiliary_base_ring, 3);
t_iv_s = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^1, a4p * v^2, a6p * v^3], tate_auxiliary_base_ring, 3);
t_istar0_ns = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^2, a4p * v^2, a6p * v^3], tate_auxiliary_base_ring, 3);
t_istar0_ss = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^2, a4p * v^2, a6p * v^4], tate_auxiliary_base_ring, 3);
t_istar1_ns = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^4], tate_auxiliary_base_ring, 3);
t_istar1_s = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5], tate_auxiliary_base_ring, 3);
t_istar2_ns = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^3, a4p * v^3, a6p * v^5], tate_auxiliary_base_ring, 3);
t_istar2_s = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^3, a4p * v^3, a6p * v^6], tate_auxiliary_base_ring, 3);
t_istar3_ns = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^3, a4p * v^4, a6p * v^6], tate_auxiliary_base_ring, 3);
t_istar3_s = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^3, a4p * v^4, a6p * v^7], tate_auxiliary_base_ring, 3);
t_istar4_ns = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^4, a4p * v^4, a6p * v^7], tate_auxiliary_base_ring, 3);
t_istar4_s = GlobalTateModel([a1p * v^1, a2p * v^1, a3p * v^4, a4p * v^4, a6p * v^8], tate_auxiliary_base_ring, 3);
t_ivstar_ns = GlobalTateModel([a1p * v^1, a2p * v^2, a3p * v^2, a4p * v^3, a6p * v^4], tate_auxiliary_base_ring, 3);
t_ivstar_s = GlobalTateModel([a1p * v^1, a2p * v^2, a3p * v^2, a4p * v^3, a6p * v^5], tate_auxiliary_base_ring, 3);
t_iiistar = GlobalTateModel([a1p * v^1, a2p * v^2, a3p * v^3, a4p * v^3, a6p * v^5], tate_auxiliary_base_ring, 3);
t_iistar = GlobalTateModel([a1p * v^1, a2p * v^2, a3p * v^3, a4p * v^4, a6p * v^5], tate_auxiliary_base_ring, 3);
t_nm = GlobalTateModel([a1p * v^1, a2p * v^2, a3p * v^3, a4p * v^4, a6p * v^6], tate_auxiliary_base_ring, 3);


#############################################################
# 1: Global Weierstrass models over concrete base space
#############################################################

w = GlobalWeierstrassModel(TestBase())
# Base.show(w)

@testset "Attributes of global Weierstrass models over concrete base spaces" begin
    @test parent(weierstrass_section_f(w)) == cox_ring(toric_base_space(w))
    @test parent(weierstrass_section_g(w)) == cox_ring(toric_base_space(w))
    @test parent(weierstrass_polynomial(w)) == cox_ring(toric_ambient_space(w))
    @test parent(discriminant(w)) == cox_ring(toric_base_space(w))
    @test dim(toric_base_space(w)) == 3
    @test dim(toric_ambient_space(w)) == 5
    @test is_smooth(toric_ambient_space(w)) == false
    @test toric_variety(cy_hypersurface(w)) == toric_ambient_space(w)
end

@testset "Error messages in global Weierstrass models over concrete base spaces" begin
    @test_throws ArgumentError GlobalWeierstrassModel(sec_f, sec_g, base)
end


#############################################################
# 2: Global Weierstrass models over generic base space
#############################################################

auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"]
w2 = GlobalWeierstrassModel(f, g, auxiliary_base_ring, 3)

@testset "Attributes of global Weierstrass models over generic base space" begin
    @test parent(weierstrass_section_f(w2)) == cox_ring(toric_base_space(w2))
    @test parent(weierstrass_section_g(w2)) == cox_ring(toric_base_space(w2))
    @test parent(weierstrass_polynomial(w2)) == cox_ring(toric_ambient_space(w2))
    @test parent(discriminant(w2)) == cox_ring(toric_base_space(w2))
    @test dim(toric_base_space(w2)) == 3
    @test dim(toric_ambient_space(w2)) == 5
    @test is_smooth(toric_ambient_space(w2)) == false
    @test toric_variety(cy_hypersurface(w2)) == toric_ambient_space(w2)
    @test length(singular_loci(w2)) == 1
end

@testset "Error messages in global Weierstrass models over generic base space" begin
    @test_throws ArgumentError GlobalWeierstrassModel(f, sec_f, auxiliary_base_ring, 3)
    @test_throws ArgumentError GlobalWeierstrassModel(f, g, auxiliary_base_ring, 0)
    @test_throws ArgumentError GlobalWeierstrassModel(f, g, auxiliary_base_ring, 4)
end


#############################################################
# 3: Global Tate models over concrete base space
#############################################################

t = GlobalTateModel(TestBase())
# Base.show(t)

@testset "Attributes of global Tate models over concrete base space" begin
    @test parent(tate_section_a1(t)) == cox_ring(toric_base_space(t))
    @test parent(tate_section_a2(t)) == cox_ring(toric_base_space(t))
    @test parent(tate_section_a3(t)) == cox_ring(toric_base_space(t))
    @test parent(tate_section_a4(t)) == cox_ring(toric_base_space(t))
    @test parent(tate_section_a6(t)) == cox_ring(toric_base_space(t))
    @test parent(tate_polynomial(t)) == cox_ring(toric_ambient_space(t))
    @test parent(discriminant(t)) == cox_ring(toric_base_space(t))
    @test dim(toric_base_space(t)) == 3
    @test dim(toric_ambient_space(t)) == 5
    @test base_fully_specified(t) == true
    @test base_fully_specified(t) == base_fully_specified(global_weierstrass_model(t))
    @test is_smooth(toric_ambient_space(t)) == false
    @test toric_variety(cy_hypersurface(t)) == toric_ambient_space(t)
end

@testset "Error messages in global Tate models over concrete base space" begin
    @test_throws ArgumentError GlobalTateModel([sec_a2, sec_a3, sec_a4, sec_a6], base)
    @test_throws ArgumentError GlobalTateModel([sec_a1, sec_a2, sec_a3, sec_a4, sec_a6], base)
end


#############################################################
# 4: Global Tate models over generic base space
#############################################################

@testset "Attributes of global Tate models over generic base space" begin
    @test parent(tate_section_a1(t_i5_s)) == cox_ring(toric_base_space(t_i5_s))
    @test parent(tate_section_a2(t_i5_s)) == cox_ring(toric_base_space(t_i5_s))
    @test parent(tate_section_a3(t_i5_s)) == cox_ring(toric_base_space(t_i5_s))
    @test parent(tate_section_a4(t_i5_s)) == cox_ring(toric_base_space(t_i5_s))
    @test parent(tate_section_a6(t_i5_s)) == cox_ring(toric_base_space(t_i5_s))
    @test parent(tate_polynomial(t_i5_s)) == cox_ring(toric_ambient_space(t_i5_s))
    @test parent(discriminant(t_i5_s)) == cox_ring(toric_base_space(t_i5_s))
    @test dim(toric_base_space(t_i5_s)) == 3
    @test dim(toric_ambient_space(t_i5_s)) == 5
    @test base_fully_specified(t_i5_s) == false
    @test base_fully_specified(t_i5_s) == base_fully_specified(global_weierstrass_model(t_i5_s))
    @test is_smooth(toric_ambient_space(t_i5_s)) == false
    @test toric_variety(cy_hypersurface(t_i5_s)) == toric_ambient_space(t_i5_s)
end

@testset "Error messages in global Tate models over generic base space" begin
    @test_throws ArgumentError GlobalTateModel([a1p], tate_auxiliary_base_ring, 3)
    @test_throws ArgumentError GlobalTateModel([a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, sec_a6], tate_auxiliary_base_ring, 3)
    @test_throws ArgumentError GlobalTateModel([a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5], tate_auxiliary_base_ring, -1)
    @test_throws ArgumentError GlobalTateModel([a1p * v^0, a2p * v^1, a3p * v^2, a4p * v^3, a6p * v^5], tate_auxiliary_base_ring, 7)
end

@testset "Singular loci of global Tate models over generic base space" begin
    @test length(singular_loci(t_i1)) == 2
    @test length(singular_loci(t_i2_ns)) == 2
    @test length(singular_loci(t_i2_s)) == 2
    @test length(singular_loci(t_i3_ns)) == 2
    @test length(singular_loci(t_i3_s)) == 2
    @test length(singular_loci(t_i4_ns)) == 2
    @test length(singular_loci(t_i4_s)) == 2
    @test length(singular_loci(t_i5_ns)) == 2
    @test length(singular_loci(t_i5_s)) == 2
    @test length(singular_loci(t_i6_ns)) == 2
    @test length(singular_loci(t_i6_s)) == 2
    @test length(singular_loci(t_i7_ns)) == 2
    @test length(singular_loci(t_i7_s)) == 2
    @test length(singular_loci(t_ii)) == 2
    @test length(singular_loci(t_iii)) == 2
    @test length(singular_loci(t_iv_ns)) == 2
    @test length(singular_loci(t_iv_s)) == 2
    @test length(singular_loci(t_istar0_ns)) == 2
    @test length(singular_loci(t_istar0_ss)) == 2
    @test length(singular_loci(t_istar0_s)) == 2
    @test length(singular_loci(t_istar1_ns)) == 2
    @test length(singular_loci(t_istar1_s)) == 2
    @test length(singular_loci(t_istar2_ns)) == 2
    @test length(singular_loci(t_istar2_s)) == 2
    @test length(singular_loci(t_istar3_ns)) == 2
    @test length(singular_loci(t_istar3_s)) == 2
    @test length(singular_loci(t_istar4_ns)) == 2
    @test length(singular_loci(t_istar4_s)) == 2
    @test length(singular_loci(t_ivstar_ns)) == 2
    @test length(singular_loci(t_ivstar_s)) == 2
    @test length(singular_loci(t_iiistar)) == 2
    @test length(singular_loci(t_iistar)) == 2
    @test length(singular_loci(t_nm)) == 2
    @test singular_loci(t_i1)[1][2:3] == ((0, 0, 1), "I_1")
    @test singular_loci(t_i2_ns)[2][2:3] == ((0, 0, 2), "Non-split I_2")
    @test singular_loci(t_i2_s)[2][2:3] == ((0, 0, 2), "Split I_2")
    @test singular_loci(t_i3_ns)[2][2:3] == ((0, 0, 3), "Non-split I_3")
    @test singular_loci(t_i3_s)[2][2:3] == ((0, 0, 3), "Split I_3")
    @test singular_loci(t_i4_ns)[2][2:3] == ((0, 0, 4), "Non-split I_4")
    @test singular_loci(t_i4_s)[2][2:3] == ((0, 0, 4), "Split I_4")
    @test singular_loci(t_i5_ns)[2][2:3] == ((0, 0, 5), "Non-split I_5")
    @test singular_loci(t_i5_s)[2][2:3] == ((0, 0, 5), "Split I_5")
    @test singular_loci(t_i6_ns)[2][2:3] == ((0, 0, 6), "Non-split I_6")
    @test singular_loci(t_i6_s)[2][2:3] == ((0, 0, 6), "Split I_6")
    @test singular_loci(t_i7_ns)[2][2:3] == ((0, 0, 7), "Non-split I_7")
    @test singular_loci(t_i7_s)[2][2:3] == ((0, 0, 7), "Split I_7")
    @test singular_loci(t_ii)[2][2:3] == ((1, 1, 2), "II")
    @test singular_loci(t_iii)[2][2:3] == ((1, 2, 3), "III")
    @test singular_loci(t_iv_ns)[2][2:3] == ((2, 2, 4), "Non-split IV")
    @test singular_loci(t_iv_s)[2][2:3] == ((2, 2, 4), "Split IV")
    @test singular_loci(t_istar0_ns)[2][2:3] == ((2, 3, 6), "Non-split, semi-split, or split I^*_0")
    @test singular_loci(t_istar0_ss)[2][2:3] == ((2, 3, 6), "Semi-split or split I^*_0")
    @test singular_loci(t_istar0_s)[2][2:3] == ((2, 3, 6), "Semi-split or split I^*_0")
    @test singular_loci(t_istar1_ns)[2][2:3] == ((2, 3, 7), "Non-split I^*_1")
    @test singular_loci(t_istar1_s)[2][2:3] == ((2, 3, 7), "Split I^*_1")
    @test singular_loci(t_istar2_ns)[2][2:3] == ((2, 3, 8), "Non-split I^*_2")
    @test singular_loci(t_istar2_s)[2][2:3] == ((2, 3, 8), "Split I^*_2")
    @test singular_loci(t_istar3_ns)[2][2:3] == ((2, 3, 9), "Non-split I^*_3")
    @test singular_loci(t_istar3_s)[2][2:3] == ((2, 3, 9), "Split I^*_3")
    @test singular_loci(t_istar4_ns)[2][2:3] == ((2, 3, 10), "Non-split I^*_4")
    @test singular_loci(t_istar4_s)[2][2:3] == ((2, 3, 10), "Split I^*_4")
    @test singular_loci(t_ivstar_ns)[2][2:3] == ((3, 4, 8), "Non-split IV^*")
    @test singular_loci(t_ivstar_s)[2][2:3] == ((3, 4, 8), "Split IV^*")
    @test singular_loci(t_iiistar)[2][2:3] == ((3, 5, 9), "III^*")
    @test singular_loci(t_iistar)[2][2:3] == ((4, 5, 10), "II^*")
    @test singular_loci(t_nm)[1][2:3] == ((4, 6, 12), "Non-minimal")
end


