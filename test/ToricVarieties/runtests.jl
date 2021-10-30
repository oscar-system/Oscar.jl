using Oscar
using Test

C = Oscar.positive_hull([1 1; -1 1])
antv = AffineNormalToricVariety(C)

@testset "Affine toric varieties" begin
    @test issmooth(antv) == false
    @test isorbifold(antv) == true
    @test length(affine_open_covering(antv)) == 1
    # @test toric_ideal_binomial_generators(antv) == [-1 -1 2]
    @test rank(torusinvariant_divisor_group(antv)) == 2
    @test rank(character_lattice(antv)) == 2
    map = map_from_character_to_principal_divisors(antv)
    @test rank(domain(map)) == 2
    @test rank(codomain(map)) == 2
    @test elementary_divisors(codomain(map_from_weil_divisors_to_class_group(antv))) == [ 2 ]
    @test elementary_divisors(class_group(antv)) == [ 2 ]
    @test ngens(cox_ring(antv)) == 2
    @test length(torusinvariant_prime_divisors(antv)) == 2
end

ntv = NormalToricVariety(C)
associated_affine_variety = AffineNormalToricVariety(ntv)

@testset "Affine toric varieties created as general normal toric varieties" begin
    @test isaffine(ntv) == true
    @test length(affine_open_covering(associated_affine_variety)) == 1
end

square = Oscar.cube(2)
nf = Oscar.normal_fan(square)
ntv2 = NormalToricVariety(nf)
ntv3 = NormalToricVariety(square)

@testset "Toric varieties from polyhedral fans" begin
    @test iscomplete(ntv2) == true
    @test iscomplete(ntv3) == true
end

P2 = toric_projective_space(2)

@testset "Projective space" begin
    @test isnormal(P2) == true
    @test isaffine(P2) == false
    @test isprojective(P2) == true
    @test issmooth(P2) == true
    @test iscomplete(P2) == true
    @test hastorusfactor(P2) == false
    @test isorbifold(P2) == true
    @test issimplicial(P2) == true
    @test ith_betti_number(P2, 0) == 1
    @test ith_betti_number(P2, 1) == 0
    @test ith_betti_number(P2, 2) == 1
    @test ith_betti_number(P2, 3) == 0
    @test ith_betti_number(P2, 4) == 1
    S = cox_ring(P2)
    @test ngens(S) == 3
    @test length(stanley_reisner_ideal(P2).gens) == 1
    @test length(irrelevant_ideal(P2).gens) == 3
end

H5 = hirzebruch_surface(5)

@testset "Hirzebruch surface" begin
    @test isnormal(H5) == true
    @test isaffine(H5) == false
    @test isprojective(H5) == true
    @test issmooth(H5) == true
    @test iscomplete(H5) == true
    @test hastorusfactor(H5) == false
    @test isorbifold(H5) == true
    @test issimplicial(H5) == true
    @test isgorenstein(H5) == true
    @test isq_gorenstein(H5) == true
    @test isfano(H5) == false
    nef_cone(H5)
    mori_cone(H5)
    @test dim(H5) == 2
    @test dim_of_torusfactor(H5) == 0
    @test euler_characteristic(H5) == 4
    @test ith_betti_number(H5, 0) == 1
    @test ith_betti_number(H5, 1) == 0
    @test ith_betti_number(H5, 2) == 2
    @test ith_betti_number(H5, 3) == 0
    @test ith_betti_number(H5, 4) == 1
    @test length(affine_open_covering(H5)) == 4
    @test fan_of_variety(H5).pm_fan.FAN_DIM == 2
    @test fan(H5).pm_fan.FAN_DIM == 2
    @test rank(torusinvariant_divisor_group(H5)) == 4
    @test rank(character_lattice(H5)) == 2
    map = map_from_character_to_principal_divisors(H5)
    @test rank(domain(map)) == 2
    @test rank(codomain(map)) == 4
    @test rank(class_group(H5)) == 2
    @test rank(codomain(map_from_weil_divisors_to_class_group(H5))) == 2
    @test ngens(cox_ring(H5)) == 4
    @test length(stanley_reisner_ideal(H5).gens) == 2
    @test length(irrelevant_ideal(H5).gens) == 4
end

@testset "delPezzo surfaces" begin
    @test_throws ArgumentError del_pezzo(-1)
    dP0 = del_pezzo(0)
    @test length(torusinvariant_prime_divisors(dP0)) == 3
    dP1 = del_pezzo(1)
    @test length(torusinvariant_prime_divisors(dP1)) == 4
    dP2 = del_pezzo(2)
    @test length(torusinvariant_prime_divisors(dP2)) == 5
    dP3 = del_pezzo(3)
    @test length(torusinvariant_prime_divisors(dP3)) == 6
    @test_throws ArgumentError del_pezzo(4)
end

blowup_variety = blowup_on_ith_minimal_torus_orbit(P2, 1)

@testset "Blowup of projective space" begin
    @test isnormal(blowup_variety) == true
    @test isaffine(blowup_variety) == false
    @test isprojective(blowup_variety) == true
    @test issmooth(blowup_variety) == true
    @test iscomplete(blowup_variety) == true
    @test hastorusfactor(blowup_variety) == false
    @test isorbifold(blowup_variety) == true
    @test issimplicial(blowup_variety) == true
    @test ith_betti_number(blowup_variety, 0) == 1
    @test ith_betti_number(blowup_variety, 1) == 0
    @test ith_betti_number(blowup_variety, 2) == 2
    @test ith_betti_number(blowup_variety, 3) == 0
    @test ith_betti_number(blowup_variety, 4) == 1
    @test euler_characteristic(blowup_variety) == 4
end

v = H5 * P2

@testset "Direct products" begin
    @test isnormal(v) == true
    @test isaffine(v) == false
    @test isprojective(v) == true
    @test issmooth(v) == true
    @test iscomplete(v) == true
    @test hastorusfactor(v) == false
    @test isorbifold(v) == true
    @test issimplicial(v) == true
    @test ith_betti_number(v, 0) == 1
    @test ith_betti_number(v, 1) == 0
    @test ith_betti_number(v, 2) == 3
    @test ith_betti_number(v, 3) == 0
    @test ith_betti_number(v, 4) == 4
    @test ith_betti_number(v, 5) == 0
    @test ith_betti_number(v, 6) == 3
    @test ith_betti_number(v, 7) == 0
    @test ith_betti_number(v, 8) == 1
end

@testset "ComparisonWithProjectiveSpace" begin
    @test isprojective_space(H5) == false
    @test isprojective_space(P2) == true
    @test isprojective_space(blowup_variety) == false
    @test isprojective_space(v) == false
    @test isprojective_space(ntv) == false
    @test isprojective_space(ntv2) == false
end

D=ToricDivisor(H5, [0,0,0,0])
D2 = DivisorOfCharacter(H5, [1,2])

@testset "Divisors" begin
    @test isprime_divisor(D) == false
    @test iscartier(D) == true
    @test isprincipal(D) == true
    @test isbasepoint_free(D) == true
    @test isample(D) == false
    @test isvery_ample(D) == false
    @test isnef(D) == true
    @test isintegral(D) == true
    @test isq_cartier(D) == true
    @test coefficients(D) == [0,0,0,0]
    @test isprime_divisor(D2) == false
    @test iscartier(D2) == true
    @test isprincipal(D2) == true
    @test isbasepoint_free(D2) == true
    @test isample(D2) == false
    @test isvery_ample(D2) == false
    @test isnef(D2) == true
    @test isintegral(D2) == true
    @test isq_cartier(D2) == true
    @test isprime_divisor(D2) == false
    @test coefficients(D2) == [1, 2, 9, -2]
end

p = polyhedron(D)

@testset "Polytopes of divisors" begin
    @test dim(p) == 0
    @test ambient_dim(p) == 2
end
