using Oscar
using Test

C = Oscar.positive_hull([1 1; -1 1])
antv = AffineNormalToricVariety(C)
f = map_from_character_lattice_to_torusinvariant_weil_divisor_group(antv)

@testset "Affine toric varieties" begin
    @test issmooth(antv) == false
    @test isorbifold(antv) == true
    @test dim(fan( antv )) == 2
    @test dim(cone(antv)) == 2
    @test length(affine_open_covering(antv)) == 1
    @test length(gens(toric_ideal(antv))) == 1
    @test rank(torusinvariant_weil_divisor_group(antv)) == 2
    @test rank(character_lattice(antv)) == 2
    @test rank(domain(f)) == 2
    @test rank(codomain(f)) == 2
    @test elementary_divisors(codomain(map_from_torusinvariant_weil_divisor_group_to_class_group(antv))) == [ 2 ]
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

r = [1 0 0; 1 0 1; 1 1 1; 1 1 0]
c = IncidenceMatrix([[1,2,3,4]])
antv2 = NormalToricVariety(PolyhedralFan(r, c))

@testset "Affine toric varieties created from fan" begin
    @test istrivial(picard_group(antv2)) == true
end

cyc = CyclicQuotientSingularity(2,1)

@testset "Cyclic quotient singularities" begin
    @test isaffine(cyc) == true
    @test issmooth( cyc ) == false
    @test issimplicial( cyc ) == true
    @test continued_fraction_hirzebruch_jung( cyc )[1] == 2
    @test dual_continued_fraction_hirzebruch_jung(cyc)[1] == 2
end

square = Oscar.cube(2)
nf = Oscar.normal_fan(square)
ntv2 = NormalToricVariety(nf)
ntv3 = NormalToricVariety(square)

@testset "Toric varieties from polyhedral fans" begin
    @test iscomplete(ntv2) == true
    @test iscomplete(ntv3) == true
    @test rank(torusinvariant_cartier_divisor_group(ntv2)) == 4
    @test rank(domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(ntv2))) == 4
end

P2 = NormalToricVariety(normal_fan(Oscar.simplex(2)))

@testset "Projective space" begin
    @test isnormal(P2) == true
    @test isaffine(P2) == false
    @test isprojective(P2) == true
    @test issmooth(P2) == true
    @test iscomplete(P2) == true
    @test hastorusfactor(P2) == false
    @test isorbifold(P2) == true
    @test issimplicial(P2) == true
    @test betti_number(P2, 0) == 1
    @test betti_number(P2, 1) == 0
    @test betti_number(P2, 2) == 1
    @test betti_number(P2, 3) == 0
    @test betti_number(P2, 4) == 1
    S = cox_ring(P2)
    @test ngens(S) == 3
    @test length(stanley_reisner_ideal(P2).gens) == 1
    @test length(irrelevant_ideal(P2).gens) == 3
end

fan_rays = [1 0; 0 1; -1 5; 0 -1]
fan_cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,1]])
H5 = NormalToricVariety(PolyhedralFan(fan_rays, fan_cones))

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
    @test is_q_gorenstein(H5) == true
    @test isfano(H5) == false
    nef_cone(H5)
    mori_cone(H5)
    @test dim(H5) == 2
    @test dim_of_torusfactor(H5) == 0
    @test euler_characteristic(H5) == 4
    @test betti_number(H5, 0) == 1
    @test betti_number(H5, 1) == 0
    @test betti_number(H5, 2) == 2
    @test betti_number(H5, 3) == 0
    @test betti_number(H5, 4) == 1
    @test length(affine_open_covering(H5)) == 4
    @test fan(H5).pm_fan.FAN_DIM == 2
    @test rank(torusinvariant_weil_divisor_group(H5)) == 4
    @test rank(character_lattice(H5)) == 2
    map = map_from_character_lattice_to_torusinvariant_weil_divisor_group(H5)
    @test rank(domain(map)) == 2
    @test rank(codomain(map)) == 4
    @test rank(class_group(H5)) == 2
    @test rank(codomain(map_from_torusinvariant_weil_divisor_group_to_class_group(H5))) == 2
    @test ngens(cox_ring(H5)) == 4
    @test length(stanley_reisner_ideal(H5).gens) == 2
    @test length(irrelevant_ideal(H5).gens) == 4
end

dP0 = NormalToricVariety(normal_fan(Oscar.simplex(2)))

fan_rays = [1 0; 0 1; -1 0; -1 -1]
fan_cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,1]])
dP1 = NormalToricVariety(PolyhedralFan(fan_rays, fan_cones))

fan_rays = [1 0; 0 1; -1 0; -1 -1; 0 -1]
fan_cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,5],[5,1]])
dP2 = NormalToricVariety(PolyhedralFan(fan_rays, fan_cones))

fan_rays = [1 0; 1 1; 0 1; -1 0; -1 -1; 0 -1]
fan_cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]])
dP3 = NormalToricVariety(PolyhedralFan(fan_rays, fan_cones))

@testset "delPezzo surfaces" begin
    @test_throws ArgumentError del_pezzo(-1)
    @test length(torusinvariant_prime_divisors(dP0)) == 3
    @test length(torusinvariant_prime_divisors(dP1)) == 4
    @test length(torusinvariant_prime_divisors(dP2)) == 5
    @test rank(torusinvariant_cartier_divisor_group(dP3)) == 6
    @test length(torusinvariant_prime_divisors(dP3)) == 6
    @test_throws ArgumentError del_pezzo(4)
    @test rank(picard_group(dP3)) == 4
    @test picard_group(dP3) == codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(dP3))
end

blowup_variety = blowup_on_ith_minimal_torus_orbit(P2, 1, "e")

@testset "Blowup of projective space" begin
    @test isnormal(blowup_variety) == true
    @test isaffine(blowup_variety) == false
    @test isprojective(blowup_variety) == true
    @test issmooth(blowup_variety) == true
    @test iscomplete(blowup_variety) == true
    @test hastorusfactor(blowup_variety) == false
    @test isorbifold(blowup_variety) == true
    @test issimplicial(blowup_variety) == true
    @test betti_number(blowup_variety, 0) == 1
    @test betti_number(blowup_variety, 1) == 0
    @test betti_number(blowup_variety, 2) == 2
    @test betti_number(blowup_variety, 3) == 0
    @test betti_number(blowup_variety, 4) == 1
    @test euler_characteristic(blowup_variety) == 4
    @test rank(picard_group(blowup_variety)) == 2
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
    @test betti_number(v, 0) == 1
    @test betti_number(v, 1) == 0
    @test betti_number(v, 2) == 3
    @test betti_number(v, 3) == 0
    @test betti_number(v, 4) == 4
    @test betti_number(v, 5) == 0
    @test betti_number(v, 6) == 3
    @test betti_number(v, 7) == 0
    @test betti_number(v, 8) == 1
end

@testset "ComparisonWithProjectiveSpace" begin
    @test is_projective_space(H5) == false
    @test is_projective_space(P2) == true
    @test is_projective_space(blowup_variety) == false
    @test is_projective_space(v) == false
    @test is_projective_space(ntv) == false
    @test is_projective_space(ntv2) == false
end

P = polarize(Polyhedron(Polymake.polytope.rand_sphere(5,60; seed=42)))
big_variety = NormalToricVariety(P)
set_coefficient_ring(big_variety, GF(13))

@testset "Additional test for Stanley-Reisner ideal" begin
    # this should run in a few seconds at most
    duration = @elapsed stanley_reisner_ideal(big_variety)
    @test duration < 10
    @test ngens(stanley_reisner_ideal(big_variety)) == 1648
end

D=ToricDivisor(H5, [0,0,0,0])
D2 = DivisorOfCharacter(H5, [1,2])

@testset "Toric divisors" begin
    @test dim(toric_variety(D)) == 2
    @test isprime(D) == false
    @test iscartier(D) == true
    @test isprincipal(D) == true
    @test is_basepoint_free(D) == true
    @test isample(D) == false
    @test is_very_ample(D) == false
    @test isnef(D) == true
    @test isintegral(D) == true
    @test is_q_cartier(D) == true
    @test coefficients(D) == [0,0,0,0]
    @test isprime(D2) == false
    @test iscartier(D2) == true
    @test isprincipal(D2) == true
    @test is_basepoint_free(D2) == true
    @test isample(D2) == false
    @test is_very_ample(D2) == false
    @test isnef(D2) == true
    @test isintegral(D2) == true
    @test is_q_cartier(D2) == true
    @test isprime(D2) == false
    @test coefficients(D2) == [1, 2, 9, -2]
    @test coefficients(D2+D2) == coefficients(2*D2)
    @test coefficients(D2-D2) == [0,0,0,0]
    @test (D == D2) == false
end

p = polyhedron(D)

@testset "Polytopes of toric divisors" begin
    @test dim(p) == 0
    @test ambient_dim(p) == 2
end

DC1 = ToricDivisorClass(H5, [0,0])
DC2 = ToricDivisorClass(H5, [1,2])

@testset "Toric divisor classes" begin
    @test istrivial(DC1) == true
    @test istrivial(2 * DC1 + DC2) == false
    @test toric_variety(DC1) === toric_variety(DC2)
    @test (divisor_class(DC1) == divisor_class(DC2)) == false
end

line_bundle = ToricLineBundle(dP3, [1,2,3,4])
line_bundle2 = ToricLineBundle(D2)

@testset "Toric line bundles" begin
    @test degree(line_bundle) == 10
    @test degree(line_bundle * line_bundle) == 20
    @test degree(line_bundle^(-1)) == -10
    @test divisor_class(line_bundle).coeff == AbstractAlgebra.matrix(ZZ, [1 2 3 4])
    @test dim(toric_variety(line_bundle)) == 2
    @test istrivial(line_bundle) == false
    @test is_basepoint_free(line_bundle) == false
    @test isample(line_bundle) == false
    @test is_very_ample(line_bundle) == false
    @test istrivial(structure_sheaf(dP3)) == true
    @test all_cohomologies(line_bundle) == [11,0,0]
    @test cohomology(line_bundle,0) == 11
end

torus_coordinate_ring = coordinate_ring_of_torus(dP3)

@testset "Characters, rational functions and global sections" begin
    @test ngens(torus_coordinate_ring.I) == 2
    @test length(basis_of_global_sections_via_rational_functions(line_bundle)) == 11
    @test length(basis_of_global_sections_via_rational_functions(line_bundle2)) == 1
    @test length(basis_of_global_sections_via_rational_functions(line_bundle2^2)) == length(basis_of_global_sections_via_homogeneous_component(line_bundle2^2))
end

P = convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1])
charges = zero_matrix(ZZ, 1, 3);
charges[1,1] = 1;
charges[1,2] = 1;
charges[1,3] = 1;

@testset "ToricVarieties from triangulations and GLSMs" begin
    @test length(NormalToricVarietyFromTriangulation(P)) == 2
    @test length(NormalToricVarietyFromGLSM(charges)) == 1
end
