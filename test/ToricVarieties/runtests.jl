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
    @test ngens(toric_ideal(ntv)) == 1
    @test dim(cone(associated_affine_variety)) == 2
end

C2 = Oscar.positive_hull([1 0])
antv2 = AffineNormalToricVariety(C2)
set_coefficient_ring(antv2, ZZ)

@testset "Affine toric varieties with torusfactor" begin
    @test_throws ArgumentError set_coordinate_names(antv2, ["u1", "u2"])
    @test_throws ArgumentError ideal_of_linear_relations(antv2)
    @test_throws ArgumentError set_coordinate_names_of_torus(antv2, ["u"])
    @test_throws ArgumentError map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(antv2)
    @test_throws ArgumentError map_from_torusinvariant_cartier_divisor_group_to_picard_group(antv2)
    @test_throws ArgumentError betti_number(antv2, 0)
    @test dim_of_torusfactor(antv2) == 1
end

C3 = Oscar.positive_hull([1 0; 0 1])
antv3 = AffineNormalToricVariety(C3)

@testset "Affine toric varieties with trivial toric ideal" begin
    @test iszero(toric_ideal(antv3)) == true
end

r = [1 0 0; 1 0 1; 1 1 1; 1 1 0]
c = IncidenceMatrix([[1,2,3,4]])
antv4 = NormalToricVariety(PolyhedralFan(r, c))

@testset "Affine toric varieties created from fan" begin
    @test istrivial(picard_group(antv4)) == true
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
    @test ngens(cox_ring(P2)) == 3
    @test length(stanley_reisner_ideal(P2).gens) == 1
    @test length(irrelevant_ideal(P2).gens) == 3
end

fan_rays = [1 0; 0 1; -1 5; 0 -1]
fan_cones = IncidenceMatrix([[1,2],[2,3],[3,4],[4,1]])
H5 = NormalToricVariety(PolyhedralFan(fan_rays, fan_cones))
nc = nef_cone(H5)
mc = mori_cone(H5)
mapping = map_from_character_lattice_to_torusinvariant_weil_divisor_group(H5)

@testset "Hirzebruch surfaces" begin
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
    @test dim(nc) == 2
    @test dim(mc) == 2
    @test rank(domain(mapping)) == 2
    @test rank(codomain(mapping)) == 4
    @test rank(class_group(H5)) == 2
    @test rank(codomain(map_from_torusinvariant_weil_divisor_group_to_class_group(H5))) == 2
    @test ngens(cox_ring(H5)) == 4
    @test length(stanley_reisner_ideal(H5).gens) == 2
    @test length(irrelevant_ideal(H5).gens) == 4
    @test isfano(hirzebruch_surface(0)) == true
    @test isfano(hirzebruch_surface(7)) == false
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

R,_ = PolynomialRing(QQ, 3)

@testset "delPezzo surfaces" begin
    @test_throws ArgumentError del_pezzo(-1)
    @test_throws ArgumentError del_pezzo(4)
    @test length(torusinvariant_prime_divisors(dP0)) == 3
    @test length(torusinvariant_prime_divisors(dP1)) == 4
    @test length(torusinvariant_prime_divisors(dP2)) == 5
    @test rank(torusinvariant_cartier_divisor_group(dP3)) == 6
    @test length(torusinvariant_prime_divisors(dP3)) == 6
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
    @test is_projective_space(antv) == false
    @test is_projective_space(H5) == false
    @test is_projective_space(P2) == true
    @test is_projective_space(del_pezzo(0)) == true
    @test is_projective_space(del_pezzo(1)) == false
    @test is_projective_space(del_pezzo(2)) == false
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
D3 = ToricDivisor(dP3, [1,0,0,0,0,0])

@testset "Toric divisors" begin
    @test_throws ArgumentError ToricDivisor(H5, [0,0,0])
    @test_throws ArgumentError D+D3
    @test_throws ArgumentError D-D3
    @test isprincipal(fmpz(2)*D+D2) == true
    @test isprincipal(2*D-D2) == true
    @test (D == D2) == false
    @test dim(toric_variety(D)) == 2
    @test isprime(D) == false
    @test iscartier(D) == true
    @test isprincipal(D) == true
    @test istrivial(D) == true
    @test is_basepoint_free(D) == true
    @test isample(D) == false
    @test is_very_ample(D) == false
    @test isnef(D) == true
    @test isintegral(D) == true
    @test is_q_cartier(D) == true
    @test isprime(D) == false
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
    @test D != D2
    @test canonical_divisor(dP3) - canonical_divisor(dP3) == trivial_divisor(dP3)
    @test anticanonical_divisor(dP3) + canonical_divisor(dP3) == trivial_divisor(dP3)
    @test isprime(D3) == true
end

p = polyhedron(D)

@testset "Polytopes of toric divisors" begin
    @test dim(p) == 0
    @test ambient_dim(p) == 2
end

DC = ToricDivisorClass(H5, [fmpz(0),fmpz(0)])
DC2 = ToricDivisorClass(H5, [1,2])
DC3 = ToricDivisorClass(dP3, [4,3,2,1])

@testset "Toric divisor classes" begin
    @test_throws ArgumentError (DC == DC3)
    @test istrivial(fmpz(2)*DC+DC2) == false
    @test istrivial(2*DC-DC2) == false
    @test (DC == DC2) == false
    @test canonical_divisor_class(dP3) - canonical_divisor_class(dP3) == trivial_divisor_class(dP3)
    @test anticanonical_divisor_class(dP3) + canonical_divisor_class(dP3) == trivial_divisor_class(dP3)
end

line_bundle = ToricLineBundle(dP3, [1,2,3,4])
line_bundle2 = ToricLineBundle(D2)
line_bundle3 = ToricLineBundle(dP1, [fmpz(1),fmpz(2)])

@testset "Toric line bundles" begin
    @test_throws ArgumentError line_bundle == line_bundle3
    @test_throws ArgumentError line_bundle * line_bundle3
    @test istrivial(line_bundle) == false
    @test istrivial(line_bundle^2) == false
    @test istrivial(inv(line_bundle)*line_bundle) == true
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
    @test contains(vanishing_sets(dP3)[1], line_bundle) == false
    @test contains(vanishing_sets(dP3)[2], line_bundle) == true
    @test contains(vanishing_sets(dP3)[3], line_bundle) == true
    @test istrivial(canonical_bundle(dP3)) == false
    @test istrivial(structure_sheaf(dP3)) == true
    @test inv(anticanonical_bundle(dP3)) == canonical_bundle(dP3)
end

torus_coordinate_ring = coordinate_ring_of_torus(dP3)

@testset "Characters, rational functions and global sections" begin
    @test_throws ArgumentError coordinate_ring_of_torus(R, dP3)
    @test_throws ArgumentError character_to_rational_function(R, dP3, [1, 2])
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
    @test length(NormalToricVarietiesFromStarTriangulations(P)) == 2
    @test length(NormalToricVarietyFromGLSM(charges)) == 1
    @test is_projective_space(NormalToricVarietyFromGLSM([[fmpz(1), fmpz(1), fmpz(1)]])[1]) == true
end

(x1,e1,x2,e3,x3,e2) = gens(cohomology_ring(dP3))
c = CohomologyClass(dP3, x1)
c2 = CohomologyClass(dP3, e1)
(u1,u2,u3,u4) = gens(cohomology_ring(dP1))
c3 = CohomologyClass(dP1, u1)

@testset "Topological intersection numbers" begin
    @test_throws ArgumentError cohomology_ring(antv)
    @test_throws ArgumentError chow_ring(antv)
    @test_throws ArgumentError c == c3
    @test_throws ArgumentError istrivial(c - c3)
    @test_throws ArgumentError istrivial(c + c3)
    @test_throws ArgumentError ideal_of_linear_relations(R, dP3)
    @test 2 * c != fmpz(3) * c2
    @test fmpq(3) * c == fmpz(3) * c
    @test ngens(ideal_of_linear_relations(v)) == 4
    @test ngens(chow_ring(v).I) == 7
    @test integrate(volume_form(v)) == 1
    @test nrows(exponents(c)) == 1
    @test length(coefficients(c)) == 1
    @test istrivial(c) == false
    @test integrate(c^2+c-3//4*c*c) == -1//4
    @test integrate(CohomologyClass(dP3,e1*e1)) == -1
    @test integrate(CohomologyClass(dP3,e2*e2)) == -1
    @test integrate(CohomologyClass(dP3,e3*e3)) == -1
    @test integrate(CohomologyClass(dP3,x1*x1)) == -1
    @test integrate(CohomologyClass(dP3,x2*x2)) == -1
    @test integrate(CohomologyClass(dP3,x3*x3)) == -1
    @test length(intersection_form(dP3)) == 21
    @test integrate(c^2+c-3//4*c*c) == -1//4
    @test integrate(c) == 0
end
