using Oscar
using Test

#####################
# (1) Tests for affine toric varieties
#####################

antv = AffineNormalToricVariety(Oscar.positive_hull([1 1; -1 1]))
antv2 = NormalToricVariety(Oscar.positive_hull([1 1; -1 1]))
antv3 = AffineNormalToricVariety(antv2)
antv4 = AffineNormalToricVariety(Oscar.positive_hull([1 0]))
antv5 = AffineNormalToricVariety(Oscar.positive_hull([1 0; 0 1]))
antv6 = NormalToricVariety([[1,0,0], [1,0,1], [1,1,1], [1,1,0]], [[1,2,3,4]])

set_coefficient_ring(antv4, ZZ)

@testset "Affine toric varieties" begin
    @test issmooth(antv) == false
    @test isorbifold(antv) == true
    @test dim(fan( antv )) == 2
    @test dim(cone(antv)) == 2
    @test length(affine_open_covering(antv)) == 1
    @test length(gens(toric_ideal(antv))) == 1
    @test rank(torusinvariant_weil_divisor_group(antv)) == 2
    @test rank(character_lattice(antv)) == 2
    @test rank(domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(antv))) == 2
    @test rank(codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(antv))) == 2
    @test elementary_divisors(codomain(map_from_torusinvariant_weil_divisor_group_to_class_group(antv))) == [ 2 ]
    @test elementary_divisors(class_group(antv)) == [ 2 ]
    @test ngens(cox_ring(antv)) == 2
    @test length(torusinvariant_prime_divisors(antv)) == 2
end

@testset "Affine toric varieties as normal toric varieties" begin
    @test isaffine(antv2) == true
    @test ngens(toric_ideal(antv2)) == 1
    @test length(affine_open_covering(antv3)) == 1
    @test dim(cone(antv3)) == 2
end

@testset "Affine toric varieties with torusfactor" begin
    @test_throws ArgumentError set_coordinate_names(antv4, ["u1", "u2"])
    @test_throws ArgumentError ideal_of_linear_relations(antv4)
    @test_throws ArgumentError set_coordinate_names_of_torus(antv4, ["u"])
    @test_throws ArgumentError map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(antv4)
    @test_throws ArgumentError map_from_torusinvariant_cartier_divisor_group_to_picard_group(antv4)
    @test_throws ArgumentError betti_number(antv4, 0)
    @test dim_of_torusfactor(antv4) == 1
end

@testset "Affine toric varieties with trivial toric ideal" begin
    @test iszero(toric_ideal(antv5)) == true
end

@testset "Affine toric varieties from fans" begin
    @test istrivial(picard_group(antv6)) == true
end





#####################
# (2) Cyclic quotient singularities
#####################

cyc = CyclicQuotientSingularity(2,1)

@testset "Cyclic quotient singularities" begin
    @test isaffine(cyc) == true
    @test issmooth( cyc ) == false
    @test issimplicial( cyc ) == true
    @test continued_fraction_hirzebruch_jung( cyc )[1] == 2
    @test dual_continued_fraction_hirzebruch_jung(cyc)[1] == 2
end





#####################
# (3) Normal toric varieties
#####################

ntv = NormalToricVariety(Oscar.normal_fan(Oscar.cube(2)))
ntv2 = NormalToricVariety(Oscar.cube(2))
ntv3 = NormalToricVarietyFromGLSM(matrix(ZZ, [[1,1,1]]))
ntv4 = NormalToricVarietiesFromStarTriangulations(convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1]))
ntv5 = NormalToricVariety(polarize(Polyhedron(Polymake.polytope.rand_sphere(5,60; seed=42))))

set_coefficient_ring(ntv5, GF(13))

@testset "Normal toric varieties from fans, triangulations and GLSMs" begin
    @test iscomplete(ntv) == true
    @test rank(torusinvariant_cartier_divisor_group(ntv)) == 4
    @test rank(domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(ntv))) == 4
    @test iscomplete(ntv2) == true
    @test length(ntv3) == 1
    @test is_projective_space(ntv3[1]) == true
    @test length(ntv4) == 2
end

@testset "Speed test for Stanley-Reisner ideal (at most a few seconds)" begin
    duration = @elapsed stanley_reisner_ideal(ntv5)
    @test duration < 10
    @test ngens(stanley_reisner_ideal(ntv5)) == 1648
end





#####################
# (4) Projective space
#####################

P2 = NormalToricVariety(normal_fan(Oscar.simplex(2)))

@testset "Projective space P2" begin
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





#####################
# (5) Hirzebruch surfaces
#####################

F0 = hirzebruch_surface(0)
F5 = NormalToricVariety([[1,0], [0,1], [-1,5], [0,-1]], [[1,2],[2,3],[3,4],[4,1]])

@testset "Hirzebruch surface F0" begin
    @test isfano(F0) == true
end

@testset "Hirzebruch surface F5" begin
    @test isnormal(F5) == true
    @test isaffine(F5) == false
    @test isprojective(F5) == true
    @test issmooth(F5) == true
    @test iscomplete(F5) == true
    @test hastorusfactor(F5) == false
    @test isorbifold(F5) == true
    @test issimplicial(F5) == true
    @test isgorenstein(F5) == true
    @test is_q_gorenstein(F5) == true
    @test isfano(F5) == false
    @test dim(F5) == 2
    @test dim_of_torusfactor(F5) == 0
    @test euler_characteristic(F5) == 4
    @test betti_number(F5, 0) == 1
    @test betti_number(F5, 1) == 0
    @test betti_number(F5, 2) == 2
    @test betti_number(F5, 3) == 0
    @test betti_number(F5, 4) == 1
    @test length(affine_open_covering(F5)) == 4
    @test fan(F5).pm_fan.FAN_DIM == 2
    @test rank(torusinvariant_weil_divisor_group(F5)) == 4
    @test rank(character_lattice(F5)) == 2
    @test dim(nef_cone(F5)) == 2
    @test dim(mori_cone(F5)) == 2
    @test rank(domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5))) == 2
    @test rank(codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5))) == 4
    @test rank(class_group(F5)) == 2
    @test rank(codomain(map_from_torusinvariant_weil_divisor_group_to_class_group(F5))) == 2
    @test ngens(cox_ring(F5)) == 4
    @test length(stanley_reisner_ideal(F5).gens) == 2
    @test length(irrelevant_ideal(F5).gens) == 4
end





#####################
# (6) del Pezzo surfaces
#####################

dP0 = NormalToricVariety(normal_fan(Oscar.simplex(2)))
dP1 = NormalToricVariety([[1,0], [0,1], [-1,0], [-1,-1]], [[1,2],[2,3],[3,4],[4,1]])
dP2 = NormalToricVariety([[1,0], [0,1], [-1,0], [-1,-1], [0,-1]], [[1,2],[2,3],[3,4],[4,5],[5,1]])
dP3 = NormalToricVariety([[1,0], [1,1], [0,1], [-1,0], [-1,-1], [0,-1]], [[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]])

@testset "Argument errors for del Pezzo surfaces" begin
    @test_throws ArgumentError del_pezzo(-1)
    @test_throws ArgumentError del_pezzo(4)
end

@testset "Attributes of del Pezzo surfaces" begin
    @test length(torusinvariant_prime_divisors(dP0)) == 3
    @test length(torusinvariant_prime_divisors(dP1)) == 4
    @test length(torusinvariant_prime_divisors(dP2)) == 5
    @test rank(torusinvariant_cartier_divisor_group(dP3)) == 6
    @test length(torusinvariant_prime_divisors(dP3)) == 6
    @test rank(picard_group(dP3)) == 4
    @test picard_group(dP3) == codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(dP3))
end





#####################
# (7) Blowup from star subdivision
#####################

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





#####################
# (8) Direct product
#####################

ntv6 = F5 * P2

@testset "Direct product of toric varieties" begin
    @test isnormal(ntv6) == true
    @test isaffine(ntv6) == false
    @test isprojective(ntv6) == true
    @test issmooth(ntv6) == true
    @test iscomplete(ntv6) == true
    @test hastorusfactor(ntv6) == false
    @test isorbifold(ntv6) == true
    @test issimplicial(ntv6) == true
    @test betti_number(ntv6, 0) == 1
    @test betti_number(ntv6, 1) == 0
    @test betti_number(ntv6, 2) == 3
    @test betti_number(ntv6, 3) == 0
    @test betti_number(ntv6, 4) == 4
    @test betti_number(ntv6, 5) == 0
    @test betti_number(ntv6, 6) == 3
    @test betti_number(ntv6, 7) == 0
    @test betti_number(ntv6, 8) == 1
end





########################
# (9) Comparison with projective space
########################

@testset "ComparisonWithProjectiveSpace" begin
    @test is_projective_space(antv) == false
    @test is_projective_space(antv2) == false
    @test is_projective_space(antv3) == false
    @test is_projective_space(antv4) == false
    @test is_projective_space(antv5) == false
    @test is_projective_space(antv6) == false
    @test is_projective_space(cyc) == false
    @test is_projective_space(ntv) == false
    @test is_projective_space(ntv2) == false
    @test is_projective_space(ntv3[1]) == true
    @test is_projective_space(ntv4[1]) == false
    @test is_projective_space(ntv4[2]) == false
    @test is_projective_space(ntv5) == false
    @test is_projective_space(P2) == true
    @test is_projective_space(F0) == false
    @test is_projective_space(F5) == false
    @test is_projective_space(dP0) == true
    @test is_projective_space(dP1) == false
    @test is_projective_space(dP2) == false
    @test is_projective_space(dP3) == false
    @test is_projective_space(blowup_variety) == false
    @test is_projective_space(ntv6) == false
end





########################
# (10) Toric divisors
########################

D=ToricDivisor(F5, [0,0,0,0])
D2 = DivisorOfCharacter(F5, [1,2])
D3 = ToricDivisor(dP3, [1,0,0,0,0,0])
D4 = canonical_divisor(dP3)
D5 = anticanonical_divisor(dP3)
D6 = trivial_divisor(dP3)

@testset "Argument errors for toric divisors" begin
    @test_throws ArgumentError ToricDivisor(F5, [0,0,0])
    @test_throws ArgumentError D+D3
    @test_throws ArgumentError D-D3
end

@testset "Toric divisor D" begin
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
    @test dim(toric_variety(D)) == 2
    @test coefficients(D) == [0,0,0,0]
    @test dim(polyhedron(D)) == 0
    @test ambient_dim(polyhedron(D)) == 2
end

@testset "Toric divisor D2" begin
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
end

@testset "Toric divisor D3" begin
    @test isprime(D3) == true
end

@testset "Arithmetics of toric divisors" begin
    @test (D == D2) == false
    @test (D4 + D5 == D6) == true
    @test isprincipal(fmpz(2)*D+D2) == true
    @test isprincipal(2*D-D2) == true
    @test coefficients(D2+D2) == coefficients(2*D2)
    @test coefficients(D2-D2) == [0,0,0,0]
end





########################
# (11) Toric divisor classes
########################

DC = ToricDivisorClass(F5, [fmpz(0),fmpz(0)])
DC2 = ToricDivisorClass(F5, [1,2])
DC3 = ToricDivisorClass(dP3, [4,3,2,1])
DC4 = canonical_divisor_class(dP3)
DC5 = anticanonical_divisor_class(dP3)
DC6 = trivial_divisor_class(dP3)

@testset "Attributes of toric divisor classes" begin
    istrivial(toric_divisor(DC2)) == false
    rank(parent(divisor_class(DC2))) == 2
    dim(toric_variety(DC2)) == 2
end

@testset "Arithmetics of toric divisor classes" begin
    @test istrivial(fmpz(2)*DC+DC2) == false
    @test istrivial(2*DC-DC2) == false
    @test (DC == DC2) == false
    @test (DC4 - DC5 == DC6) == false
    @test (DC == DC3) == false
end





########################
# (12) Toric line bundles
########################

l = ToricLineBundle(dP3, [1,2,3,4])
l2 = ToricLineBundle(D2)
l3 = ToricLineBundle(dP1, [fmpz(1),fmpz(2)])
l4 = canonical_bundle(dP3)
l5 = anticanonical_bundle(dP3)
l6 = ToricLineBundle(dP3, trivial_divisor(dP3))
l7 = structure_sheaf(dP3)

@testset "Argument errors for toric line bundles" begin
    @test_throws ArgumentError l * l3
end

@testset "Properties and attributes of toric line bundles" begin
    @test istrivial(l) == false
    @test is_basepoint_free(l) == false
    @test isample(l) == false
    @test is_very_ample(l) == false
    @test degree(l) == 10
    @test divisor_class(l).coeff == AbstractAlgebra.matrix(ZZ, [1 2 3 4])
    @test dim(toric_variety(l)) == 2
end

@testset "Arithmetics of torlc line bundles" begin
    @test (l == l3) == false
    @test (l4 * l5 == l7) == true
    @test (l * l6 * inv(l) == l7) == true
    @test degree(l^(-1)) == -10
    @test degree(l*l) == 20
end





################################
# (13) Line bundle cohomologies and vanishing sets
################################

vs = vanishing_sets(dP3)
R,_ = PolynomialRing(QQ, 3)

@testset "Line bundle cohomologies with cohomCalg" begin
    @test cohomology(l,0) == 11
    @test all_cohomologies(l) == [11,0,0]
    @test cohomology(l2,0) == 1
    @test all_cohomologies(l2) == [1,0,0]
end

@testset "Toric vanishingSets of dP3" begin
    @test is_projective_space(toric_variety(vs[1])) == false
    @test contains(vs[1], l) == false
    @test contains(vs[1], l2) == false
    @test contains(vs[2], l) == true
    @test contains(vs[3], l) == true
    @test contains(vs[1], l7) == false
    @test contains(vs[2], l7) == true
    @test contains(vs[3], l7) == true
    @test length(polyhedra(vs[1])) == 1
    @test cohomology_index(vs[1]) == 0
end

@testset "Argument errors for characters and rational functions" begin
    @test_throws ArgumentError coordinate_ring_of_torus(R, dP3)
    @test_throws ArgumentError character_to_rational_function(R, dP3, [1, 2])
end

@testset "Basis of global sections from rational functions and homogeneous components" begin
    @test ngens(coordinate_ring_of_torus(dP3).I) == 2
    @test length(basis_of_global_sections_via_rational_functions(l)) == 11
    @test length(basis_of_global_sections(l)) == 11
    @test length(basis_of_global_sections(l5^2)) == 19
    @test length(basis_of_global_sections_via_rational_functions(l4)) == 0
    @test length(basis_of_global_sections(l4)) == 0
end





#########################
# (14) Topological intersection numbers
#########################

(u1,u2,u3,u4) = gens(cohomology_ring(dP1))
(x1,e1,x2,e3,x3,e2) = gens(cohomology_ring(dP3))
c = CohomologyClass(dP3, x1)
c2 = CohomologyClass(dP3, e1)
c3 = CohomologyClass(dP1, u1)

@testset "Argument errors for cohomology classes" begin
    @test_throws ArgumentError cohomology_ring(antv)
    @test_throws ArgumentError chow_ring(antv)
    @test_throws ArgumentError istrivial(c - c3)
    @test_throws ArgumentError istrivial(c + c3)
    @test_throws ArgumentError ideal_of_linear_relations(R, dP3)
end

@testset "Cohomology ring, Chow ring and volume form of direct product space" begin
    @test ngens(ideal_of_linear_relations(ntv6)) == 4
    @test ngens(chow_ring(ntv6).I) == 7
    @test istrivial(volume_form(ntv6)) == false
end

@testset "Properties, attributes and arithmetics of cohomology classes" begin
    @test istrivial(c) == false
    @test nrows(exponents(c)) == 1
    @test length(coefficients(c)) == 1
    @test fmpq(3) * c == fmpz(3) * c
    @test 2 * c != fmpz(3) * c2
    @test (c == c3) == false
end

@testset "Computing topological intersection numbers" begin
    @test integrate(CohomologyClass(dP3,e1*e1)) == -1
    @test integrate(CohomologyClass(dP3,e2*e2)) == -1
    @test integrate(CohomologyClass(dP3,e3*e3)) == -1
    @test integrate(CohomologyClass(dP3,x1*x1)) == -1
    @test integrate(CohomologyClass(dP3,x2*x2)) == -1
    @test integrate(CohomologyClass(dP3,x3*x3)) == -1
    @test integrate(c) == 0
    @test integrate(c^2+c-3//4*c*c) == -1//4
    @test length(intersection_form(dP3)) == 21
    @test integrate(c^2+c-3//4*c*c) == -1//4
end
