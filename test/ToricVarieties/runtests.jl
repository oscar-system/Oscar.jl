using Oscar
using Test

##################################################
# (1) Tests for affine toric varieties
##################################################

antv = AffineNormalToricVariety(Oscar.positive_hull([1 1; -1 1]))
antv2 = NormalToricVariety(Oscar.positive_hull([1 1; -1 1]))
antv3 = AffineNormalToricVariety(antv2)
antv4 = AffineNormalToricVariety(Oscar.positive_hull([1 0]))
antv5 = AffineNormalToricVariety(Oscar.positive_hull([1 0; 0 1]))
antv6 = NormalToricVariety([[1,0,0], [1,0,1], [1,1,1], [1,1,0]], [[1,2,3,4]])

set_coordinate_names(antv4, ["u"])
set_coordinate_names_of_torus(antv4, ["u1","u2"])

@testset "Affine toric varieties" begin
    @test is_smooth(antv) == false
    @test is_orbifold(antv) == true
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
    @test is_affine(antv2) == true
    @test ngens(toric_ideal(antv2)) == 1
    @test length(affine_open_covering(antv3)) == 1
    @test dim(cone(antv3)) == 2
end

@testset "Affine toric varieties with torusfactor" begin
    @test_throws ArgumentError set_coordinate_names(antv4, ["u1", "u2"])
    @test_throws ArgumentError set_coordinate_names_of_torus(antv4, ["u"])
    @test_throws ArgumentError ideal_of_linear_relations(antv4)
    @test_throws ArgumentError map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(antv4)
    @test_throws ArgumentError map_from_torusinvariant_cartier_divisor_group_to_picard_group(antv4)
    @test_throws ArgumentError betti_number(antv4, 0)
    @test dim_of_torusfactor(antv4) == 1
end

@testset "Errors from changes to coordinates and coefficients once variety is finalized" begin
    @test ngens(cox_ring(antv4)) == 1
    @test is_finalized(antv4) == true
    @test_throws ErrorException set_coordinate_names(antv4, ["u"])
    @test_throws ErrorException set_coordinate_names_of_torus(antv4, ["u1", "u2"])
end

@testset "Affine toric varieties with trivial toric ideal" begin
    @test iszero(toric_ideal(antv5)) == true
end

@testset "Affine toric varieties from fans" begin
    @test is_trivial(picard_group(antv6)) == true
end





##################################################
# (2) Cyclic quotient singularities
##################################################

cyc = CyclicQuotientSingularity(2,1)

@testset "Cyclic quotient singularities" begin
    @test is_affine(cyc) == true
    @test is_smooth( cyc ) == false
    @test is_simplicial( cyc ) == true
    @test continued_fraction_hirzebruch_jung( cyc )[1] == 2
    @test dual_continued_fraction_hirzebruch_jung(cyc)[1] == 2
end





##################################################
# (3) Normal toric varieties
##################################################

ntv = NormalToricVariety(Oscar.normal_fan(Oscar.cube(2)))
set_coordinate_names(ntv, ["x1","x2","y1","y2"])
ntv2 = NormalToricVariety(Oscar.cube(2))
ntv3 = NormalToricVarietyFromGLSM(matrix(ZZ, [[1,1,1]]))
ntv4 = NormalToricVarietiesFromStarTriangulations(convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1]))
ntv5 = NormalToricVariety(polarize(Polyhedron(Polymake.polytope.rand_sphere(5,60; seed=42))))

@testset "Normal toric varieties from fans, triangulations and GLSMs" begin
    @test is_complete(ntv) == true
    @test rank(torusinvariant_cartier_divisor_group(ntv)) == 4
    @test rank(domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(ntv))) == 4
    @test is_complete(ntv2) == true
    @test length(ntv3) == 1
    @test is_projective_space(ntv3[1]) == true
    @test length(ntv4) == 2
end

@testset "Speed test for Stanley-Reisner ideal (at most a few seconds)" begin
    duration = @elapsed stanley_reisner_ideal(ntv5)
    @test duration < 10
    @test ngens(stanley_reisner_ideal(ntv5)) == 1648
end





##################################################
# (4) Projective space
##################################################

P2 = NormalToricVariety(normal_fan(Oscar.simplex(2)))
P2v2 = projective_space(NormalToricVariety,2)

@testset "Projective space P2" begin
    @test is_normal(P2) == true
    @test is_affine(P2) == false
    @test is_projective(P2) == true
    @test is_smooth(P2) == true
    @test is_complete(P2) == true
    @test has_torusfactor(P2) == false
    @test is_orbifold(P2) == true
    @test is_simplicial(P2) == true
    @test betti_number(P2, 0) == 1
    @test betti_number(P2, 1) == 0
    @test betti_number(P2, 2) == 1
    @test betti_number(P2, 3) == 0
    @test betti_number(P2, 4) == 1
    @test ngens(cox_ring(P2)) == 3
    @test length(stanley_reisner_ideal(P2).gens) == 1
    @test length(irrelevant_ideal(P2).gens) == 3
end

@testset "Test constructor for Hirzebruch surfaces" begin
  @test vcat([map_from_torusinvariant_weil_divisor_group_to_class_group(P2v2)(x).coeff for x in gens(torusinvariant_weil_divisor_group(P2v2))]) == matrix(ZZ, [[1],[1],[1]])
  @test coordinate_names(P2v2) == ["x1","x2","x3"]
end




##################################################
# (5) Hirzebruch surfaces
##################################################

F0 = hirzebruch_surface(0)
F5 = NormalToricVariety([[1,0], [0,1], [-1,5], [0,-1]], [[1,2],[2,3],[3,4],[4,1]])
F5v2 = hirzebruch_surface(5)

@testset "Hirzebruch surface F0" begin
    @test is_fano(F0) == true
end

@testset "Hirzebruch surface F5" begin
    @test is_normal(F5) == true
    @test is_affine(F5) == false
    @test is_projective(F5) == true
    @test is_smooth(F5) == true
    @test is_complete(F5) == true
    @test has_torusfactor(F5) == false
    @test is_orbifold(F5) == true
    @test is_simplicial(F5) == true
    @test is_gorenstein(F5) == true
    @test is_q_gorenstein(F5) == true
    @test is_fano(F5) == false
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

@testset "Test constructor for Hirzebruch surfaces" begin
  @test vcat([map_from_torusinvariant_weil_divisor_group_to_class_group(F5v2)(x).coeff for x in gens(torusinvariant_weil_divisor_group(F5v2))]) == matrix(ZZ, [[0,1],[1,0],[0,1],[1,5]])
  @test coordinate_names(F5v2) == ["t1","x1","t2","x2"]
end



##################################################
# (6) del Pezzo surfaces
##################################################

dP0 = NormalToricVariety(normal_fan(Oscar.simplex(2)))
dP1 = NormalToricVariety([[1,0], [0,1], [-1,0], [-1,-1]], [[1,2],[2,3],[3,4],[4,1]])
dP2 = NormalToricVariety([[1,0], [0,1], [-1,0], [-1,-1], [0,-1]], [[1,2],[2,3],[3,4],[4,5],[5,1]])
dP3 = NormalToricVariety([[1,0], [1,1], [0,1], [-1,0], [-1,-1], [0,-1]], [[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]])
set_coordinate_names(dP3,["x1","e1","x2","e3","x3","e2"])
dP3v2 = del_pezzo(3)

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

@testset "Test constructor for delPezzo surfaces" begin
  @test vcat([map_from_torusinvariant_weil_divisor_group_to_class_group(dP3v2)(x).coeff for x in gens(torusinvariant_weil_divisor_group(dP3v2))]) == matrix(ZZ, [[1,1,1,0],[1,1,0,1],[1,0,1,1],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]])
  @test coordinate_names(dP3v2) == ["x1","x2","x3","e1","e2","e3"]
end




##################################################
# (7) Blowup from star subdivision
##################################################

blowup_variety = blowup_on_ith_minimal_torus_orbit(P2, 1, "e")

@testset "Blowup of projective space" begin
    @test is_normal(blowup_variety) == true
    @test is_affine(blowup_variety) == false
    @test is_projective(blowup_variety) == true
    @test is_smooth(blowup_variety) == true
    @test is_complete(blowup_variety) == true
    @test has_torusfactor(blowup_variety) == false
    @test is_orbifold(blowup_variety) == true
    @test is_simplicial(blowup_variety) == true
    @test betti_number(blowup_variety, 0) == 1
    @test betti_number(blowup_variety, 1) == 0
    @test betti_number(blowup_variety, 2) == 2
    @test betti_number(blowup_variety, 3) == 0
    @test betti_number(blowup_variety, 4) == 1
    @test euler_characteristic(blowup_variety) == 4
    @test rank(picard_group(blowup_variety)) == 2
end





##################################################
# (8) Direct product
##################################################

ntv6 = F5 * P2

@testset "Direct product of toric varieties" begin
    @test is_normal(ntv6) == true
    @test is_affine(ntv6) == false
    @test is_projective(ntv6) == true
    @test is_smooth(ntv6) == true
    @test is_complete(ntv6) == true
    @test has_torusfactor(ntv6) == false
    @test is_orbifold(ntv6) == true
    @test is_simplicial(ntv6) == true
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





##################################################
# (9) Comparison with projective space
##################################################

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





##################################################
# (10) Toric divisors
##################################################

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
    @test is_prime(D) == false
    @test is_cartier(D) == true
    @test is_principal(D) == true
    @test is_trivial(D) == true
    @test is_basepoint_free(D) == true
    @test is_ample(D) == false
    @test is_very_ample(D) == false
    @test is_nef(D) == true
    @test is_integral(D) == true
    @test is_q_cartier(D) == true
    @test is_prime(D) == false
    @test dim(toric_variety(D)) == 2
    @test coefficients(D) == [0,0,0,0]
    @test dim(polyhedron(D)) == 0
    @test ambient_dim(polyhedron(D)) == 2
end

@testset "Toric divisor D2" begin
    @test is_prime(D2) == false
    @test is_cartier(D2) == true
    @test is_principal(D2) == true
    @test is_basepoint_free(D2) == true
    @test is_ample(D2) == false
    @test is_very_ample(D2) == false
    @test is_nef(D2) == true
    @test is_integral(D2) == true
    @test is_q_cartier(D2) == true
    @test is_prime(D2) == false
    @test coefficients(D2) == [1, 2, 9, -2]
end

@testset "Toric divisor D3" begin
    @test is_prime(D3) == true
end

@testset "Arithmetics of toric divisors" begin
    @test (D == D2) == false
    @test (D4 + D5 == D6) == true
    @test is_principal(fmpz(2)*D+D2) == true
    @test is_principal(2*D-D2) == true
    @test coefficients(D2+D2) == coefficients(2*D2)
    @test coefficients(D2-D2) == [0,0,0,0]
end





##################################################
# (11) Toric divisor classes
##################################################

DC = ToricDivisorClass(F5, [fmpz(0),fmpz(0)])
DC2 = ToricDivisorClass(F5, [1,2])
DC3 = ToricDivisorClass(dP3, [4,3,2,1])
DC4 = canonical_divisor_class(dP3)
DC5 = anticanonical_divisor_class(dP3)
DC6 = trivial_divisor_class(dP3)

@testset "Attributes of toric divisor classes" begin
    is_trivial(toric_divisor(DC2)) == false
    rank(parent(divisor_class(DC2))) == 2
    dim(toric_variety(DC2)) == 2
end

@testset "Arithmetics of toric divisor classes" begin
    @test is_trivial(fmpz(2)*DC+DC2) == false
    @test is_trivial(2*DC-DC2) == false
    @test (DC == DC2) == false
    @test (DC4 - DC5 == DC6) == false
    @test (DC == DC3) == false
end





##################################################
# (12) Toric line bundles
##################################################

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
    @test is_trivial(l) == false
    @test is_basepoint_free(l) == false
    @test is_ample(l) == false
    @test is_very_ample(l) == false
    @test degree(l) == 10
    @test divisor_class(l).coeff == AbstractAlgebra.matrix(ZZ, [1 2 3 4])
    @test dim(toric_variety(l)) == 2
end

@testset "Arithmetics of toric line bundles" begin
    @test (l == l3) == false
    @test (l4 * l5 == l7) == true
    @test (l * l6 * inv(l) == l7) == true
    @test degree(l^(-1)) == -10
    @test degree(l*l) == 20
end





##################################################
# (13) Line bundle cohomologies and vanishing sets
##################################################

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





##################################################
# (14) Topological intersection numbers
##################################################

(u1,u2,u3,u4) = gens(cohomology_ring(dP1))
(x1,e1,x2,e3,x3,e2) = gens(cohomology_ring(dP3))
c = CohomologyClass(dP3, x1)
c2 = CohomologyClass(dP3, e1)
c3 = CohomologyClass(dP1, u1)

@testset "Argument errors for cohomology classes" begin
    @test_throws ArgumentError cohomology_ring(antv)
    @test_throws ArgumentError chow_ring(antv)
    @test_throws ArgumentError is_trivial(c - c3)
    @test_throws ArgumentError is_trivial(c + c3)
    @test_throws ArgumentError ideal_of_linear_relations(R, dP3)
end

@testset "Cohomology ring, Chow ring and volume form of direct product space" begin
    @test ngens(ideal_of_linear_relations(ntv6)) == 4
    @test ngens(chow_ring(ntv6).I) == 7
    @test is_trivial(volume_form(ntv6)) == false
end

@testset "Properties, attributes and arithmetics of cohomology classes" begin
    @test is_trivial(c) == false
    @test nrows(exponents(c)) == 3
    @test length(coefficients(c)) == 3
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





##################################################
# (15) Closd subvarieties
##################################################

(x1, x2, y1, y2) = gens(cox_ring(ntv));
sv1 = ClosedSubvarietyOfToricVariety(ntv, [x1])
sv2 = ClosedSubvarietyOfToricVariety(ntv, [x1^2+x1*x2+x2^2,y2])

@testset "Test error messages for closed subvarieties" begin
    @test_throws ArgumentError ClosedSubvarietyOfToricVariety(ntv, [x1 - y1])
    @test_throws ArgumentError ClosedSubvarietyOfToricVariety(antv6, [gens(cox_ring(antv6))[1]])
end

@testset "Properties and attributes of closd subvarieties" begin
    @test is_empty(sv1) == false
    @test radical(sv1) == defining_ideal(sv1)
    @test dim(toric_variety(sv1)) == 2
end




##################################################
# (16) Algebraic cycles
##################################################

ac1 = RationalEquivalenceClass(D2)
ac2 = RationalEquivalenceClass(DC3)
ac3 = RationalEquivalenceClass(l4)
ac4 = RationalEquivalenceClass(c3)
ac5 = RationalEquivalenceClass(CohomologyClass(dP3, x3))

@testset "Error of rational equivalence classes" begin
  @test_throws ArgumentError RationalEquivalenceClass(antv, [1,2,3])
  @test_throws ArgumentError RationalEquivalenceClass(toric_variety(ac1),[1,2,3])
  @test_throws ArgumentError ac1 + ac3
  @test_throws ArgumentError ac1 - ac3
  @test_throws ArgumentError ac1 * ac3
end

@testset "Arithmetics of rational equivalence cycles" begin
    @test (ac1 == ac2) == false
    @test is_trivial(ac2 - ac3) == false
    @test is_trivial(ac1) == true
end

@testset "Attributes of rational equivalence cycles" begin
    @test dim(toric_variety(ac1)) == 2
    @test polynomial(ac1) == 0
    @test parent(representant(ac1)) == cox_ring(toric_variety(ac1))
    @test is_trivial(cohomology_class(3*ac4)) == false
    @test is_trivial(RationalEquivalenceClass(sv2)) == false
    @test length(coefficients(ac1)) == 0
    @test length(components(ac2-ac3)) == 3
end

@testset "Intersection of rational equivalence cycles" begin
    @test is_trivial(ac2*ac2) == false
    @test is_trivial(ac2*ac2*ac2) == true
end
