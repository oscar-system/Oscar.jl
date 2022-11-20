@testset "Line bundle cohomologies and vanishing sets" begin
    dP3 = NormalToricVariety([[1, 0], [1, 1], [0, 1], [-1, 0], [-1, -1], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]])
    F5 = NormalToricVariety([[1, 0], [0, 1], [-1, 5], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    D2 = DivisorOfCharacter(F5, [1, 2])

    l = ToricLineBundle(dP3, [1, 2, 3, 4])
    l2 = ToricLineBundle(D2)
    l4 = canonical_bundle(dP3)
    l5 = anticanonical_bundle(dP3)
    l7 = structure_sheaf(dP3)

    vs = vanishing_sets(dP3)
    R,_ = PolynomialRing(QQ, 3)

    @testset "Cohomology with cohomCalg" begin
        @test cohomology(l, 0) == 11
        @test all_cohomologies(l) == [11, 0, 0]
        @test cohomology(l2, 0) == 1
        @test all_cohomologies(l2) == [1, 0, 0]
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

    @testset "Should fail" begin
        @test_throws ArgumentError coordinate_ring_of_torus(R, dP3)
        @test_throws ArgumentError character_to_rational_function(R, dP3, [1, 2])
    end

    @testset "Global sections from rational functions and homogeneous components" begin
        @test ngens(coordinate_ring_of_torus(dP3).I) == 2
        @test length(basis_of_global_sections_via_rational_functions(l)) == 11
        @test length(basis_of_global_sections(l)) == 11
        @test length(basis_of_global_sections(l5^2)) == 19
        @test length(basis_of_global_sections_via_rational_functions(l4)) == 0
        @test length(basis_of_global_sections(l4)) == 0
    end

end
