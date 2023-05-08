using Oscar
using Test

@testset "Line bundle cohomologies and vanishing sets (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    P2 = projective_space(NormalToricVariety, 2; set_attributes)
    dP3 = del_pezzo_surface(NormalToricVariety, 3; set_attributes)
    F5 = hirzebruch_surface(NormalToricVariety, 5; set_attributes)
    
    l = toric_line_bundle(dP3, [1, 2, 3, 4])
    l2 = toric_line_bundle(divisor_of_character(F5, [1, 2]))
    l3 = toric_line_bundle(P2, [1])
    
    vs = vanishing_sets(dP3)
    vs2 = vanishing_sets(P2)
    R,_ = polynomial_ring(QQ, 3)
    
    @testset "Cohomology with cohomCalg on dP3" begin
        @test cohomology(l, 0) == 0
        @test all_cohomologies(l) == [0, 16, 0]
    end
    
    @testset "Cohomology with cohomCalg on F5" begin
        @test cohomology(l2, 0) == 1
        @test all_cohomologies(l2) == [1, 0, 0]
    end
    
    @testset "Toric vanishing sets of dP3" begin
        @test is_projective_space(toric_variety(vs[1])) == false
        @test length(polyhedra(vs[1])) == 1
        @test cohomology_index(vs[1]) == 0
        @test contains(vs[1], structure_sheaf(dP3)) == false
        @test contains(vs[2], structure_sheaf(dP3)) == true
        @test contains(vs[3], structure_sheaf(dP3)) == true
        @test contains(vs[1], l) == true
        @test contains(vs[2], l) == false
        @test contains(vs[3], l) == true
    end
    
    @testset "Toric vanishing sets of P2" begin
        @test contains(vs2[1], l3) == false
        @test contains(vs2[2], l3) == true
        @test contains(vs2[3], l3) == true
    end
    
    @testset "Membership of line bundles in vanishing sets of different varieties" begin
        @test contains(vs[1], l2) == false
        @test contains(vs[1], l3) == false
        @test contains(vs2[2], l) == false
        @test contains(vs2[2], l2) == false
    end
    
    @testset "Should fail" begin
        @test_throws ArgumentError coordinate_ring_of_torus(R, dP3)
        @test_throws ArgumentError character_to_rational_function(R, dP3, [1, 2])
    end
    
    @testset "Global sections from rational functions and homogeneous components" begin
        @test ngens(coordinate_ring_of_torus(dP3).I) == 2
        @test length(basis_of_global_sections(anticanonical_bundle(dP3)^2)) == 19
        @test length(basis_of_global_sections_via_rational_functions(canonical_bundle(dP3))) == 0
        @test length(basis_of_global_sections(canonical_bundle(dP3))) == 0
        @test length(basis_of_global_sections_via_rational_functions(l)) == 0
        @test length(basis_of_global_sections(l)) == 0
    end
end
