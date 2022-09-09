@testset "describe() for permutation groups" begin
   @test describe(symmetric_group(1)) == "1"
   @test describe(symmetric_group(2)) == "C2"
   @test describe(symmetric_group(3)) == "S3"
   @test describe(symmetric_group(4)) == "S4"
   @test describe(symmetric_group(5)) == "S5"

   @test describe(alternating_group(1)) == "1"
   @test describe(alternating_group(2)) == "1"
   @test describe(alternating_group(3)) == "C3"
   @test describe(alternating_group(4)) == "A4"
   @test describe(alternating_group(5)) == "A5"
end

@testset "describe() for pc groups" begin
   @test describe(dihedral_group(2)) == "C2"
   @test describe(dihedral_group(4)) == "C2 x C2"
   @test describe(dihedral_group(6)) == "S3"
   @test describe(dihedral_group(8)) == "D8"
   @test describe(dihedral_group(10)) == "D10"
end

@testset "describe() for matrix groups" begin
   @test describe(general_linear_group(2,2)) == "S3"        # FIXME: correct but not ideal
   @test describe(general_linear_group(2,3)) == "GL(2,3)"
   @test describe(general_linear_group(2,4)) == "GL(2,4)"
   @test describe(general_linear_group(3,2)) == "PSL(3,2)"  # FIXME: correct but not ideal

   @test describe(special_linear_group(2,2)) == "S3"        # FIXME: correct but not ideal
   @test describe(special_linear_group(2,3)) == "SL(2,3)"
   @test describe(special_linear_group(2,4)) == "A5"        # FIXME: correct but not ideal
   @test describe(special_linear_group(3,2)) == "PSL(3,2)"  # FIXME: correct but not ideal
   @test describe(special_linear_group(3,3)) == "PSL(3,3)"  # FIXME: correct but not ideal
   @test describe(special_linear_group(4,2)) == "A8"        # FIXME: correct but not ideal
end

@testset "describe() for projective groups" begin
   # TODO: PGL / projective_general_linear_group constructor is missing
   # TODO: PSL / projective_general_linear_group constructor is missing
end

@testset "describe() for free groups" begin
   @test describe(free_group(0)) == "1"
   @test describe(free_group(1)) == "Z"
   @test describe(free_group(2)) == "a free group of rank 2"
   @test describe(free_group(3)) == "a free group of rank 3"

   F = free_group(4)
   subs = [sub(F, gens(F)[1:n])[1] for n in 0:4]
   @test describe(subs[1]) == "1"
   @test describe(subs[2]) == "Z"
   @test describe(subs[3]) == "a free group of rank 2"
   @test describe(subs[4]) == "a free group of rank 3"
   @test describe(subs[5]) == "a free group of rank 4"
end

@testset "describe() for finitely presented groups" begin

   F = free_group(1)
   @test describe(direct_product(F, F)) == "Z x Z"

   F = free_group(2)
   G = quo(F, [gen(F,2)])[1]
   @test describe(G) == "Z"
   G = quo(F, [gen(F,2)/gen(F,1)])[1]
   @test describe(G) == "Z"

   F = free_group(3)
   G, _ = quo(F, [gen(F,1)^2,gen(F,2)^2])
   @test describe(G) == "a finitely presented infinite group"
   is_abelian(G)
   @test describe(G) == "a finitely presented non-abelian infinite group"

   @test describe(derived_subgroup(free_group(2))[1]) == "a non-finitely generated group"

end


@testset "describe() for some real world examples" begin

   P = SimplicialComplex([[1,2,4,9],[1,2,4,15],[1,2,6,14],[1,2,6,15],[1,2,9,14],[1,3,4,12],[1,3,4,15],[1,3,7,10],[1,3,7,12],[1,3,10,15],[1,4,9,12],[1,5,6,13],[1,5,6,14],[1,5,8,11],[1,5,8,13],[1,5,11,14],[1,6,13,15],[1,7,8,10],[1,7,8,11],[1,7,11,12],[1,8,10,13],[1,9,11,12],[1,9,11,14],[1,10,13,15],[2,3,5,10],[2,3,5,11],[2,3,7,10],[2,3,7,13],[2,3,11,13],[2,4,9,13],[2,4,11,13],[2,4,11,15],[2,5,8,11],[2,5,8,12],[2,5,10,12],[2,6,10,12],[2,6,10,14],[2,6,12,15],[2,7,9,13],[2,7,9,14],[2,7,10,14],[2,8,11,15],[2,8,12,15],[3,4,5,14],[3,4,5,15],[3,4,12,14],[3,5,10,15],[3,5,11,14],[3,7,12,13],[3,11,13,14],[3,12,13,14],[4,5,6,7],[4,5,6,14],[4,5,7,15],[4,6,7,11],[4,6,10,11],[4,6,10,14],[4,7,11,15],[4,8,9,12],[4,8,9,13],[4,8,10,13],[4,8,10,14],[4,8,12,14],[4,10,11,13],[5,6,7,13],[5,7,9,13],[5,7,9,15],[5,8,9,12],[5,8,9,13],[5,9,10,12],[5,9,10,15],[6,7,11,12],[6,7,12,13],[6,10,11,12],[6,12,13,15],[7,8,10,14],[7,8,11,15],[7,8,14,15],[7,9,14,15],[8,12,14,15],[9,10,11,12],[9,10,11,16],[9,10,15,16],[9,11,14,16],[9,14,15,16],[10,11,13,16],[10,13,15,16],[11,13,14,16],[12,13,14,15],[13,14,15,16]])
   pi_1 = fundamental_group(P)
   @test describe(pi_1) == "SL(2,5)"

   @test describe(fundamental_group(torus())) == "Z x Z"

end
