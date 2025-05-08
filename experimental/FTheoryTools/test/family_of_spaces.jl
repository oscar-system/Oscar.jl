coord_ring, _ = QQ[:f, :g, :Kbar, :u]
grading = [4 6 1 0; 0 0 0 1]
d = 3
f = family_of_spaces(coord_ring, grading, d)

@testset "Attributes of FamilyOfSpaces" begin
  @test coordinate_ring(f) == coord_ring
  @test weights(f) == grading
  @test dim(f) == d
  @test ngens(irrelevant_ideal(f)) == 4
  @test ngens(ideal_of_linear_relations(f)) == 2
end
