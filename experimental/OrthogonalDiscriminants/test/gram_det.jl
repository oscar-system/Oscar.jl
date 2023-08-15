@testset "dimensions and Gram determinants of Specht modules" begin
  @test [dimension_specht_module(p) for p in partitions(5)] ==
        [1, 4, 5, 6, 5, 4, 1]
  @test [gram_determinant_specht_module(p) for p in partitions(5)] ==
        [[],
         [[5, 1]],
         [[2, 1], [3, 4]],
         [[2, 6], [5, 3]],
         [[2, 14], [3, 1]],
         [[2, 4], [3, 4], [5, 3]],
         [[2, 3], [3, 1], [5, 1]]]
end
