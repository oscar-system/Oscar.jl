@testset "PartitionProperties test" begin
    @testset "pair" begin
        @test is_pair(set_partition([1, 2, 2, 1], [3, 3]))
        @test !is_pair(set_partition([1, 2, 2, 1], [4, 3]))
        @test !is_pair(set_partition([1, 1, 1], [3, 3, 3, 3]))
        @test is_pair(set_partition([2, 1], [1, 2]))
    end
    @testset "balanced" begin
        @test is_balanced(set_partition([1, 2, 2, 1], [3, 3]))
        @test is_balanced(set_partition([1, 2, 3], [3, 2, 1]))
        @test !is_balanced(set_partition([1, 1, 3], [3, 3]))
    end
    @testset "non-crossing" begin
        @test is_non_crossing(set_partition([1, 2, 3, 4, 5, 5, 4, 3, 2, 1], [1, 6, 6, 1]))
        @test !is_non_crossing(set_partition([1, 2, 3, 4, 5, 5, 4, 3, 2, 1], [1, 6, 6, 2]))
    end
end
