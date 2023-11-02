@testset "Partition.jl test" begin
    @testset "Partition constructors" begin
        @test Partition([2, 3], [8, 9]).upper_points == [1, 2] && Partition([2, 3], [8, 9]).lower_points == [3, 4]
        @test ColoredPartition(Partition([2, 4], [8, 9]), [1, 0], [0, 1]).partition.upper_points == [1, 2] && ColoredPartition(Partition([2, 4], [8, 9]), [1, 0], [0, 1]).partition.lower_points == [3, 4]
        @test SpatialPartition(Partition([89, 99], [1, 2]), 2).partition.upper_points == [1, 2] && SpatialPartition(Partition([89, 99], [1, 2]), 2).partition.lower_points == [3, 4]
    end
    @testset "Operations on classical partitions" begin
        @test involution(Partition([2, 3], [4, 4])) == Partition([1, 1], [2, 3])
        @test tensor_product(Partition([1, 2, 3], [2, 1, 2]), Partition([3, 3], [])) == Partition([1, 2, 3, 4, 4], [2, 1, 2])
        @test_throws ErrorException composition(Partition([1, 2, 3], [3]), Partition([1, 2], [1, 2]))
        @test composition(Partition([1, 2, 2, 2], [2, 3, 4]), Partition([1, 2], [3, 3, 1, 2])) == Partition([1, 1], [1, 2, 3])
        @test composition(Partition([1, 2, 3, 4, 1, 3], [3, 4]), Partition([1, 2, 3], [1, 1, 2, 2, 3, 1])) == Partition([1, 1, 1], [1, 1])
        @test_throws ErrorException rotation(Partition([], [1]), true, true) 
        @test rotation(Partition([1, 2], [2, 1]), false, false) == Partition([1, 2, 1], [2]) 
        @test vertical_reflection(Partition([2, 3, 2, 2], [2, 3])) == Partition([1, 1, 2, 1], [2, 1])
    end
    @testset "Operations on colored partitions" begin
        @test involution(ColoredPartition(Partition([99, 1, 4], [4, 99]), [0, 1], [1, 0])) == ColoredPartition(Partition([1, 2], [2, 3, 1]), [1, 0], [0, 1])
        @test tensor_product(ColoredPartition(Partition([1, 2], [2, 1]), [0, 0], [0, 0]), ColoredPartition(Partition([1], [1, 1]), [1], [0, 0])) == ColoredPartition(Partition([1, 2, 3], [2, 1, 3, 3]), [0, 0, 1], [0, 0, 0, 0])
        @test composition(ColoredPartition(Partition([1, 2], [2, 1]), [0, 1], [1, 0]), ColoredPartition(Partition([1, 2], [2, 1]), [1, 0], [0, 1])) == ColoredPartition(Partition([1, 2], [1, 2]), [1, 0], [1, 0])
        @test_throws ErrorException composition(ColoredPartition(Partition([1, 2], [2, 1]), [1, 1], [1, 0]), ColoredPartition(Partition([1, 2], [2, 1]), [1, 0], [0, 1]))
        @test rotation(ColoredPartition(Partition([1, 2], [3, 2]), [1, 1], [1, 1]), true, false) == ColoredPartition(Partition([1, 2, 3], [3]), [0, 1, 1], [1])
    end
    @testset "Operations on spatial partitions" begin
        @test involution(SpatialPartition(Partition([99, 1, 4, 4], [4, 99]), 2)) == SpatialPartition(Partition([1, 2], [2, 3, 1, 1]), 2)
        @test tensor_product(SpatialPartition(Partition([1, 2], [2, 1]), 1), SpatialPartition(Partition([1], [1, 1]), 1)) == SpatialPartition(Partition([1, 2, 3], [2, 1, 3, 3]), 1)
        @test_throws ErrorException tensor_product(SpatialPartition(Partition([1, 2], [2, 1]), 2), SpatialPartition(Partition([1], [1, 1]), 1))
        @test composition(SpatialPartition(Partition([1, 2], [2, 1]), 2), SpatialPartition(Partition([1, 2], [2, 1]), 2)) == SpatialPartition(Partition([1, 2], [1, 2]), 2)
        @test_throws ErrorException composition(SpatialPartition(Partition([1, 2], [2, 1]), 2), SpatialPartition(Partition([1, 2], [2, 1]), 1))
    end
end