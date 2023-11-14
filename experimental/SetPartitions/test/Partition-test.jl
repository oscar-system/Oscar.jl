@testset "SetPartition.jl test" begin
    @testset "Partition constructors" begin
        @test SetPartition([2, 3], [8, 9]).upper_points == [1, 2] && SetPartition([2, 3], [8, 9]).lower_points == [3, 4]
        @test ColoredPartition(SetPartition([2, 4], [8, 9]), [1, 0], [0, 1]).partition.upper_points == [1, 2] && ColoredPartition(SetPartition([2, 4], [8, 9]), [1, 0], [0, 1]).partition.lower_points == [3, 4]
        @test SpatialPartition(SetPartition([89, 99], [1, 2]), 2).partition.upper_points == [1, 2] && SpatialPartition(SetPartition([89, 99], [1, 2]), 2).partition.lower_points == [3, 4]
    end
    @testset "Operations on classical partitions" begin
        @test involution(SetPartition([2, 3], [4, 4])) == SetPartition([1, 1], [2, 3])
        @test tensor_product(SetPartition([1, 2, 3], [2, 1, 2]), SetPartition([3, 3], [])) == SetPartition([1, 2, 3, 4, 4], [2, 1, 2])
        @test_throws ErrorException composition(SetPartition([1, 2, 3], [3]), SetPartition([1, 2], [1, 2]))
        @test composition(SetPartition([1, 2, 2, 2], [2, 3, 4]), SetPartition([1, 2], [3, 3, 1, 2])) == SetPartition([1, 1], [1, 2, 3])
        @test composition(SetPartition([1, 2, 3, 4, 1, 3], [3, 4]), SetPartition([1, 2, 3], [1, 1, 2, 2, 3, 1])) == SetPartition([1, 1, 1], [1, 1])
        @test_throws ErrorException rotation(SetPartition([], [1]), true, true) 
        @test rotation(SetPartition([1, 2], [2, 1]), false, false) == SetPartition([1, 2, 1], [2]) 
        @test vertical_reflection(SetPartition([2, 3, 2, 2], [2, 3])) == SetPartition([1, 1, 2, 1], [2, 1])
    end
    @testset "Operations on colored partitions" begin
        @test involution(ColoredPartition(SetPartition([99, 1, 4], [4, 99]), [0, 1], [1, 0])) == ColoredPartition(SetPartition([1, 2], [2, 3, 1]), [1, 0], [0, 1])
        @test tensor_product(ColoredPartition(SetPartition([1, 2], [2, 1]), [0, 0], [0, 0]), ColoredPartition(SetPartition([1], [1, 1]), [1], [0, 0])) == ColoredPartition(SetPartition([1, 2, 3], [2, 1, 3, 3]), [0, 0, 1], [0, 0, 0, 0])
        @test composition(ColoredPartition(SetPartition([1, 2], [2, 1]), [0, 1], [1, 0]), ColoredPartition(SetPartition([1, 2], [2, 1]), [1, 0], [0, 1])) == ColoredPartition(SetPartition([1, 2], [1, 2]), [1, 0], [1, 0])
        @test_throws ErrorException composition(ColoredPartition(SetPartition([1, 2], [2, 1]), [1, 1], [1, 0]), ColoredPartition(SetPartition([1, 2], [2, 1]), [1, 0], [0, 1]))
        @test rotation(ColoredPartition(SetPartition([1, 2], [3, 2]), [1, 1], [1, 1]), true, false) == ColoredPartition(SetPartition([1, 2, 3], [3]), [0, 1, 1], [1])
    end
    @testset "Operations on spatial partitions" begin
        @test involution(SpatialPartition(SetPartition([99, 1, 4, 4], [4, 99]), 2)) == SpatialPartition(SetPartition([1, 2], [2, 3, 1, 1]), 2)
        @test tensor_product(SpatialPartition(SetPartition([1, 2], [2, 1]), 1), SpatialPartition(SetPartition([1], [1, 1]), 1)) == SpatialPartition(SetPartition([1, 2, 3], [2, 1, 3, 3]), 1)
        @test_throws ErrorException tensor_product(SpatialPartition(SetPartition([1, 2], [2, 1]), 2), SpatialPartition(SetPartition([1], [1, 1]), 1))
        @test composition(SpatialPartition(SetPartition([1, 2], [2, 1]), 2), SpatialPartition(SetPartition([1, 2], [2, 1]), 2)) == SpatialPartition(SetPartition([1, 2], [1, 2]), 2)
        @test_throws ErrorException composition(SpatialPartition(SetPartition([1, 2], [2, 1]), 2), SpatialPartition(SetPartition([1, 2], [2, 1]), 1))
    end
end
