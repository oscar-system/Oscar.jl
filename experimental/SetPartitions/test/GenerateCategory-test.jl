@testset "GenerateCategory test" begin
    @testset "Classical Partitions" begin
        nonc2 = construct_category([Partition([1], [1]), Partition([], [1, 1])], 10)
        @test length(nonc2) <= 462
        for partition in nonc2
            @test is_pair(partition) && is_noncrossing(partition) && PartitionsJulia.size(partition) == 10
        end

        nonc = construct_category([Partition([1], [1]), Partition([], [1, 1]), Partition([1], [1, 1])], 6)
        @test length(nonc) <= 924
        for partition in nonc
            @test is_noncrossing(partition) && PartitionsJulia.size(partition) == 6
        end

        balanced = construct_category([Partition([1], [1]), Partition([], [1, 1]), Partition([1, 2, 3], [3, 2, 1]), Partition([1, 1], [1, 1])], 6)
        @test length(balanced) <= 217
        for partition in balanced
            if !is_balanced(partition)
                println(partition)
            end
            @test is_balanced(partition) && PartitionsJulia.size(partition) == 6
        end

        nonc12 = construct_category([Partition([1], [1]), Partition([], [1, 1]), Partition([1], [2])], 6)
        @test length(nonc12) <= 210
        for partition in nonc12
            @test is_noncrossing(partition) && PartitionsJulia.size(partition) == 6
        end
    end

    @testset "Colored Partitions" begin
        nonc2 = construct_category([ColoredPartition(Partition([1], [1]), [0], [0]), ColoredPartition(Partition([1], [1]), [1], [1]), ColoredPartition(Partition([], [1, 1]), [], [0, 1]), ColoredPartition(Partition([], [1, 1]), [], [1, 0])], 6)

        for partition in nonc2
            @test is_pair(partition) && is_noncrossing(partition) && partition isa ColoredPartition && partition != ColoredPartition(Partition([1], [1]), [0], [1]) && PartitionsJulia.size(partition) == 6
        end

        nonc2 = construct_category([ColoredPartition(Partition([1], [1]), [0], [0]), ColoredPartition(Partition([1], [1]), [1], [1]), ColoredPartition(Partition([], [1, 1]), [], [0, 1]), ColoredPartition(Partition([], [1, 1]), [], [1, 0]), ColoredPartition(Partition([1], [1]), [0], [1])], 4)

        for partition in nonc2
            @test is_pair(partition) && is_noncrossing(partition) && partition isa ColoredPartition && PartitionsJulia.size(partition) == 4
        end
    end

    @testset "Spatial Partitions" begin
        P2 = construct_category([SpatialPartition(Partition([1, 2], [1, 3, 3, 2]), 2), SpatialPartition(Partition([], [1, 1]), 2), SpatialPartition(Partition([1, 2], [1, 2]), 2), SpatialPartition(Partition([], [1, 2, 1, 2]), 2)], 6, false, 8)
        @test length(P2) <= 105
        for partition in P2
            @test is_pair(partition) && partition isa SpatialPartition && PartitionsJulia.size(partition) == 6
        end

        nonc2classic = construct_category([Partition([1], [1]), Partition([], [1, 1])], 8)
        nonc2 = construct_category([SpatialPartition(Partition([1], [1]), 1), SpatialPartition(Partition([], [1, 1]), 1)], 8)
        @test length(nonc2) <= 126
        for partition in nonc2
            @test is_pair(partition) && is_noncrossing(partition.partition) && PartitionsJulia.size(partition) == 8 && partition.partition in nonc2classic
        end
    end
end