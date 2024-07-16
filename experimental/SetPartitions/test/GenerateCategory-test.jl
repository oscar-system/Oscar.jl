@testset "GenerateCategory test" begin
    @testset "Classical Partitions" begin
        nonc2 = construct_category([set_partition([1], [1]), set_partition([], [1, 1])], 10)
        @test length(nonc2) <= 462
        for partition in nonc2
            @test is_pair(partition) && is_non_crossing(partition) && size(partition) == 10
        end

        nonc = construct_category([set_partition([1], [1]), 
                                    set_partition([], [1, 1]), 
                                    set_partition([1], [1, 1])], 6)
        @test length(nonc) <= 924
        for partition in nonc
            @test is_non_crossing(partition) && size(partition) == 6
        end

        balanced = construct_category([set_partition([1], [1]), 
                                        set_partition([], [1, 1]), 
                                        set_partition([1, 2, 3], [3, 2, 1]), 
                                        set_partition([1, 1], [1, 1])], 6)
        @test length(balanced) <= 217
        for partition in balanced
            if !is_balanced(partition)
                println(partition)
            end
            @test is_balanced(partition) && size(partition) == 6
        end

        nonc12 = construct_category([set_partition([1], [1]), 
                                        set_partition([], [1, 1]), 
                                        set_partition([1], [2])], 6)
        @test length(nonc12) <= 210
        for partition in nonc12
            @test is_non_crossing(partition) && size(partition) == 6
        end
    end

    @testset "Colored Partitions" begin
        nonc2 = construct_category([colored_partition([1], [1], [0], [0]), 
                                    colored_partition([1], [1], [1], [1]), 
                                    colored_partition([], [1, 1], [], [0, 1]), 
                                    colored_partition([], [1, 1], [], [1, 0])]
                                    , 6)

        for partition in nonc2
            @test is_pair(partition) && 
                is_non_crossing(partition) && 
                partition isa ColoredPartition && 
                partition != colored_partition([1], [1], [0], [1]) &&
                size(partition) == 6
        end

        nonc2 = construct_category([colored_partition([1], [1], [0], [0]), 
                                    colored_partition([1], [1], [1], [1]), 
                                    colored_partition([], [1, 1], [], [0, 1]), 
                                    colored_partition([], [1, 1], [], [1, 0]), 
                                    colored_partition([1], [1], [0], [1])], 4)

        for partition in nonc2
            @test is_pair(partition) && 
                is_non_crossing(partition) && 
                partition isa ColoredPartition && 
                size(partition) == 4
        end
    end

    @testset "Spatial Partitions" begin
        P2 = construct_category([spatial_partition([1, 2], [1, 3, 3, 2], 2), 
                                spatial_partition([], [1, 1], 2), 
                                spatial_partition([1, 2], [1, 2], 2), 
                                spatial_partition([], [1, 2, 1, 2], 2)], 
                                6, false, 8)
        @test length(P2) <= 105
        for partition in P2
            @test is_pair(partition) && 
                partition isa SpatialPartition && 
                size(partition) == 6
        end

        nonc2classic = construct_category([set_partition([1], [1]), 
                                            set_partition([], [1, 1])], 8)
        nonc2 = construct_category([spatial_partition([1], [1], 1), 
                                    spatial_partition([], [1, 1], 1)], 8)
        @test length(nonc2) <= 126
        for partition in nonc2
            @test is_pair(partition) && 
                is_non_crossing(set_partition(partition)) && 
                size(partition) == 8 && 
                set_partition(partition) in nonc2classic
        end
        function spatial_rotation_example(p::SpatialPartition)
            if !isempty(upper_points(p))
                return spatial_partition(rotate_top_left(set_partition(p)), levels(p))
            elseif !isempty(lower_points(p))
                return spatial_partition(rotate_bottom_left(set_partition(p)), levels(p))
            else
                return p
            end
        end
        spatial_category_costum_rotation = construct_category([
            spatial_partition([1, 1], [], 1)], 
            2, false, 0, spatial_rotation_example)
        @test spatial_partition([1], [1], 1) in 
            spatial_category_costum_rotation
    end
end
