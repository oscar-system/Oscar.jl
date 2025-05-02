@testset "SetPartition.jl test" begin
    @testset "Partition constructors" begin
        @test upper_points(set_partition([2, 3], [8, 9])) == [1, 2] 
        @test lower_points(set_partition([2, 3], [8, 9])) == [3, 4]
        
        @test upper_points(colored_partition([2, 4], [8, 9], [1, 0], [0, 1])) == [1, 2]
        @test lower_points(colored_partition([2, 4], [8, 9], [1, 0], [0, 1])) == [3, 4]
        @test_throws ArgumentError colored_partition([1, 2], [3], [1, 0], [])
        @test_throws ArgumentError colored_partition([1, 2], [5], [2, 0], [1])

        @test upper_points(spatial_partition([89, 99], [1, 2], 2)) == [1, 2]
        @test lower_points(spatial_partition([89, 99], [1, 2], 2)) == [3, 4]
        @test_throws ArgumentError spatial_partition([1, 2, 3], [2, 3], 2)
        @test_throws ArgumentError spatial_partition([1, 2], [2, 1], -1)
    end
    @testset "Operations on classical partitions" begin
        @test involution(set_partition([2, 3], [4, 4])) == set_partition([1, 1], [2, 3])
        
        @test tensor_product(set_partition([1, 2, 3], [2, 1, 2]), 
                    set_partition([3, 3], [])) == set_partition([1, 2, 3, 4, 4], [2, 1, 2])
        
        @test_throws ArgumentError compose(set_partition([1, 2, 3], [3]), 
                    set_partition([1, 2], [1, 2]))
        
        @test compose(set_partition([1, 2, 2, 2], [2, 3, 4]), 
                    set_partition([1, 2], [3, 3, 1, 2])) == set_partition([1, 1], [1, 2, 3])
        
        @test compose(set_partition([1, 2, 3, 4, 1, 3], [3, 4]), 
                    set_partition([1, 2, 3], [1, 1, 2, 2, 3, 1])) == 
                    set_partition([1, 1, 1], [1, 1])
        
        @test_throws ArgumentError rotate_top_left(set_partition([], [1])) 
        
        @test rotate_bottom_left(set_partition([1, 2], [3, 2])) == 
                set_partition([1, 2, 3], [3])
        @test rotate_top_left(set_partition([1, 2], [3, 2])) == 
                set_partition([1], [2, 3, 1])
        @test rotate_top_right(set_partition([1, 2], [3, 2])) == 
                set_partition([1], [2, 3, 3])
        @test rotate_bottom_right(set_partition([1, 2], [3, 2])) == 
                set_partition([1, 2, 2], [3])
        
        @test reflect_vertical(set_partition([2, 3, 2, 2], [2, 3])) == 
                    set_partition([1, 1, 2, 1], [2, 1])
    end
    @testset "Operations on colored partitions" begin
        @test involution(colored_partition([99, 1, 4], [4, 99], [0, 1, 1], [1, 0])) == 
                colored_partition([1, 2], [2, 3, 1], [1, 0], [0, 1, 1])

        @test tensor_product(colored_partition([1, 2], [2, 1], [0, 0], [0, 0]), 
                             colored_partition([1], [1, 1],[1], [0, 0])) == 
                colored_partition([1, 2, 3], [2, 1, 3, 3], [0, 0, 1], [0, 0, 0, 0])

        @test compose(colored_partition([1, 2], [2, 1], [0, 1], [1, 0]), 
                      colored_partition([1, 2], [2, 1], [1, 0], [0, 1])) == 
                colored_partition([1, 2], [1, 2], [1, 0], [1, 0])

        @test_throws ArgumentError compose(colored_partition([1, 2], [2, 1], [1, 1], [1, 0]), 
                                           colored_partition([1, 2], [2, 1], [1, 0], [0, 1]))

        @test rotate_bottom_left(colored_partition([1, 2], [3, 2], [1, 1], [1, 1])) == 
                colored_partition([1, 2, 3], [3], [0, 1, 1], [1])
        @test rotate_top_left(colored_partition([1, 2], [3, 2], [1, 1], [1, 1])) == 
                colored_partition([1], [2, 3, 1], [1], [0, 1, 1])
        @test rotate_top_right(colored_partition([1, 2], [3, 2], [1, 1], [1, 1])) == 
                colored_partition([1], [2, 3, 3], [1], [1, 1, 0])
        @test rotate_bottom_right(colored_partition([1, 2], [3, 2], [1, 1], [1, 1])) == 
                colored_partition([1, 2, 2], [3], [1, 1, 0], [1])
    end
    @testset "Operations on spatial partitions" begin
        @test involution(spatial_partition([99, 1, 4, 4], [4, 99], 2)) == 
                         spatial_partition([1, 2], [2, 3, 1, 1], 2)

        @test tensor_product(spatial_partition([1, 2], [2, 1], 1), 
                             spatial_partition([1], [1, 1], 1)) == 
                spatial_partition([1, 2, 3], [2, 1, 3, 3], 1)
                            
        @test_throws ArgumentError tensor_product(spatial_partition([1, 2], [2, 1], 2), 
                                                  spatial_partition([1], [1, 1], 1))

        @test compose(spatial_partition([1, 2], [2, 1], 2), 
                      spatial_partition([1, 2], [2, 1], 2)) == 
                spatial_partition([1, 2], [1, 2], 2)

        @test_throws ArgumentError compose(spatial_partition([1, 2], [2, 1], 2), 
                                           spatial_partition([1, 2], [2, 1], 1))
    end
    @testset "Functions for Partitioned Permutations" begin
        @test is_dominated_by(set_partition([1, 2, 3], [2, 1, 3, 3]), set_partition([1, 1, 2], [1, 1, 2, 2]))
        @test !is_dominated_by(set_partition([1, 2, 3], [2, 1, 3, 3]), set_partition([1, 1, 2], [1, 1, 2, 3]))
    
        @test cycle_partition(perm(symmetric_group(5), [2, 1, 3, 5, 4])) == set_partition([1, 1, 2, 3, 3], Int64[])
        @test cycle_partition(perm(symmetric_group(6), [2, 4, 6, 1, 3, 5])) == set_partition([1, 1, 2, 1, 2, 2], Int64[])
    
        @test join(set_partition([1, 2, 3], [2, 1, 4, 4]), set_partition([1, 1, 1], [2, 3, 1, 1])) == set_partition([1, 1, 1], [1, 1, 1, 1])
        @test_throws ArgumentError join(set_partition([1, 2], [2, 1, 4, 4]), set_partition([1, 1, 1], [2, 3, 1, 1]))
        @test_throws ArgumentError join(set_partition([1, 2], [2, 1, 4, 4]), set_partition([1, 1], [2, 3, 1]))
    end
end
