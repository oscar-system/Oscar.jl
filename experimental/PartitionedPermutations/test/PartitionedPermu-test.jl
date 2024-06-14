@testset "PartitionedPermutation" begin
    @testset "PartitionedPermutation Constructor" begin
        @test permutation(partitioned_permutation(perm(symmetric_group(4), [2, 1, 3, 4]), [1, 1, 2, 2])) == perm(symmetric_group(4), [2, 1, 3, 4])
        @test partition(partitioned_permutation(perm(symmetric_group(4), [2, 1, 3, 4]), [1, 1, 2, 2])) == set_partition([1, 1, 2, 2], Int64[])

        @test_throws ArgumentError partitioned_permutation(perm(symmetric_group(3), [2, 1, 3]), [1, 1, 2, 2])
        @test_throws ArgumentError partitioned_permutation(perm(symmetric_group(4), [2, 1, 3, 4]), [1, 2, 3, 3])
    end
    @testset "PartitionedPermutation length" begin
        @test length(partitioned_permutation(perm(symmetric_group(4), [2, 1, 3, 4]), [1, 1, 2, 2])) == 4
        @test adjusted_length(partitioned_permutation(perm(symmetric_group(4), [2, 1, 3, 4]), [1, 1, 2, 2])) == 3
    end
    @testset "EnumeratePartitionedPermutations.jl" begin
        @test length(enumerate_partitioned_permutations(3)) == 13
        @test length(enumerate_partitioned_permutations(5)) == 501
        @test length(enumerate_partitioned_permutations(7)) == 37633
    end
    @testset "PartitionedPermutation Product" begin
        @test partitioned_permutation(perm(symmetric_group(6), [2, 1, 3, 6, 4, 5]), [1, 1, 1, 1, 1, 1]) * partitioned_permutation(perm(symmetric_group(6), [6, 1, 2, 4, 5, 3]), [1, 1, 1, 1, 1, 1]) == 
            partitioned_permutation(perm(symmetric_group(6), [1, 2, 3, 4, 5, 6]), [1, 2, 3, 4, 5, 6])
        @test partitioned_permutation(perm(symmetric_group(3), [1, 2, 3]), [1, 2, 3]) * partitioned_permutation(perm(symmetric_group(3), [2, 1, 3]), [1, 1, 3]) ==
            partitioned_permutation(perm(symmetric_group(3), [2, 1, 3]), [1, 1, 2])
        @test length(factor(partitioned_permutation(perm(symmetric_group(3), [2, 1, 3]), [1, 1, 2]))) == 2
        @test length(factor(partitioned_permutation(perm(symmetric_group(4), [2, 1, 4, 3]), [1, 1, 2, 2]))) <
              length(factor(partitioned_permutation(perm(symmetric_group(4), [1, 2, 3, 4]), [1, 2, 3, 4])))
        @test length(factor(partitioned_permutation(perm(symmetric_group(3), [2, 1, 3]), [1, 1, 2]))) <
              length(factor(partitioned_permutation(perm(symmetric_group(3), [1, 2, 3]), [1, 2, 3])))
    end
end
