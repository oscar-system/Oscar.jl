@testset "Functions for SetPartitions" begin
    @test SetPartition([1, 2, 3], [2, 1, 3, 3]) <= SetPartition([1, 1, 2], [1, 1, 2, 2])
    @test !(SetPartition([1, 2, 3], [2, 1, 3, 3]) <= SetPartition([1, 1, 2], [1, 1, 2, 3]))

    @test cycle_partition(Perm([2, 1, 3, 5, 4])) == SetPartition([1, 1, 2, 3, 3], Int64[])
    @test cycle_partition(Perm([2, 4, 6, 1, 3, 5])) == SetPartition([1, 1, 2, 1, 2, 2], Int64[])

    @test join(SetPartition([1, 2, 3], [2, 1, 4, 4]), SetPartition([1, 1, 1], [2, 3, 1, 1])) == SetPartition([1, 1, 1], [1, 1, 1, 1])
    @test_throws ArgumentError join(SetPartition([1, 2], [2, 1, 4, 4]), SetPartition([1, 1, 1], [2, 3, 1, 1]))
    @test_throws ArgumentError join(SetPartition([1, 2], [2, 1, 4, 4]), SetPartition([1, 1], [2, 3, 1]))
end
@testset "PartitionedPermutation" begin
    @testset "PartitionedPermutation Constructor" begin
        @test PartitionedPermutation(Perm([2, 1, 3, 4]), [1, 1, 2, 2]).p == Perm([2, 1, 3, 4])
        @test PartitionedPermutation(Perm([2, 1, 3, 4]), [1, 1, 2, 2]).V == SetPartition([1, 1, 2, 2], Int64[])

        @test_throws ArgumentError PartitionedPermutation(Perm([2, 1, 3]), [1, 1, 2, 2])
        @test_throws ArgumentError PartitionedPermutation(Perm([2, 1, 3, 4]), [1, 2, 3, 3])
    end
    @testset "PartitionedPermutation length" begin
        @test length(PartitionedPermutation(Perm([2, 1, 3, 4]), [1, 1, 2, 2])) == 4
        @test length2(PartitionedPermutation(Perm([2, 1, 3, 4]), [1, 1, 2, 2])) == 3
    end
    @testset "EnumeratePartitionedPermutations.jl" begin
        @test length(enumerate_partitioned_perm(3)) == 13
        @test length(enumerate_partitioned_perm(5)) == 501
        @test length(enumerate_partitioned_perm(7)) == 37633
    end
    @testset "PartitionedPermutation Product" begin
        @test PartitionedPermutation(Perm([2, 1, 3, 6, 4, 5]), [1, 1, 1, 1, 1, 1]) * PartitionedPermutation(Perm([6, 1, 2, 4, 5, 3]), [1, 1, 1, 1, 1, 1]) == 
            PartitionedPermutation(Perm([1, 2, 3, 4, 5, 6]),[1, 2, 3, 4, 5, 6])
        @test PartitionedPermutation(Perm([1, 2, 3]), [1, 2, 3]) * PartitionedPermutation(Perm([2, 1, 3]), [1, 1, 3]) ==
            PartitionedPermutation(Perm([2, 1, 3]), [1, 1, 2])
        @test length(factorization_partitioned_permutation(PartitionedPermutation(Perm([2, 1, 3]), [1, 1, 2]))) == 2
        @test length(factorization_partitioned_permutation(PartitionedPermutation(Perm([2, 1, 4, 3]), [1, 1, 2, 2]))) <
              length(factorization_partitioned_permutation(PartitionedPermutation(Perm([1, 2, 3, 4]), [1, 2, 3, 4])))
        @test length(factorization_partitioned_permutation(PartitionedPermutation(Perm([2, 1, 3]), [1, 1, 2]))) <
              length(factorization_partitioned_permutation(PartitionedPermutation(Perm([1, 2, 3]), [1, 2, 3])))
    end
end
