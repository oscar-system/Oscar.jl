@testset "compositions" begin
  # Check some stupid cases
  @test number_of_compositions(0, 0) == 1
  @test number_of_compositions(0, 1) == 0
  @test number_of_compositions(1, 0) == 0
  @test number_of_compositions(0) == 1

  # First few number of compositions from https://oeis.org/A011782
  nums = [1, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648, 4294967296, 8589934592]

  # Check if number_of_compositions is correct in these examples
  @test [number_of_compositions(n) for n in 0:length(nums) - 1] == nums

  # Complete check of compositions for small cases
  for n in 0:5
    # We'll piece together all compositions of n by the compositions
    # of n into k parts for 0 <= k <= n
    allcomps = []
    for k in 0:n
      C = collect(compositions(n, k))
      @test length(C) == number_of_compositions(n, k)

      # Check if each composition consists of k parts and sums up to n
      for c in C
        @test length(c) == k
        @test sum(c) == n
        @test all(is_positive, c)
      end

      # Check if all compositions are distinct
      @test C == unique(C)

      append!(allcomps, C)
    end

    # All compositions need to be distinct
    @test allcomps == unique(allcomps)

    # Number of compositions needs to be correct
    @test length(allcomps) == number_of_compositions(n)

    # Finally, check compositions(n) function
    allcomps2 = collect(compositions(n))
    @test allcomps2 == unique(allcomps2)
    @test Set(allcomps) == Set(allcomps2)
  end
end

@testset "Ascending compositions" begin
  for n in 0:20
    C = collect(ascending_compositions(n))
    @test length(C) == number_of_partitions(n)
    @test C == unique(C)
    for lambda in C
      @test sum(lambda) == n
      @test all(i -> lambda[i] <= lambda[i + 1], 1:length(lambda) - 1)
    end
  end
end
