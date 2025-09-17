@testset "Compositions with integer type $T" for T in [Int, Int8, ZZRingElem]
  # Check some stupid cases
  @test number_of_compositions(T(0), T(0)) == 1
  @test number_of_compositions(T(0), T(1)) == 0
  @test number_of_compositions(T(1), T(0)) == 0
  @test number_of_compositions(T(0)) == 1

  # First few number of compositions from https://oeis.org/A011782
  nums = [1, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648, 4294967296, 8589934592]

  # Check if number_of_compositions is correct in these examples
  @test [number_of_compositions(T(n)) for n in 0:length(nums) - 1] == nums

  # Complete check of compositions for small cases
  for n in 0:5
    # We'll piece together all compositions of n by the compositions
    # of n into k parts for 0 <= k <= n
    allcomps = []
    for k in 0:n
      C = @inferred collect(compositions(T(n), T(k)))
      @test C isa Vector{Oscar.Composition{T}}
      @test length(C) == number_of_compositions(T(n), T(k))

      # Check if each composition consists of k parts and sums up to n
      for c in C
        @test length(c) == k
        @test sum(c) == n
        @test all(is_positive, c)
      end

      # Check if all compositions are distinct
      @test allunique(C)

      append!(allcomps, C)
    end

    # All compositions need to be distinct
    @test allunique(allcomps)

    # Number of compositions needs to be correct
    @test length(allcomps) == number_of_compositions(T(n))

    # Finally, check compositions(n) function
    allcomps2 = @inferred collect(compositions(T(n)))
    @test allcomps2 isa Vector{Oscar.Composition{T}}
    @test allunique(allcomps2)
    @test Set(allcomps) == Set(allcomps2)
  end

  # Test the case k > n
  C = @inferred collect(compositions(T(2), T(3)))
  @test C isa Vector{Oscar.Composition{T}}
  @test isempty(C)
end

@testset "Ascending compositions with integer type $T" for T in [Int, Int8, ZZRingElem]
  for n in 0:20
    C = @inferred collect(ascending_compositions(T(n)))
    @test length(C) == number_of_partitions(T(n))
    @test allunique(C)
    for lambda in C
      @test sum(lambda) == n
      @test issorted(lambda)
    end
  end
end
