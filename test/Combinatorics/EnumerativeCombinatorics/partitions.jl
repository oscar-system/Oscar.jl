@testset "Integer Partitions" begin

  ############################################################################
  # Partition constructor
  ############################################################################
  @testset "Partition constructor" begin
    @test partition(Int8,2,2,1) == partition(Int8[2,2,1])
    @test typeof(partition(Int8,2,2,1)) == typeof(partition(Int8[2,2,1]))
    @test partition(Int8,1) == partition(Int8[1])
    @test typeof(partition(Int8,1)) == typeof(partition(Int8[1]))
    @test partition(Int8,) == partition(Int8[])
    @test typeof(partition(Int8,)) == typeof(partition(Int8[]))

    @test partition(2,2,1) == partition([2,2,1])
    @test partition(1) == partition([1])
    @test partition() == partition([])
    @test partition(Int8) == partition(Int8[])
  end

  ############################################################################
  # number_of_partitions(n)
  ############################################################################
  @testset "number_of_partitions(n) for integer type $T" for T in [Int, Int8, ZZRingElem]
    # From https://oeis.org/A000041
    @test [ number_of_partitions(T(i)) for i in 0:49 ] ==
      ZZRingElem[ 1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627, 792, 1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, 12310, 14883, 17977, 21637, 26015, 31185, 37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754, 147273, 173525 ]


    # For some random large numbers, checked with Sage
    # Partitions(991).cardinality()
    if T !== Int8
      @test number_of_partitions(T(991)) == ZZ(16839773100833956878604913215477)
    end

    @test number_of_partitions(T(-1)) == ZZ(0)
  end

  ############################################################################
  # partitions(n)
  ############################################################################
  @testset "partitions($n) for integer type $T" for n in 0:10, T in [Int, Int8, ZZRingElem]
    P = @inferred collect(partitions(T(n)))
    @test P isa Vector{Oscar.Partition{T}}

    # Check that the number of partitions is correct
    # Note that number_of_partitions(n) is computed independently of partitions(n)
    @test length(P) == number_of_partitions(T(n))

    # Check that all partitions are distinct
    @test allunique(P)

    # Check that partitions are really partitions of n
    for lambda in P
      @test sum(lambda) == n
    end
  end

  ############################################################################
  # number_of_partitions(n,k)
  ############################################################################
  @testset "number_of_partitions(n,k) for integer type $T" for T in [Int, Int8, ZZRingElem]
    @test number_of_partitions(T(0), T(0)) == 1
    @test number_of_partitions(T(1), T(0)) == 0
    @test number_of_partitions(T(1), T(1)) == 1
    @test number_of_partitions(T(0), T(1)) == 0
    @test number_of_partitions(T(2), T(3)) == 0

    # From https://oeis.org/A008284
    @test [ number_of_partitions(T(n), T(k)) for n in 1:14 for k in 1:n ] ==
      ZZRingElem[ 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 1, 1, 1, 3, 3, 2, 1, 1, 1, 3, 4, 3, 2, 1, 1, 1, 4, 5, 5, 3, 2, 1, 1, 1, 4, 7, 6, 5, 3, 2, 1, 1, 1, 5, 8, 9, 7, 5, 3, 2, 1, 1, 1, 5, 10, 11, 10, 7, 5, 3, 2, 1, 1, 1, 6, 12, 15, 13, 11, 7, 5, 3, 2, 1, 1, 1, 6, 14, 18, 18, 14, 11, 7, 5, 3, 2, 1, 1, 1, 7, 16, 23, 23, 20, 15, 11, 7, 5, 3, 2, 1, 1 ]

    # For some random large numbers, checked with Sage
    # Partitions(1991,length=170).cardinality()
    if T !== Int8
      @test number_of_partitions(T(1991), T(170)) == ZZ(22381599503916828837298114953756766080813312)
      @test number_of_partitions(T(1991), T(1000)) == ZZ(16839773100833956878604913215477)
      @test number_of_partitions(T(1991), T(670)) == ZZ(3329965216307826492368402165868892548)
      @test number_of_partitions(T(1991), T(1991)) == ZZ(1)
      @test number_of_partitions(T(1991), T(1)) == ZZ(1)
    end

    # From Knuth (2011), p. 25.
    @test sum([number_of_partitions(T(30), T(i)) for i in 0:10]) == 3590

    @test number_of_partitions(T(-1), T(0)) == ZZ(0)
    @test number_of_partitions(T(0), T(-1)) == ZZ(0)
  end

  ############################################################################
  # partitions(n,k)
  ############################################################################
  @testset "partitions($n, $k) for integer type $T" for n in 0:10, k in 0:n+1, T in [Int, Int8, ZZRingElem]
    P = @inferred collect(partitions(T(n), T(k)))
    @test P isa Vector{Oscar.Partition{T}}

    # Create the same by filtering all partitions
    Q = @inferred collect(partitions(T(n)))
    @test Q isa Vector{Oscar.Partition{T}}
    filter!(Q -> length(Q) == k, Q)

    # Check that P and Q coincide (up to reordering)
    @test length(P) == length(Q)
    @test Set(P) == Set(Q)

    # Compare length with number_of_partitions(n,k)
    @test length(P) == number_of_partitions(T(n), T(k))
  end

  ############################################################################
  # partitions(n,k,lb,ub)
  ############################################################################
  @testset "partitions($n, $k, lb, ub)" for n in 0:10, k in 0:n+1
    @testset "partitions($n, $k, $lb, $ub) for integer type $T" for lb in 0:n, ub in lb:n, T in [Int, Int8, ZZRingElem]
      P = @inferred collect(partitions(T(n), T(k), T(lb), T(ub)))
      @test P isa Vector{Oscar.Partition{T}}

      # Create the same by filtering all partitions
      Q = @inferred collect(partitions(T(n), T(k)))
      @test Q isa Vector{Oscar.Partition{T}}
      filter!( Q->all(>=(lb),Q), Q)
      filter!( Q->all(<=(ub),Q), Q)

      # Check that P and Q coincide (up to reordering)
      @test length(P) == length(Q)
      @test Set(P) == Set(Q)
    end
  end

  ############################################################################
  # partitions(n,k,lb,ub; only_distinct_parts=true)
  ############################################################################
  @testset "partitions($n, $k, lb, ub; only_distinct_parts=true)" for n in 0:10, k in 0:n+1
    @testset "partitions($n, $k, $lb, $ub; only_distinct_parts=true) for integer type $T" for lb in 0:n, ub in lb:n, T in [Int, Int16, ZZRingElem]
      P = @inferred collect(partitions(T(n), T(k), T(lb), T(ub); only_distinct_parts=true))
      @test P isa Vector{Oscar.Partition{T}}

      # Create the same by filtering all partitions
      Q = @inferred collect(partitions(T(n), T(k), T(lb), T(ub)))
      @test Q isa Vector{Oscar.Partition{T}}
      filter!(allunique, Q)

      # Check that P and Q coincide (up to reordering)
      @test length(P) == length(Q)
      @test issetequal(P, Q)
    end
  end

  ############################################################################
  # partitions(n,k,v,mu)
  ############################################################################
  @testset "partitions(n, k, v, mu) for integer type $T" for T in [Int, Int16, ZZRingElem]
    # Check argument errors
    @test_throws ArgumentError partitions(T(5), 2, T[1,2], [1,2,3])
    @test_throws ArgumentError partitions(T(-1), 2, T[1,2,3], [1,2,3])
    @test_throws ArgumentError partitions(T(5), 0, T[1,2,3], [1,2])
    @test_throws ArgumentError partitions(T(6), 3, T[3,2,1], [1,2,3]) #req v strictly increasing
    @test_throws ArgumentError partitions(T(6), 3, T[3,2,1], [0,2,3]) #mu > 0
    @test_throws ArgumentError partitions(T(6), 3, T[0,2,1], [1,2,3]) #v > 0

    # Issues from https://github.com/oscar-system/Oscar.jl/issues/2043
    @test length(collect(partitions(T(17), 3, T[1, 4], [1,4]))) == 0
    @test collect(partitions(T(17), 5, T[1, 4], [1, 4])) == [ partition(T[4, 4, 4, 4, 1]) ]
    @test length(collect(partitions(T(17), 6, T[1, 2], [1, 7]))) == 0
    @test length(collect(partitions(T(20), 5, T[1, 2, 3], [1, 3, 6]))) == 0

    # Issues UT found
    @test length(collect(partitions(T(1), 1, T[1], [1]))) == 1
    @test length(collect(partitions(T(100), 7, T[1, 2, 5, 10, 20, 50], [2, 2, 2, 2, 2, 2]))) == 1

    # Special cases
    @testset "Special cases n=$n, k=$k" for n in 0:20, k in 0:n+1
      P = @inferred collect(partitions(T(n), k, T[i for i in 1:n], [n for i in 1:n]))
      @test P isa Vector{Oscar.Partition{T}}
      Q = collect(partitions(T(n), T(k)))
      @test length(P) == length(Q)
      @test Set(P) == Set(Q)

      P = @inferred collect(partitions(T(n), k, T[i for i in 1:n], [1 for i in 1:n]))
      @test P isa Vector{Oscar.Partition{T}}
      Q = collect(partitions(T(n), T(k), T(1), T(n); only_distinct_parts=true))
      @test length(P) == length(Q)
      @test Set(P) == Set(Q)
    end

    # From https://www.maa.org/frank-morgans-math-chat-293-ways-to-make-change-for-a-dollar
    @test length(collect(partitions(100, [1, 5, 10, 25, 50]))) == 292
    @test length(collect(partitions(200, [1, 5, 10, 25, 50, 100]))) == 2728

    # From Knu11, Exercise 11 on page 408
    @test length(collect(partitions(100, [1, 2, 5, 10, 20, 50], [2, 2, 2, 2, 2, 2]))) == 6
    @test length(collect(partitions(100, [1, 2, 5, 10, 20, 50]))) == 4562

    # From https://oeis.org/A000008
    @test [ length(collect(partitions(n, [1,2,5,10]))) for n in 0:60 ] ==
      [ 1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 11, 12, 15, 16, 19, 22, 25, 28, 31, 34,
      40, 43, 49, 52, 58, 64, 70, 76, 82, 88, 98, 104, 114, 120, 130, 140,
      150, 160, 170, 180, 195, 205, 220, 230, 245, 260, 275, 290, 305, 320,
      341, 356, 377, 392, 413, 434, 455, 476, 497, 518, 546 ]

    @test collect(partitions(0, Int[])) == [partition(Int[])]
    @test collect(partitions(1, Int[])) == Partition{Int}[]
    @test collect(partitions(0, Int[], Int[])) == [partition(Int[])]
    @test collect(partitions(1, Int[], Int[])) == Partition{Int}[]
  end

  ############################################################################
  # dominates(P)
  ############################################################################
  @testset "dominates(P)" begin
    @test dominates(partition([4,2]), partition([3,2,1])) == true
    @test dominates(partition([4,1,1]), partition([3,3])) == false
    @test dominates(partition([3,3]), partition([4,1,1])) == false
    @test dominates(partition([5]), partition([2,2,2])) == false

    # Ful97, page 26
    @test dominates( partition(3,1), partition(2,2) ) == true
    @test dominates( partition(4,1), partition(3,3) ) == false
  end

  ############################################################################
  # conjugate(P)
  ############################################################################
  @testset "conjugate(P)" begin
    # From Knu11, page 394
    @test conjugate(partition(8,8,8,7,2,1,1)) == partition([7, 5, 4, 4, 4, 4, 4, 3])

    # From Ful97, page 2
    @test conjugate(partition([6,4,3,1])) == partition([4, 3, 3, 2, 1, 1])

    @test conjugate(partition([])) == partition([])

    # Check if conjugate is an involution
    for p in partitions(10)
      @test conjugate(conjugate(p)) == p
    end
  end
end
