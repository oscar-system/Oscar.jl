@testset "Integer Partitions" begin

  ############################################################################
  # Partition constructor
  ############################################################################
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

  ############################################################################
  # number_of_partitions(n)
  ############################################################################

  # From https://oeis.org/A000041
  @test [ number_of_partitions(i) for i in 0:49 ] ==
    ZZRingElem[ 1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627, 792, 1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, 12310, 14883, 17977, 21637, 26015, 31185, 37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754, 147273, 173525 ]


  # For some random large numbers, checked with Sage
  # Partitions(991).cardinality()
  @test number_of_partitions(991) == ZZ(16839773100833956878604913215477)

  ############################################################################
  # partitions(n)
  ############################################################################
  for n = 0:20
    P = collect(partitions(n))

    # Check that the number of partitions is correct
    # Note that number_of_partitions(n) is computed independently of partitions(n)
    @test length(P) == number_of_partitions(n)

    # Check that all partitions are distinct
    @test P == unique(P)

    # Check that partitions are really partitions of n
    for lambda in P
      @test sum(lambda) == n
    end
  end

  ############################################################################
  # number_of_partitions(n,k)
  ############################################################################
  @test number_of_partitions(0,0) == 1
  @test number_of_partitions(1,0) == 0
  @test number_of_partitions(1,1) == 1
  @test number_of_partitions(0,1) == 0
  @test number_of_partitions(2,3) == 0

  # From https://oeis.org/A008284
  @test [ number_of_partitions(n,k) for n in 1:14 for k in 1:n ] == 
    ZZRingElem[ 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 1, 1, 1, 3, 3, 2, 1, 1, 1, 3, 4, 3, 2, 1, 1, 1, 4, 5, 5, 3, 2, 1, 1, 1, 4, 7, 6, 5, 3, 2, 1, 1, 1, 5, 8, 9, 7, 5, 3, 2, 1, 1, 1, 5, 10, 11, 10, 7, 5, 3, 2, 1, 1, 1, 6, 12, 15, 13, 11, 7, 5, 3, 2, 1, 1, 1, 6, 14, 18, 18, 14, 11, 7, 5, 3, 2, 1, 1, 1, 7, 16, 23, 23, 20, 15, 11, 7, 5, 3, 2, 1, 1 ]

  # For some random large numbers, checked with Sage
  # Partitions(1991,length=170).cardinality()
  @test number_of_partitions(1991,170) == ZZ(22381599503916828837298114953756766080813312)
  @test number_of_partitions(1991,1000) == ZZ(16839773100833956878604913215477)
  @test number_of_partitions(1991,670) == ZZ(3329965216307826492368402165868892548)
  @test number_of_partitions(1991,1991) == ZZ(1)
  @test number_of_partitions(1991,1) == ZZ(1)

  # From Knuth (2011), p. 25.
  @test sum([number_of_partitions(30, i) for i in 0:10]) == 3590

  ############################################################################
  # partitions(n,k)
  ############################################################################
  for n in 0:20
    for k = 0:n+1
      P = collect(partitions(n,k))

      # Create the same by filtering all partitions
      Q = collect(partitions(n))
      filter!( Q->length(Q) == k, Q)

      # Check that P and Q coincide (up to reordering)
      @test length(P) == length(Q)
      @test Set(P) == Set(Q)

      # Compare length with number_of_partitions(n,k)
      @test length(P) == number_of_partitions(n,k)
    end
  end

  ############################################################################
  # partitions(n,k,l1,l2)
  ############################################################################
  for n in 0:20
    for k = 0:n+1
      for l1 = 0:n
        for l2 = l1:n
          P = collect(partitions(n,k,l1,l2))

          # Create the same by filtering all partitions
          Q = collect(partitions(n,k))
          filter!( Q->all(>=(l1),Q), Q)
          filter!( Q->all(<=(l2),Q), Q)

          # Check that P and Q coincide (up to reordering)
          @test length(P) == length(Q)
          @test Set(P) == Set(Q)
        end
      end
    end
  end

  ############################################################################
  # partitions(n,k,l1,l2; only_distinct_parts=true)
  ############################################################################
  for n in 0:20
    for k = 0:n+1
      for l1 = 0:n
        for l2 = l1:n
          P = collect(partitions(n, k, l1, l2; only_distinct_parts=true))

          # Create the same by filtering all partitions
          Q = collect(partitions(n, k, l1, l2))
          filter!( Q->Q==unique(Q), Q )

          # Check that P and Q coincide (up to reordering)
          @test length(P) == length(Q)
          @test Set(P) == Set(Q)
        end
      end
    end
  end

  ############################################################################
  # partitions(m,n,v,mu)
  ############################################################################

  # Check argument errors
  @test_throws ArgumentError partitions(5,2, [1,2], [1,2,3])
  @test_throws ArgumentError partitions(-1,2, [1,2,3], [1,2,3])
  @test_throws ArgumentError partitions(5,0,[1,2,3], [1,2])
  @test_throws ArgumentError partitions(6,3,[3,2,1],[1,2,3]) #req v strictly increasing
  @test_throws ArgumentError partitions(6,3,[3,2,1],[0,2,3]) #mu > 0
  @test_throws ArgumentError partitions(6,3,[0,2,1],[1,2,3]) #v > 0

  # Issues from https://github.com/oscar-system/Oscar.jl/issues/2043
  @test length(collect(partitions(17, 3, [1, 4], [1,4]))) == 0
  @test collect(partitions(17, 5, [1, 4], [1, 4])) == [ partition(4, 4, 4, 4, 1) ]
  @test length(collect(partitions(17, 6, [1, 2], [1, 7]))) == 0
  @test length(collect(partitions(20, 5, [1, 2, 3], [1, 3, 6]))) == 0

  # Issues UT found
  @test length(collect(partitions(1, 1, [1], [1]))) == 1
  @test length(collect(partitions(100, 7, [1, 2, 5, 10, 20, 50], [2, 2, 2, 2, 2, 2]))) == 1

  # Special cases
  for n in 0:20
    for k in 0:n + 1
      P = collect(partitions(n, k, [i for i in 1:n], [n for i in 1:n]))
      Q = collect(partitions(n, k))
      @test length(P) == length(Q)
      @test Set(P) == Set(Q)

      P = collect(partitions(n, k, [i for i in 1:n], [1 for i in 1:n]))
      Q = collect(partitions(n, k, 1, n; only_distinct_parts=true))
      @test length(P) == length(Q)
      @test Set(P) == Set(Q)
    end
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


  ############################################################################
  # dominates(P)
  ############################################################################
  @test dominates(partition([4,2]), partition([3,2,1])) == true
  @test dominates(partition([4,1,1]), partition([3,3])) == false
  @test dominates(partition([3,3]), partition([4,1,1])) == false
  @test dominates(partition([5]), partition([2,2,2])) == false

  # Ful97, page 26
  @test dominates( partition(3,1), partition(2,2) ) == true
  @test dominates( partition(4,1), partition(3,3) ) == false

  ############################################################################
  # conjugate(P)
  ############################################################################

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
