@testset "Multipartitions for integer type $T" for T in [Int, ZZRingElem]
  #constructors
  @test multipartition([T[3,2], T[1]]) == multipartition([Partition(T[3,2]), Partition(T[1])])
  @test multipartition([T[3,2], T[]]) == multipartition([Partition(T[3,2]),Partition(T[])])

  # multi-partitions
  @testset "multipartitions($n, $r)" for n in 0:10, r in 1:5
    MP = multipartitions(T(n),r)

    @test MP == unique(MP)
    @test all(mp -> length(mp) == r, MP)
    for mp in MP
      @test sum(mp) == n
    end
  end

  #counting
  @testset "Counting multipartitions($n, $r)" for n in 0:5, r in 1:n+1
      @test length(multipartitions(T(n),r)) == number_of_multipartitions(T(n),r)
  end

end
