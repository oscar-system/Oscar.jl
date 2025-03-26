@testset "Multipartitions" for T in [Int, ZZRingElem]
  #constructors
  @test multipartition([T[3,2], T[1]]) == multipartition([Partition(T[3,2]), Partition(T[1])])
  @test multipartition([T[3,2], T[]]) == multipartition([Partition(T[3,2]),Partition(T[])])

  # multi-partitions
  check = true
  for n in 0:10
    for r in 1:5
      MP = multipartitions(T(n),r)
      # check that all multipartitions are distinct
      if MP != unique(MP)
        check = false
        break
      end
      # check that multipartititons are really multipartitions of n
      for mp in MP
        if sum(mp) != n
          check = false
          break
        end
      end
      # check that all multipartitions have r parts
      if length(MP) != 0 && unique([ length(mp) for mp in MP ]) != [r]
        check = false
        break
      end
    end
  end
  @test check == true

  #counting
  for n in 0:5
    for k in 1:n+1
      @test length(multipartitions(T(n),k)) == number_of_multipartitions(T(n),k)
    end
  end

end
