@testset "Ascending partitions" begin
  for n in 0:20
    # go through all algorithms
    for a in [ :ks, :m ]
      P = Oscar.ascending_partitions(n,algorithm=a)
      # check that number of partitions is correct
      @test length(P) == number_of_partitions(n)
      # check that all partitions are distinct
      @test P == unique(P)
      # check that partitions are really partitions of n
      for lambda in P
        @test sum(lambda) == n
      end
    end
  end
end
