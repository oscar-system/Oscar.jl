@testset "Young Tableaux" begin
  # reading_word
  @test reading_word(YoungTableau([ [1,2,5,7], [3,4], [6]])) == [6,3,4,1,2,5,7]
  @test reading_word(YoungTableau([ [1], [2], [3]])) == [3,2,1]
  @test reading_word(YoungTableau([[1,2,3]])) == [1,2,3]
  @test reading_word(YoungTableau(Array{Int,1}[])) == Int[]

  # weight
  @test weight(YoungTableau([[1,2,3],[1,2],[1]])) == [3,2,1]
  @test weight(YoungTableau([[1,2,3,4,5]])) == [1,1,1,1,1]
  @test weight(YoungTableau([[1],[1],[1]])) == [3]
  @test weight(YoungTableau(Array{Int,1}[])) == Int[]

  # is_standard
  @test is_standard(YoungTableau([[1,2,4,7,8],[3,5,6,9],[10]])) == true
  @test is_standard(YoungTableau([[1,2],[3,4]])) == true
  @test is_standard(YoungTableau([[1,3],[2,4]])) == true
  @test is_standard(YoungTableau([[1,4],[2,4]])) == false
  @test is_standard(YoungTableau([[1,2],[4]])) == false
  @test is_standard(YoungTableau([[1,3,2],[4]])) == false

  # is_semistandard
  @test is_semistandard(YoungTableau([[1,2,4,7,8],[3,5,6,9],[10]])) == true
  @test is_semistandard(YoungTableau([[1,2],[3,4]])) == true
  @test is_semistandard(YoungTableau([[1,3],[2,4]])) == true
  @test is_semistandard(YoungTableau([[1,4],[2,4]])) == false
  @test is_semistandard(YoungTableau([[1,2],[4]])) == true
  @test is_semistandard(YoungTableau([[1,2,2],[3]])) == true
  @test is_semistandard(YoungTableau([[1,2,3],[1,4]])) == false
  @test is_semistandard(YoungTableau([[1,2,1],[2,4]])) == false

  # semistandard_tableaux(shape::Array{T,1}, max_val=sum(shape)::Integer)
  shapes = [[3,2,1],[3,3,1],[2,2,2]]
  for s in shapes
    SST = semistandard_tableaux(s)
    #check that all tableaux are distinct
    @test SST == unique(SST)

    #check that all tableaux are semistandard_tableaux
    for tab in SST
      @test is_semistandard(tab)
    end
  end
  @test isempty(semistandard_tableaux([3,2,1],2))

  # semistandard_tableaux(s::Array{T,1}, weight::Array{T,1})
  shapes = [[5,3,1,1],[4,3,2,1],[2,2,2,2,2]]
  weights = [[1,1,1,1,1,1,1,1,1,1],[3,0,2,0,0,5],[4,3,2,1]]
  for s in shapes
    for w in weights
      SST = semistandard_tableaux(s,w)
      #check that all tableaux are distinct
      @test SST == unique(SST)
      #check that all tableaux are semistandard_tableaux
      for tab in SST
        @test is_semistandard(tab)
      end
      #check that all tableaux have the correct shape
      for tab in SST
        @test shape(tab) == s
      end
      #check that all tableaux have the correct weight
      for tab in SST
        @test weight(tab) == w
      end
    end
  end
  @test semistandard_tableaux(Int[], Int[]) == [YoungTableau(Array{Int,1}[])]

  #semistandard_tableaux(box_num, max_val)
  BoxNum = 0:5
  MaxVal = 1:6
  for box_num in BoxNum
    for max_val in MaxVal
      SST = semistandard_tableaux(box_num, max_val)
      #check that all tableaux are distinct
      @test SST == unique(SST)
      #check that all tableaux are semistandard_tableaux
      for tab in SST
        @test is_semistandard(tab)
      end
      #check that all tableaux have box_num boxes
      for tab in SST
        @test sum(shape(tab)) == box_num
      end
      #check that all tableaux have values ≤ max_val
      for tab in SST
        for i in 1:length(tab)
          @test tab[i][end] <= max_val
        end
      end
    end
  end

  # number_of_standard_tableaux
  # standard_tableaux(s::Partition)
  for i = 1:10
    for s in partitions(i)
      ST = standard_tableaux(s)
      #check that all tableaux are distinct
      @test ST == unique(ST)
      #check that all tableaux are standard_tableaux
      for tab in ST
        @test is_standard(tab)
      end
      #check that all tableaux where found
      @test length(ST) == number_of_standard_tableaux(s)
    end
  end
  @test standard_tableaux(partition(Int[])) == [YoungTableau(Array{Int,1}[])]
  @test standard_tableaux([3,2,1]) == standard_tableaux(partition([3,2,1]))

  # standard_tableaux(n::Integer)
  for n = 0:10
    ST = standard_tableaux(n)
    #check that all tableaux are distinct
    @test ST == unique(ST)
    #check that all tableaux are standard_tableaux
    for tab in ST
      @test is_standard(tab)
    end
    #check that all tableaux have n boxes
    for tab in ST
      @test sum(shape(tab)) == n
    end
  end

  # hook_length
  @test hook_length(partition([1]),1,1) == 1
  @test hook_length(partition([4,3,1,1]),1,1) == 7
  @test hook_length(YoungTableau([[1,2,3,4],[5,6,7],[8],[9]]),1,1) == 7

  # hook_lengths
  @test hook_lengths(partition([4,3,1,1])) == YoungTableau([[7,4,3,1],[5,2,1],[2],[1]])
  @test hook_lengths(partition([1])) == YoungTableau([[1]])
  @test hook_lengths(partition([])) == YoungTableau(Array{Int,1}[])

  # schensted
  @test schensted([6,2,7,3,5,4,1]) == (YoungTableau([[1,3,4],[2,7],[5],[6]]),YoungTableau([[1,3,5],[2,4],[6],[7]]))
  @test schensted([5,2,7,1,3,8,6,4]) == (YoungTableau([[1,3,4],[2,6,8],[5,7]]),YoungTableau([[1,3,6],[2,5,7],[4,8]]))
  @test schensted([1]) == (YoungTableau([[1]]),YoungTableau([[1]]))
  @test schensted(Int[]) == (YoungTableau(Array{Int,1}[]),YoungTableau(Array{Int,1}[]))

  # bump!
  tab = YoungTableau(Array{Int,1}[])
  tab2 = YoungTableau(Array{Int,1}[])
  Q = YoungTableau(Array{Int,1}[])
  for x in [1,2,1,1,3,4,1,1]
    bump!(tab,x)
    bump!(tab2, x, Q, x)
  end
  @test tab == YoungTableau([[1,1,1,1,1],[2,3,4]])
  @test tab2 == YoungTableau([[1,1,1,1,1],[2,3,4]])

end
