@testset "Young Tableaux" begin
  # reading_word
  @test reading_word(young_tableau([ [1,2,5,7], [3,4], [6]])) == [6,3,4,1,2,5,7]
  @test reading_word(young_tableau([ [1], [2], [3]])) == [3,2,1]
  @test reading_word(young_tableau([[1,2,3]])) == [1,2,3]
  @test reading_word(young_tableau(Array{Int,1}[])) == Int[]

  # weight
  @test weight(young_tableau([[1,2,3],[2,3],[3]])) == [1,2,3]
  @test weight(young_tableau([[1,2,3,4,5]])) == [1,1,1,1,1]
  @test_throws  ArgumentError  weight(young_tableau([[1],[1],[1]]))
  @test weight(young_tableau(Array{Int,1}[])) == Int[]

  # is_standard
  @test is_standard(young_tableau([[1,2,4,7,8],[3,5,6,9],[10]])) == true
  @test is_standard(young_tableau([[1,2],[3,4]])) == true
  @test is_standard(young_tableau([[1,3],[2,4]])) == true
  @test is_standard(young_tableau([[1,4],[2,4]])) == false
  @test is_standard(young_tableau([[1,2],[4]])) == false
  @test is_standard(young_tableau([[1,3,2],[4]])) == false
  @test is_standard(young_tableau([[-1]])) == false

  # is_semistandard
  @test is_semistandard(young_tableau([[1,2,4,7,8],[3,5,6,9],[10]])) == true
  @test is_semistandard(young_tableau([[1,2],[3,4]])) == true
  @test is_semistandard(young_tableau([[1,3],[2,4]])) == true
  @test is_semistandard(young_tableau([[1,4],[2,4]])) == false
  @test is_semistandard(young_tableau([[1,2],[4]])) == true
  @test is_semistandard(young_tableau([[1,2,2],[3]])) == true
  @test is_semistandard(young_tableau([[1,2,3],[1,4]])) == false
  @test is_semistandard(young_tableau([[1,2,1],[2,4]])) == false
  @test is_semistandard(young_tableau([[-1]])) == true

  @testset "Generating tableaux with integer type $T" for T in [Int8, Int]
    # semistandard_tableaux(shape::Array{T,1}, max_val=sum(shape)::Integer)
    shapes = [T[3, 2, 1], T[3, 3, 1], T[2, 2, 2]]
    for s in shapes
      SST = @inferred collect(semistandard_tableaux(s))
      @test SST isa Vector{Oscar.YoungTableau{T}}
      #check that all tableaux are distinct
      @test allunique(SST)

      #check that all tableaux are semistandard_tableaux
      for tab in SST
        @test is_semistandard(tab)
      end
    end
    @test isempty(semistandard_tableaux(T[3, 2, 1], T(2)))

    # semistandard_tableaux(s::Array{T,1}, weight::Array{T,1})
    shapes = [T[5, 3, 1, 1], T[4, 3, 2, 1], T[2, 2, 2, 2, 2]]
    weights = [T[1, 1, 1, 1, 1, 1, 1, 1, 1, 1], T[3, 0, 2, 0, 0, 5], T[4, 3, 2, 1]]
    for s in shapes
      for w in weights
        SST = @inferred collect(semistandard_tableaux(s, w))
        @test SST isa Vector{Oscar.YoungTableau{T}}
        #check that all tableaux are distinct
        @test allunique(SST)
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
    @test collect(semistandard_tableaux(T[], T[])) == [young_tableau(Vector{T}[])]
    @test length(collect(semistandard_tableaux(partition(T[5, 3, 1, 1]), T[4, 3, 2, 1]))) == 2

    #semistandard_tableaux(box_num, max_val)
    BoxNum = T(0):T(5)
    MaxVal = T(1):T(6)
    for box_num in BoxNum
      for max_val in MaxVal
        SST = @inferred collect(semistandard_tableaux(box_num, max_val))
        @test SST isa Vector{Oscar.YoungTableau{T}}
        #check that all tableaux are distinct
        @test allunique(SST)
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
      for s in partitions(T(i))
        ST = @inferred collect(standard_tableaux(s))
        @test ST isa Vector{Oscar.YoungTableau{T}}
        #check that all tableaux are distinct
        @test allunique(ST)
        #check that all tableaux are standard_tableaux
        for tab in ST
          @test is_standard(tab)
        end
        #check that all tableaux where found
        @test length(ST) == number_of_standard_tableaux(s)
      end
    end
    @test collect(standard_tableaux(partition(T[]))) == [young_tableau(Vector{T}[])]
    @test collect(standard_tableaux(T[3, 2, 1])) == collect(standard_tableaux(partition(T[3, 2, 1])))

    # standard_tableaux(n::Integer)
    for n = 0:10
      ST = @inferred collect(standard_tableaux(T(n)))
      @test ST isa Vector{Oscar.YoungTableau{T}}
      #check that all tableaux are distinct
      @test allunique(ST)
      #check that all tableaux are standard_tableaux
      for tab in ST
        @test is_standard(tab)
      end
      #check that all tableaux have n boxes
      for tab in ST
        @test sum(shape(tab)) == n
      end
    end
  end # testest "Generating tableaux"

  # hook_length
  @test hook_length(partition([1]),1,1) == 1
  @test hook_length(partition([4,3,1,1]),1,1) == 7
  @test hook_length(young_tableau([[1,2,3,4],[5,6,7],[8],[9]]),1,1) == 7

  # hook_lengths
  @test hook_lengths(partition([4,3,1,1])) == young_tableau([[7,4,3,1],[5,2,1],[2],[1]])
  @test hook_lengths(partition([1])) == young_tableau([[1]])
  @test hook_lengths(partition([])) == young_tableau(Array{Int,1}[])

  # schensted
  @test schensted([6,2,7,3,5,4,1]) == (young_tableau([[1,3,4],[2,7],[5],[6]]),young_tableau([[1,3,5],[2,4],[6],[7]]))
  @test schensted([5,2,7,1,3,8,6,4]) == (young_tableau([[1,3,4],[2,6,8],[5,7]]),young_tableau([[1,3,6],[2,5,7],[4,8]]))
  @test schensted([1]) == (young_tableau([[1]]),young_tableau([[1]]))
  @test schensted(Int[]) == (young_tableau(Array{Int,1}[]),young_tableau(Array{Int,1}[]))

  # bump!
  tab = young_tableau(Array{Int,1}[])
  tab2 = young_tableau(Array{Int,1}[])
  Q = young_tableau(Array{Int,1}[])
  for x in [1,2,1,1,3,4,1,1]
    bump!(tab,x)
    bump!(tab2, x, Q, x)
  end
  @test tab == young_tableau([[1,1,1,1,1],[2,3,4]])
  @test tab2 == young_tableau([[1,1,1,1,1],[2,3,4]])

end


Oscar.@_AuxDocTest "Print Young Tableaux", (fix = false),
raw"""
    some (potentially mathematically meaningless) "tableau"

```jldoctest young_tableau_printing.test
julia> using Oscar

julia> young_tableau(Vector{Int}[])
Empty Young tableau

julia> young_tableau([[1, 2, 3, 4], [5, 6], [7], [8]])
+---+---+---+---+
| 1 | 2 | 3 | 4 |
+---+---+---+---+
| 5 | 6 |
+---+---+
| 7 |
+---+
| 8 |
+---+

julia> young_tableau([[1, 2, 3, 4], [5, 6], [7, 8, 9]], check = false)
+---+---+---+---+
| 1 | 2 | 3 | 4 |
+---+---+---+---+
| 5 | 6 |
+---+---+---+
| 7 | 8 | 9 |
+---+---+---+

julia> young_tableau([Int[], [1], Int[], Int[], [1, 2, 3], Int[]], check = false)
+
|
+---+
| 1 |
+---+
|
|
|
+---+---+---+
| 1 | 2 | 3 |
+---+---+---+
|
+

julia> young_tableau([[1, 2, 3, 4, 5, 6, 7, 8, 9], [10]])
+----+----+----+----+----+----+----+----+----+
|  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 |
+----+----+----+----+----+----+----+----+----+
| 10 |
+----+

julia> old_unicode = Oscar.allow_unicode(true; temporary=true);

julia> young_tableau(Vector{Int}[])
Empty Young tableau

julia> young_tableau([[1, 2, 3, 4], [5, 6], [7], [8]])
┌───┬───┬───┬───┐
│ 1 │ 2 │ 3 │ 4 │
├───┼───┼───┴───┘
│ 5 │ 6 │
├───┼───┘
│ 7 │
├───┤
│ 8 │
└───┘

julia> young_tableau([[1, 2, 3, 4], [5, 6], [7, 8, 9]], check = false)
┌───┬───┬───┬───┐
│ 1 │ 2 │ 3 │ 4 │
├───┼───┼───┴───┘
│ 5 │ 6 │
├───┼───┼───┐
│ 7 │ 8 │ 9 │
└───┴───┴───┘

julia> young_tableau([Int[], [1], Int[], Int[], [1, 2, 3], Int[]], check = false)
┌
│
├───┐
│ 1 │
├───┘
│
│
│
├───┬───┬───┐
│ 1 │ 2 │ 3 │
├───┴───┴───┘
│
└

julia> young_tableau([[1, 2, 3, 4, 5, 6, 7, 8, 9], [10]])
┌────┬────┬────┬────┬────┬────┬────┬────┬────┐
│  1 │  2 │  3 │  4 │  5 │  6 │  7 │  8 │  9 │
├────┼────┴────┴────┴────┴────┴────┴────┴────┘
│ 10 │
└────┘

julia> Oscar.allow_unicode(old_unicode; temporary=true);
```
"""
