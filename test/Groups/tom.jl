# Make sure that the tests start without cached tables of marks
# and without cached character tables.
# (Running the tests will store information that changes some test outputs,
# thus running the tests twice needs these calls.)
GAP.Globals.UnloadTableOfMarksData()
empty!(Oscar.tables_of_marks_by_id)
GAP.Globals.UnloadCharacterTableData()
empty!(Oscar.character_tables_by_id)

Oscar.@_AuxDocTest "show and print tables of marks", (fix = false),
raw"""
```jldoctest tables_of_marks.test
julia> using Oscar

julia> t_a4 = table_of_marks(alternating_group(4));

julia> t_a5 = table_of_marks("A5");
```

`print` shows an abbrev. form

```jldoctest tables_of_marks.test
julia> print(t_a4)
table of marks of Alt(4)

julia> print(t_a5)
table of marks of A5
```

`show` uses the abbrev. form for nested objects

```jldoctest tables_of_marks.test
julia> show([t_a4])
Oscar.GAPGroupTableOfMarks[table of marks of Alt(4)]

julia> show([t_a5])
Oscar.GAPGroupTableOfMarks[table of marks of A5]
```

terse printing
```jldoctest tables_of_marks.test
julia> print(AbstractAlgebra.terse(stdout), t_a4)
table of marks of a group

julia> print(AbstractAlgebra.terse(stdout), t_a5)
table of marks of a group
```

default `show`

```jldoctest tables_of_marks.test
julia> show(stdout, MIME("text/plain"), t_a4)
Table of marks of alternating group of degree 4 and order 12

1: 12        
2:  6 2      
3:  4 . 1    
4:  3 3 . 3  
5:  1 1 1 1 1
```

LaTeX format

```jldoctest tables_of_marks.test
julia> show(stdout, MIME("text/latex"), t_a4)
Table of marks of alternating group of degree 4 and order 12

$\begin{array}{rrrrrr}
1: & 12 &  &  &  &  \\ 
2: & 6 & 2 &  &  &  \\ 
3: & 4 & . & 1 &  &  \\ 
4: & 3 & 3 & . & 3 &  \\ 
5: & 1 & 1 & 1 & 1 & 1 \\
\end{array}
$
```

Show some separating lines, in the screen format ...

```jldoctest tables_of_marks.test
julia> show(IOContext(stdout, :separators_col => [0,5],
                                 :separators_row => [0,5]),
                   MIME("text/plain"), t_a5)
A5

--+----------+-------
1:|60        |       
2:|30 2      |       
3:|20 . 2    |       
4:|15 3 . 3  |       
5:|12 . . . 2|       
--+----------+-------
6:|10 2 1 . .|1      
7:| 6 2 . . 1|. 1    
8:| 5 1 2 1 .|. . 1  
9:| 1 1 1 1 1|1 1 1 1
```

... and in LaTeX format

```jldoctest tables_of_marks.test
julia> show(IOContext(stdout, :separators_col => [0,5],
                                 :separators_row => [0,5]),
                   MIME("text/latex"), t_a5)
A5

$\begin{array}{r|rrrrr|rrrr}
\hline
1: & 60 &  &  &  &  &  &  &  &  \\ 
2: & 30 & 2 &  &  &  &  &  &  &  \\ 
3: & 20 & . & 2 &  &  &  &  &  &  \\ 
4: & 15 & 3 & . & 3 &  &  &  &  &  \\ 
5: & 12 & . & . & . & 2 &  &  &  &  \\ 
\hline
6: & 10 & 2 & 1 & . & . & 1 &  &  &  \\ 
7: & 6 & 2 & . & . & 1 & . & 1 &  &  \\ 
8: & 5 & 1 & 2 & 1 & . & . & . & 1 &  \\ 
9: & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\
\end{array}
$
```

distribute the table into column portions, in the screen format ...

```jldoctest tables_of_marks.test
julia> show(IOContext(stdout, :separators_col => [0],
                            :separators_row => [0],
                            :portions_col => [4,5]),
              MIME("text/plain"), t_a5)
A5

--+--------
1:|60      
2:|30 2    
3:|20 . 2  
4:|15 3 . 3
5:|12 . . .
6:|10 2 1 .
7:| 6 2 . .
8:| 5 1 2 1
9:| 1 1 1 1

--+---------
1:|         
2:|         
3:|         
4:|         
5:|2        
6:|. 1      
7:|1 . 1    
8:|. . . 1  
9:|1 1 1 1 1
```

... and in LaTeX format

```jldoctest tables_of_marks.test
julia> show(IOContext(stdout, :separators_col => [0],
                            :separators_row => [0],
                            :portions_col => [4,5]),
              MIME("text/latex"), t_a5)
A5

$\begin{array}{r|rrrr}
\hline
1: & 60 &  &  &  \\ 
2: & 30 & 2 &  &  \\ 
3: & 20 & . & 2 &  \\ 
4: & 15 & 3 & . & 3 \\ 
5: & 12 & . & . & . \\ 
6: & 10 & 2 & 1 & . \\ 
7: & 6 & 2 & . & . \\ 
8: & 5 & 1 & 2 & 1 \\ 
9: & 1 & 1 & 1 & 1 \\
\end{array}

\begin{array}{r|rrrrr}
\hline
1: &  &  &  &  &  \\ 
2: &  &  &  &  &  \\ 
3: &  &  &  &  &  \\ 
4: &  &  &  &  &  \\ 
5: & 2 &  &  &  &  \\ 
6: & . & 1 &  &  &  \\ 
7: & 1 & . & 1 &  &  \\ 
8: & . & . & . & 1 &  \\ 
9: & 1 & 1 & 1 & 1 & 1 \\
\end{array}
$
```

distribute the table into row portions,
in the screen format (perhaps not relevant) ...

```jldoctest tables_of_marks.test
julia> show(IOContext(stdout, :separators_col => [0],
                                 :separators_row => [0],
                                 :portions_row => [4,5]),
                   MIME("text/plain"), t_a5)
A5

--+------------------
1:|60                
2:|30 2              
3:|20 . 2            
4:|15 3 . 3          

5:|12 . . . 2        
6:|10 2 1 . . 1      
7:| 6 2 . . 1 . 1    
8:| 5 1 2 1 . . . 1  
9:| 1 1 1 1 1 1 1 1 1
```

... and in LaTeX format (may be interesting)

```jldoctest tables_of_marks.test
julia> show(IOContext(stdout, :separators_col => [0],
                                 :separators_row => [0],
                                 :portions_row => [4,5]),
                   MIME("text/latex"), t_a5)
A5

$\begin{array}{r|rrrrrrrrr}
\hline
1: & 60 &  &  &  &  &  &  &  &  \\ 
2: & 30 & 2 &  &  &  &  &  &  &  \\ 
3: & 20 & . & 2 &  &  &  &  &  &  \\ 
4: & 15 & 3 & . & 3 &  &  &  &  &  \\
\end{array}

\begin{array}{r|rrrrrrrrr}
5: & 12 & . & . & . & 2 &  &  &  &  \\ 
6: & 10 & 2 & 1 & . & . & 1 &  &  &  \\ 
7: & 6 & 2 & . & . & 1 & . & 1 &  &  \\ 
8: & 5 & 1 & 2 & 1 & . & . & . & 1 &  \\ 
9: & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\
\end{array}
$
```

Test the case where a group has a custom name where the first character
should not be turned into lowercase
```jldoctest tables_of_marks.test
julia> table_of_marks(SL(2,2))
Table of marks of SL(2,2)

1: 6      
2: 3 1    
3: 2 . 2  
4: 1 1 1 1
```
"""


@testset "create tables of marks" begin
  @testset "library of tables of marks" begin
    @test table_of_marks("J5") === nothing
  end
end

@testset "access fields in tables of marks" begin
  # table without group
  t = table_of_marks("A5")
  @test GapObj(t) === t.GAPTable

  # table with `GAPGroup` group
  g = symmetric_group(4)
  t = table_of_marks(g)
  @test GapObj(t) === t.GAPTable
  @test group(t) === t.group === g
  @test Oscar.isomorphism_to_GAP_group(t) === t.isomorphism

  # table with `FinGenAbGroup` group
  g = abelian_group([2, 4])
  t = table_of_marks(g)
  @test GapObj(t) === t.GAPTable
  @test group(t) === t.group === g
  @test Oscar.isomorphism_to_GAP_group(t) === t.isomorphism

  d = Dict(t => 0)
  @test t in keys(d)
end

@testset "attributes of tables of marks" begin
  t = table_of_marks("A5")
  @test identifier(t) == "A5"
  @test order(t) == 60
  m = matrix(t)
  @test [m[9,i] for i in 1:9] == collect(t[end])
  orders = orders_class_representatives(t)
  @test orders == [1, 2, 3, 4, 5, 6, 10, 12, 60]
  @test class_lengths(t) == [1, 15, 10, 5, 6, 10, 6, 5, 1]

  @test all(i -> order(representative(t, i)) == orders[i],
            1:length(t))
  @test_throws ArgumentError representative(t, length(t)+1)
end

@testset "marks vectors" begin
  g = symmetric_group(4)
  t = table_of_marks(g)
  @test nrows(t) == 11
  @test ncols(t) == 11
  @test_throws ErrorException t[12]
  tr = t[end]
  @test parent(tr) === t
  @test tr in t
  @test tr != t[1]
  @test tr == t[end]
  @test coordinates(sum(t)) == ones(ZZRingElem, length(t))
  chi = t[2]
  @test chi isa Oscar.GAPGroupMarksVector
  @test group(chi) === g
  @test GapObj([chi], recursive = true)[1] == GapObj(values(chi))
  @test chi[4] == t[2,4]
  @test [chi[i] for i in 1:2] == values(chi)
  @test [2*chi[i] for i in 1:2] == values(2 * chi)
  @test [2*chi[i] for i in 1:2] == values(chi + chi)
  @test [chi[i]^2 for i in 1:2] == values(chi * chi)
  @test chi * t[end] == chi
  @test [chi[i]^2 for i in 1:2] == values(chi^2)
  @test [-chi[i] for i in 1:2] == values(-chi)
  @test [] == values(chi-chi)
  @test values(t[end]) == values(one(chi))
  @test [] == values(zero(chi))
  @test length(chi) == ncols(t)
  @test 4 in chi
  res = coordinates(chi)
  @test res == [0, 1]
  @test res isa Vector{ZZRingElem}
  @test coordinates(Int, chi) isa Vector{Int}
  @test allunique(t)

  d = Dict(chi => 0)
  @test chi in keys(d)
end

@testset "restriction to the character table" begin
  # compute from prescribed group
  # - start with char. table
  g = symmetric_group(4)
  tom = table_of_marks(g)
  tbl = character_table(g)
  tbl2 = character_table(tom)
  @test tbl === tbl2
  @test table_of_marks(tbl) === tom
  @test restrict(tom[end], tbl) == trivial_character(tbl)

  # - start with table of marks
  g = symmetric_group(3)
  tbl = character_table(g)
  tom = table_of_marks(g)
  tom2 = table_of_marks(tbl)
  @test tom === tom2
  @test character_table(tom) === tbl
  @test restrict(tom[end], tbl) == trivial_character(tbl)

  # fetch from library
  # - start with char. table
  tom = table_of_marks("A5")
  tbl = character_table("A5")
  tbl2 = character_table(tom)
  @test tbl === tbl2
  @test table_of_marks(tbl) === tom
  @test restrict(tom[end], tbl) == trivial_character(tbl)

  # - start with table of marks
  tbl = character_table("A6")
  tom = table_of_marks("A6")
  tom2 = table_of_marks(tbl)
  @test tom === tom2
  @test character_table(tom) === tbl
  @test restrict(tom[end], tbl) == trivial_character(tbl)

  # restrict to a Brauer table
  modtbl = mod(tbl, 2)
  @test restrict(tom[end], modtbl) == trivial_character(modtbl)
end
