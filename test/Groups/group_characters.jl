@testset "show and print character tables" begin
  io = IOBuffer();

  t_a4 = character_table(alternating_group(4))
  t_a5 = character_table("A5")

  # `print` shows an abbrev. form
  print(io, t_a4)
  @test String(take!(io)) == "character_table(Alt( [ 1 .. 4 ] ))"

  # default `show`
  Oscar.with_unicode() do
    show(io, t_a4)
  end
  @test String(take!(io)) ==
  """
  Alt( [ 1 .. 4 ] )
  
   2  2  2       .       .
   3  1  .       1       1
                          
     1a 2a      3a      3b
  2P 1a 1a      3b      3a
  3P 1a 2a      1a      1a
                          
  χ₁  1  1       1       1
  χ₂  1  1 -ζ₃ - 1      ζ₃
  χ₃  1  1      ζ₃ -ζ₃ - 1
  χ₄  3 -1       .       .
  """

  show(io, t_a4)
  @test String(take!(io)) ==
  """
  Alt( [ 1 .. 4 ] )
  
    2  2  2        .        .
    3  1  .        1        1
                             
      1a 2a       3a       3b
   2P 1a 1a       3b       3a
   3P 1a 2a       1a       1a
                             
  X_1  1  1        1        1
  X_2  1  1 -z_3 - 1      z_3
  X_3  1  1      z_3 -z_3 - 1
  X_4  3 -1        .        .
  """

  # LaTeX format
  show(io, MIME("text/latex"), t_a4)
  @test String(take!(io)) ==
  """
  \$Alt( [ 1 .. 4 ] )

  \\begin{array}{rrrrr}
  2 & 2 & 2 & . & . \\\\
  3 & 1 & . & 1 & 1 \\\\
   &  &  &  &  \\\\
   & 1a & 2a & 3a & 3b \\\\
  2P & 1a & 1a & 3b & 3a \\\\
  3P & 1a & 2a & 1a & 1a \\\\
   &  &  &  &  \\\\
  \\chi_{1} & 1 & 1 & 1 & 1 \\\\
  \\chi_{2} & 1 & 1 & -\\zeta_{3} - 1 & \\zeta_{3} \\\\
  \\chi_{3} & 1 & 1 & \\zeta_{3} & -\\zeta_{3} - 1 \\\\
  \\chi_{4} & 3 & -1 & . & . \\\\
  \\end{array}
  \$"""

  # show a legend of irrationalities instead of self-explanatory values,
  # in the screen format ...
  Oscar.with_unicode() do
    show(IOContext(io, :with_legend => true), t_a4)
  end
  @test String(take!(io)) ==
  """
  Alt( [ 1 .. 4 ] )

   2  2  2  .  .
   3  1  .  1  1
                
     1a 2a 3a 3b
  2P 1a 1a 3b 3a
  3P 1a 2a 1a 1a
                
  χ₁  1  1  1  1
  χ₂  1  1  A  A̅
  χ₃  1  1  A̅  A
  χ₄  3 -1  .  .

  A = -ζ₃ - 1
  A̅ = ζ₃
  """

  show(IOContext(io, :with_legend => true), t_a4)
  @test String(take!(io)) ==
  """
  Alt( [ 1 .. 4 ] )
  
    2  2  2  .  .
    3  1  .  1  1
                 
      1a 2a 3a 3b
   2P 1a 1a 3b 3a
   3P 1a 2a 1a 1a
                 
  X_1  1  1  1  1
  X_2  1  1  A /A
  X_3  1  1 /A  A
  X_4  3 -1  .  .
  
  A = -z_3 - 1
  /A = z_3
  """

  # ... and in LaTeX format
  show(IOContext(io, :with_legend => true), MIME("text/latex"), t_a4)
  @test String(take!(io)) ==
  """
  \$Alt( [ 1 .. 4 ] )

  \\begin{array}{rrrrr}
  2 & 2 & 2 & . & . \\\\
  3 & 1 & . & 1 & 1 \\\\
   &  &  &  &  \\\\
   & 1a & 2a & 3a & 3b \\\\
  2P & 1a & 1a & 3b & 3a \\\\
  3P & 1a & 2a & 1a & 1a \\\\
   &  &  &  &  \\\\
  \\chi_{1} & 1 & 1 & 1 & 1 \\\\
  \\chi_{2} & 1 & 1 & A & \\overline{A} \\\\
  \\chi_{3} & 1 & 1 & \\overline{A} & A \\\\
  \\chi_{4} & 3 & -1 & . & . \\\\
  \\end{array}

  A = -\\zeta_{3} - 1
  \\overline{A} = \\zeta_{3}
  \$"""

  # show the screen format for a table with real and non-real irrationalities
  Oscar.with_unicode() do
    show(IOContext(io, :with_legend => true), character_table("L2(11)"))
  end
  @test String(take!(io)) ==
  """
  L2(11)

    2  2  2  1  .  .  1   .   .
    3  1  1  1  .  .  1   .   .
    5  1  .  .  1  1  .   .   .
   11  1  .  .  .  .  .   1   1
                               
      1a 2a 3a 5a 5b 6a 11a 11b
   2P 1a 1a 3a 5b 5a 3a 11b 11a
   3P 1a 2a 1a 5b 5a 2a 11a 11b
   5P 1a 2a 3a 1a 1a 6a 11a 11b
  11P 1a 2a 3a 5a 5b 6a  1a  1a
                               
   χ₁  1  1  1  1  1  1   1   1
   χ₂  5  1 -1  .  .  1   B   B̅
   χ₃  5  1 -1  .  .  1   B̅   B
   χ₄ 10 -2  1  .  .  1  -1  -1
   χ₅ 10  2  1  .  . -1  -1  -1
   χ₆ 11 -1 -1  1  1 -1   .   .
   χ₇ 12  .  .  A A*  .   1   1
   χ₈ 12  .  . A*  A  .   1   1

  A = -ζ₅³ - ζ₅² - 1
  A* = ζ₅³ + ζ₅²
  B = ζ₁₁⁹ + ζ₁₁⁵ + ζ₁₁⁴ + ζ₁₁³ + ζ₁₁
  B̅ = -ζ₁₁⁹ - ζ₁₁⁵ - ζ₁₁⁴ - ζ₁₁³ - ζ₁₁ - 1
  """

  show(IOContext(io, :with_legend => true), character_table("L2(11)"))
  @test String(take!(io)) ==
  """
  L2(11)
  
    2  2  2  1  .  .  1   .   .
    3  1  1  1  .  .  1   .   .
    5  1  .  .  1  1  .   .   .
   11  1  .  .  .  .  .   1   1
                               
      1a 2a 3a 5a 5b 6a 11a 11b
   2P 1a 1a 3a 5b 5a 3a 11b 11a
   3P 1a 2a 1a 5b 5a 2a 11a 11b
   5P 1a 2a 3a 1a 1a 6a 11a 11b
  11P 1a 2a 3a 5a 5b 6a  1a  1a
                               
  X_1  1  1  1  1  1  1   1   1
  X_2  5  1 -1  .  .  1   B  /B
  X_3  5  1 -1  .  .  1  /B   B
  X_4 10 -2  1  .  .  1  -1  -1
  X_5 10  2  1  .  . -1  -1  -1
  X_6 11 -1 -1  1  1 -1   .   .
  X_7 12  .  .  A A*  .   1   1
  X_8 12  .  . A*  A  .   1   1
  
  A = -z_5^3 - z_5^2 - 1
  A* = z_5^3 + z_5^2
  B = z_11^9 + z_11^5 + z_11^4 + z_11^3 + z_11
  /B = -z_11^9 - z_11^5 - z_11^4 - z_11^3 - z_11 - 1
  """

  # show some separating lines, in the screen format ...
  Oscar.with_unicode() do
    show(IOContext(io, :separators_col => [0,5],
                       :separators_row => [0,5]), t_a5)
  end
  @test String(take!(io)) ==
  """
  A5

   2│ 2  2  .             .             .│
   3│ 1  .  1             .             .│
   5│ 1  .  .             1             1│
    │                                    │
    │1a 2a 3a            5a            5b│
  2P│1a 1a 3a            5b            5a│
  3P│1a 2a 1a            5b            5a│
  5P│1a 2a 3a            1a            1a│
    │                                    │
  ──┼────────────────────────────────────┼
  χ₁│ 1  1  1             1             1│
  χ₂│ 3 -1  . ζ₅³ + ζ₅² + 1    -ζ₅³ - ζ₅²│
  χ₃│ 3 -1  .    -ζ₅³ - ζ₅² ζ₅³ + ζ₅² + 1│
  χ₄│ 4  .  1            -1            -1│
  χ₅│ 5  1 -1             .             .│
  ──┼────────────────────────────────────┼
  """

  show(IOContext(io, :separators_col => [0,5],
                     :separators_row => [0,5]), t_a5)
  @test String(take!(io)) ==
  """
  A5
  
    2| 2  2  .                 .                 .|
    3| 1  .  1                 .                 .|
    5| 1  .  .                 1                 1|
     |                                            |
     |1a 2a 3a                5a                5b|
   2P|1a 1a 3a                5b                5a|
   3P|1a 2a 1a                5b                5a|
   5P|1a 2a 3a                1a                1a|
     |                                            |
  ---+--------------------------------------------+
  X_1| 1  1  1                 1                 1|
  X_2| 3 -1  . z_5^3 + z_5^2 + 1    -z_5^3 - z_5^2|
  X_3| 3 -1  .    -z_5^3 - z_5^2 z_5^3 + z_5^2 + 1|
  X_4| 4  .  1                -1                -1|
  X_5| 5  1 -1                 .                 .|
  ---+--------------------------------------------+
  """

  # ... and in LaTeX format
  show(IOContext(io, :separators_col => [0,5],
                     :separators_row => [0,5]),
                     MIME("text/latex"), t_a5)
  @test String(take!(io)) ==
  """
  \$A5

  \\begin{array}{r|rrrrr|}
  2 & 2 & 2 & . & . & . \\\\
  3 & 1 & . & 1 & . & . \\\\
  5 & 1 & . & . & 1 & 1 \\\\
   &  &  &  &  &  \\\\
   & 1a & 2a & 3a & 5a & 5b \\\\
  2P & 1a & 1a & 3a & 5b & 5a \\\\
  3P & 1a & 2a & 1a & 5b & 5a \\\\
  5P & 1a & 2a & 3a & 1a & 1a \\\\
   &  &  &  &  &  \\\\
  \\hline
  \\chi_{1} & 1 & 1 & 1 & 1 & 1 \\\\
  \\chi_{2} & 3 & -1 & . & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 & -\\zeta_{5}^{3} - \\zeta_{5}^{2} \\\\
  \\chi_{3} & 3 & -1 & . & -\\zeta_{5}^{3} - \\zeta_{5}^{2} & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 \\\\
  \\chi_{4} & 4 & . & 1 & -1 & -1 \\\\
  \\chi_{5} & 5 & 1 & -1 & . & . \\\\
  \\hline
  \\end{array}
  \$"""

  # distribute the table into column portions, in the screen format ...
  Oscar.with_unicode() do
    show(IOContext(io, :separators_col => [0],
                       :separators_row => [0],
                       :portions_col => [2,3]), t_a5)
  end
  @test String(take!(io)) ==
  """
  A5

   2│ 2  2
   3│ 1  .
   5│ 1  .
    │     
    │1a 2a
  2P│1a 1a
  3P│1a 2a
  5P│1a 2a
    │     
  ──┼─────
  χ₁│ 1  1
  χ₂│ 3 -1
  χ₃│ 3 -1
  χ₄│ 4  .
  χ₅│ 5  1

   2│ .             .             .
   3│ 1             .             .
   5│ .             1             1
    │                              
    │3a            5a            5b
  2P│3a            5b            5a
  3P│1a            5b            5a
  5P│3a            1a            1a
    │                              
  ──┼──────────────────────────────
  χ₁│ 1             1             1
  χ₂│ . ζ₅³ + ζ₅² + 1    -ζ₅³ - ζ₅²
  χ₃│ .    -ζ₅³ - ζ₅² ζ₅³ + ζ₅² + 1
  χ₄│ 1            -1            -1
  χ₅│-1             .             .
  """

  show(IOContext(io, :separators_col => [0],
                     :separators_row => [0],
                     :portions_col => [2,3]), t_a5)
  @test String(take!(io)) ==
  """
  A5
  
    2| 2  2
    3| 1  .
    5| 1  .
     |     
     |1a 2a
   2P|1a 1a
   3P|1a 2a
   5P|1a 2a
     |     
  ---+-----
  X_1| 1  1
  X_2| 3 -1
  X_3| 3 -1
  X_4| 4  .
  X_5| 5  1
  
    2| .                 .                 .
    3| 1                 .                 .
    5| .                 1                 1
     |                                      
     |3a                5a                5b
   2P|3a                5b                5a
   3P|1a                5b                5a
   5P|3a                1a                1a
     |                                      
  ---+--------------------------------------
  X_1| 1                 1                 1
  X_2| . z_5^3 + z_5^2 + 1    -z_5^3 - z_5^2
  X_3| .    -z_5^3 - z_5^2 z_5^3 + z_5^2 + 1
  X_4| 1                -1                -1
  X_5|-1                 .                 .
  """

  # ... and in LaTeX format
  show(IOContext(io, :separators_col => [0],
                     :separators_row => [0],
                     :portions_col => [2,3]), MIME("text/latex"), t_a5)
  @test String(take!(io)) ==
  """
  \$A5

  \\begin{array}{r|rr}
  2 & 2 & 2 \\\\
  3 & 1 & . \\\\
  5 & 1 & . \\\\
   &  &  \\\\
   & 1a & 2a \\\\
  2P & 1a & 1a \\\\
  3P & 1a & 2a \\\\
  5P & 1a & 2a \\\\
   &  &  \\\\
  \\hline
  \\chi_{1} & 1 & 1 \\\\
  \\chi_{2} & 3 & -1 \\\\
  \\chi_{3} & 3 & -1 \\\\
  \\chi_{4} & 4 & . \\\\
  \\chi_{5} & 5 & 1 \\\\
  \\end{array}

  \\begin{array}{r|rrr}
  2 & . & . & . \\\\
  3 & 1 & . & . \\\\
  5 & . & 1 & 1 \\\\
   &  &  &  \\\\
   & 3a & 5a & 5b \\\\
  2P & 3a & 5b & 5a \\\\
  3P & 1a & 5b & 5a \\\\
  5P & 3a & 1a & 1a \\\\
   &  &  &  \\\\
  \\hline
  \\chi_{1} & 1 & 1 & 1 \\\\
  \\chi_{2} & . & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 & -\\zeta_{5}^{3} - \\zeta_{5}^{2} \\\\
  \\chi_{3} & . & -\\zeta_{5}^{3} - \\zeta_{5}^{2} & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 \\\\
  \\chi_{4} & 1 & -1 & -1 \\\\
  \\chi_{5} & -1 & . & . \\\\
  \\end{array}
  \$"""

  # distribute the table into row portions,
  # in the screen format (perhaps not relevant) ...
  Oscar.with_unicode() do
    show(IOContext(io, :separators_col => [0],
                       :separators_row => [0],
                       :portions_row => [2,3]), t_a5)
  end
  @test String(take!(io)) ==
  """
  A5

   2│ 2  2  .             .             .
   3│ 1  .  1             .             .
   5│ 1  .  .             1             1
    │                                    
    │1a 2a 3a            5a            5b
  2P│1a 1a 3a            5b            5a
  3P│1a 2a 1a            5b            5a
  5P│1a 2a 3a            1a            1a
    │                                    
  ──┼────────────────────────────────────
  χ₁│ 1  1  1             1             1
  χ₂│ 3 -1  . ζ₅³ + ζ₅² + 1    -ζ₅³ - ζ₅²

   2│ 2  2  .             .             .
   3│ 1  .  1             .             .
   5│ 1  .  .             1             1
    │                                    
    │1a 2a 3a            5a            5b
  2P│1a 1a 3a            5b            5a
  3P│1a 2a 1a            5b            5a
  5P│1a 2a 3a            1a            1a
    │                                    
  ──┼────────────────────────────────────
  χ₃│ 3 -1  .    -ζ₅³ - ζ₅² ζ₅³ + ζ₅² + 1
  χ₄│ 4  .  1            -1            -1
  χ₅│ 5  1 -1             .             .
  """

  show(IOContext(io, :separators_col => [0],
                     :separators_row => [0],
                     :portions_row => [2,3]), t_a5)
  @test String(take!(io)) ==
  """
  A5
  
    2| 2  2  .                 .                 .
    3| 1  .  1                 .                 .
    5| 1  .  .                 1                 1
     |                                            
     |1a 2a 3a                5a                5b
   2P|1a 1a 3a                5b                5a
   3P|1a 2a 1a                5b                5a
   5P|1a 2a 3a                1a                1a
     |                                            
  ---+--------------------------------------------
  X_1| 1  1  1                 1                 1
  X_2| 3 -1  . z_5^3 + z_5^2 + 1    -z_5^3 - z_5^2
  
    2| 2  2  .                 .                 .
    3| 1  .  1                 .                 .
    5| 1  .  .                 1                 1
     |                                            
     |1a 2a 3a                5a                5b
   2P|1a 1a 3a                5b                5a
   3P|1a 2a 1a                5b                5a
   5P|1a 2a 3a                1a                1a
     |                                            
  ---+--------------------------------------------
  X_3| 3 -1  .    -z_5^3 - z_5^2 z_5^3 + z_5^2 + 1
  X_4| 4  .  1                -1                -1
  X_5| 5  1 -1                 .                 .
  """

  # ... and in LaTeX format (may be interesting)
  show(IOContext(io, :separators_col => [0],
                     :separators_row => [0],
                     :portions_row => [2,3]), MIME("text/latex"), t_a5)
  @test String(take!(io)) ==
  """\$A5

  \\begin{array}{r|rrrrr}
  2 & 2 & 2 & . & . & . \\\\
  3 & 1 & . & 1 & . & . \\\\
  5 & 1 & . & . & 1 & 1 \\\\
   &  &  &  &  &  \\\\
   & 1a & 2a & 3a & 5a & 5b \\\\
  2P & 1a & 1a & 3a & 5b & 5a \\\\
  3P & 1a & 2a & 1a & 5b & 5a \\\\
  5P & 1a & 2a & 3a & 1a & 1a \\\\
   &  &  &  &  &  \\\\
  \\hline
  \\chi_{1} & 1 & 1 & 1 & 1 & 1 \\\\
  \\chi_{2} & 3 & -1 & . & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 & -\\zeta_{5}^{3} - \\zeta_{5}^{2} \\\\
  \\end{array}

  \\begin{array}{r|rrrrr}
  2 & 2 & 2 & . & . & . \\\\
  3 & 1 & . & 1 & . & . \\\\
  5 & 1 & . & . & 1 & 1 \\\\
   &  &  &  &  &  \\\\
   & 1a & 2a & 3a & 5a & 5b \\\\
  2P & 1a & 1a & 3a & 5b & 5a \\\\
  3P & 1a & 2a & 1a & 5b & 5a \\\\
  5P & 1a & 2a & 3a & 1a & 1a \\\\
   &  &  &  &  &  \\\\
  \\hline
  \\chi_{3} & 3 & -1 & . & -\\zeta_{5}^{3} - \\zeta_{5}^{2} & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 \\\\
  \\chi_{4} & 4 & . & 1 & -1 & -1 \\\\
  \\chi_{5} & 5 & 1 & -1 & . & . \\\\
  \\end{array}
  \$"""

  # show indicators in the screen format ...
  Oscar.with_unicode() do
    show(IOContext(io, :indicator => [2]), t_a4)
  end
  @test String(take!(io)) ==
  """
  Alt( [ 1 .. 4 ] )

      2  2  2       .       .
      3  1  .       1       1
                             
        1a 2a      3a      3b
     2P 1a 1a      3b      3a
     3P 1a 2a      1a      1a
      2                      
  χ₁  +  1  1       1       1
  χ₂  o  1  1 -ζ₃ - 1      ζ₃
  χ₃  o  1  1      ζ₃ -ζ₃ - 1
  χ₄  +  3 -1       .       .
  """

  # ... and in LaTeX format
  show(IOContext(io, :indicator => [2]), MIME("text/latex"), t_a4)
  @test String(take!(io)) ==
  """\$Alt( [ 1 .. 4 ] )

  \\begin{array}{rrrrrr}
   & 2 & 2 & 2 & . & . \\\\
   & 3 & 1 & . & 1 & 1 \\\\
   &  &  &  &  &  \\\\
   &  & 1a & 2a & 3a & 3b \\\\
   & 2P & 1a & 1a & 3b & 3a \\\\
   & 3P & 1a & 2a & 1a & 1a \\\\
   & 2 &  &  &  &  \\\\
  \\chi_{1} & + & 1 & 1 & 1 & 1 \\\\
  \\chi_{2} & o & 1 & 1 & -\\zeta_{3} - 1 & \\zeta_{3} \\\\
  \\chi_{3} & o & 1 & 1 & \\zeta_{3} & -\\zeta_{3} - 1 \\\\
  \\chi_{4} & + & 3 & -1 & . & . \\\\
  \\end{array}
  \$"""

  # ordinary table:
  # show character field degrees in the screen format ...
  Oscar.with_unicode() do
    show(IOContext(io, :character_field => true), t_a4)
  end
  @test String(take!(io)) ==
  """
  Alt( [ 1 .. 4 ] )
  
      2  2  2       .       .
      3  1  .       1       1
                             
        1a 2a      3a      3b
     2P 1a 1a      3b      3a
     3P 1a 2a      1a      1a
      d                      
  χ₁  1  1  1       1       1
  χ₂  2  1  1 -ζ₃ - 1      ζ₃
  χ₃  2  1  1      ζ₃ -ζ₃ - 1
  χ₄  1  3 -1       .       .
  """

  # ... and in LaTeX format
  show(IOContext(io, :character_field => true), MIME("text/latex"), t_a4)
  @test String(take!(io)) ==
  """\$Alt( [ 1 .. 4 ] )

  \\begin{array}{rrrrrr}
   & 2 & 2 & 2 & . & . \\\\
   & 3 & 1 & . & 1 & 1 \\\\
   &  &  &  &  &  \\\\
   &  & 1a & 2a & 3a & 3b \\\\
   & 2P & 1a & 1a & 3b & 3a \\\\
   & 3P & 1a & 2a & 1a & 1a \\\\
   & d &  &  &  &  \\\\
  \\chi_{1} & 1 & 1 & 1 & 1 & 1 \\\\
  \\chi_{2} & 2 & 1 & 1 & -\\zeta_{3} - 1 & \\zeta_{3} \\\\
  \\chi_{3} & 2 & 1 & 1 & \\zeta_{3} & -\\zeta_{3} - 1 \\\\
  \\chi_{4} & 1 & 3 & -1 & . & . \\\\
  \\end{array}
  \$"""

  # Brauer table:
  # show character field degrees in the screen format ...
  Oscar.with_unicode() do
    show(IOContext(io, :character_field => true), mod(t_a4, 2))
  end
  @test String(take!(io)) ==
  """
  Alt( [ 1 .. 4 ] ) mod 2
  
      2  2       .       .
      3  1       1       1
                          
        1a      3a      3b
     2P 1a      3b      3a
     3P 1a      1a      1a
      d                   
  χ₁  1  1       1       1
  χ₂  2  1 -ζ₃ - 1      ζ₃
  χ₃  2  1      ζ₃ -ζ₃ - 1
  """

  # ... and in LaTeX format
  show(IOContext(io, :character_field => true), MIME("text/latex"), mod(t_a4, 2))
  @test String(take!(io)) ==
  """\$Alt( [ 1 .. 4 ] ) mod 2
  
  \\begin{array}{rrrrr}
   & 2 & 2 & . & . \\\\
   & 3 & 1 & 1 & 1 \\\\
   &  &  &  &  \\\\
   &  & 1a & 3a & 3b \\\\
   & 2P & 1a & 3b & 3a \\\\
   & 3P & 1a & 1a & 1a \\\\
   & d &  &  &  \\\\
  \\chi_{1} & 1 & 1 & 1 & 1 \\\\
  \\chi_{2} & 2 & 1 & -\\zeta_{3} - 1 & \\zeta_{3} \\\\
  \\chi_{3} & 2 & 1 & \\zeta_{3} & -\\zeta_{3} - 1 \\\\
  \\end{array}
  \$"""
end

@testset "create character tables" begin
  @testset "library of character tables" begin
    @test character_table("J5") == nothing
    @test_throws ErrorException character_table("A5mod2", 3)
  end

  @test character_table(alternating_group(5), 2) == nothing
end

@testset "characters" begin
  g = symmetric_group(4)
  t = character_table(g)
  @test nrows(t) == 5
  @test ncols(t) == 5
  @test_throws ErrorException t[6]
  tr = trivial_character(t)
  @test tr in t
  @test tr != t[1]
  @test tr == t[5]
  @test tr == trivial_character(g)
  chi = t[2]
  @test chi isa Oscar.GAPGroupClassFunction
  @test chi[4] == t[2,4]
  @test [chi[i] for i in 1:5] == values(chi)
  @test [2*chi[i] for i in 1:5] == values(chi + chi)
  @test [chi[i]^2 for i in 1:5] == values(chi * chi)
  @test [chi[i]^2 for i in 1:5] == values(chi^2)
  @test [-chi[i] for i in 1:5] == values(-chi)
  @test [0*chi[i] for i in 1:5] == values(chi-chi)
  @test [1 for i in 1:5] == values(one(chi))
  @test [0 for i in 1:5] == values(zero(chi))
  @test degree(chi) == chi[1]
  @test length(chi) == ncols(t)
  @test -1 in chi
  chi = filter(x -> x[1] == 2, collect(t))[1]
  @test [chi(representative(c)) for c in conjugacy_classes(g)] == values(chi)
  @test all(is_irreducible, t)

  scp = scalar_product(t[1], t[1])
  @test scp == 1
  @test scp isa fmpq
  for T in [fmpz, fmpq, Int64, QQAbElem]
    scpT = scalar_product(T, t[1],t[1])
    @test scpT == scp
    @test scpT isa T
  end
  scp = scalar_product(t[1], t[2])
  @test scp == 0
  @test scp isa fmpq

  # conjugate characters
  h = pcore(g, 2)[1]
  th = character_table(h)
  chi = th[2]
  ggens = gens(g)
  x = ggens[1]*ggens[2]
  psi = chi^x
  @test all(y -> chi(x*y*x^-1) == psi(y), collect(h))

  # a corresponding Brauer table (possible because the group is solvable)
  tmod2 = mod(t, 2);
  @test tmod2 isa Oscar.GAPGroupCharacterTable
  chi = tmod2[1]
  @test degree(chi) == chi[1]
  @test length(chi) == ncols(tmod2)
  @test chi[2] == tmod2[1,2]
  tr = trivial_character(tmod2)
  @test ! (-1 in tr)

  mat = decomposition_matrix(tmod2)
  @test nrows(mat) == nrows(t)
  @test ncols(mat) == nrows(tmod2)

  g = symmetric_group(4)
  t = character_table(g)

  g = symmetric_group(5)
  t = character_table(g)
  @test mod(t, 2) == nothing

  g = general_linear_group(2, 3)
  h = derived_subgroup(g)[1]
  p = 2
  t = mod(character_table(h), p)
  chi = t[2]
  x = gen(g, 1)
  @test ! (x in h)
  psi = chi^x
  @test all(y -> mod(order(y), p) == 0 || chi(x*y*x^-1) == psi(y), collect(h))
end

@testset "Galois conjugacy of characters" begin
  g = alternating_group(4)
  t = character_table(g)
  @test all(x -> conj(x) == QQAbAutomorphism(5)(x), t)
  @test all(x -> x == QQAbAutomorphism(4)(x), t)
end

@testset "induction of characters" begin
  g = symmetric_group(4)
  t = character_table(g)
  indcyc = induced_cyclic(t)
  @test sort!([degree(chi) for chi in indcyc]) == [6, 8, 12, 12, 24]
  @test all(x -> scalar_product(trivial_character(t), x) == 1, indcyc)

  h = derived_subgroup(g)[1]
  ind = [chi^t for chi in character_table(h)]
  @test sort!([degree(chi) for chi in ind]) == [2, 2, 2, 6]
  @test [scalar_product(trivial_character(t), x) for x in ind] == [1, 0, 0, 0]
end

@testset "natural characters" begin
  G = symmetric_group(4)
  chi = Oscar.natural_character(G)
  @test degree(chi) == 4
  psi = chi
  @test scalar_product(chi, psi) == 2
  @test scalar_product(chi, chi) == 2
  @test scalar_product(chi, trivial_character(G)) == 1

  K, a = CyclotomicField(5, "a")
  L, b = CyclotomicField(3, "b")

  inputs = [
    #[ matrix(ZZ, [0 1 0; -1 0 0; 0 0 -1]) ],
    [ matrix(QQ, [0 1 0; -1 0 0; 0 0 -1]) ],
    [ matrix(K, [a 0; 0 a]) ],
    [ matrix(L, 2, 2, [b, 0, -b - 1, 1]), matrix(L, 2, 2, [1, b + 1, 0, b]) ],
  ]

  @testset "... over ring $(base_ring(mats[1]))" for mats in inputs
    G = matrix_group(mats)
    chi = Oscar.natural_character(G)
    @test degree(chi) == degree(G)
  end
end

@testset "character fields" begin
  for id in [ "C5", "A5" ]   # cyclotomic and non-cyclotomic number fields
    for chi in character_table(id)
      F, phi = character_field(chi)
      for i in 1:length(chi)
        x = chi[i]
        xF = preimage(phi, x)
        @test parent(xF) == F
        @test phi(xF) == x
      end
    end
  end
end

@testset "Schur index" begin
  t = character_table("2.A5")
  @test map(schur_index, collect(t)) == [1,1,1,1,1,2,2,2,2]

  # Test a character table with group.
  g = small_group(192, 1022)
  h = derived_subgroup(g)[1];
  s = character_table(h);
  t = character_table(g);
  trivial_character(s)^t;  # side-effect: stores a class fusion
  @test length(names_of_fusion_sources(t)) > 0
  if hasproperty(GAP.Globals, :SchurIndexByCharacter)
    # We can compute the values.
    @test sort!(map(schur_index, collect(t))) == append!(repeat([1],15), repeat([2],4))
  else
    # We can fail.
    for chi in t
      try
        schur_index(chi)
      catch(e)
        msg = sprint(showerror, e)
        @test msg == "cannot determine the Schur index with the currently used criteria"
      end
    end
  end

  # For a character table without group, we can fail.
  t = character_table("S6")
  for chi in t
    try
      schur_index(chi)
    catch(e)
      msg = sprint(showerror, e)
      @test msg == "cannot determine the Schur index with the currently used criteria"
    end
  end
end

@testset "specialized generic tables" begin
  inputs = [(:Cyclic, 3),
            (:Dihedral, 8),
            (:Symmetric, 4),
            (:Alternating, 4),
            (:WeylB, 3),
            (:WeylD, 3),
            (:DoubleCoverSymmetric, 5),
            (:DoubleCoverAlternating, 5),
            (:GL2, 3),
            (:SL2odd, 7),
            (:SL2even, 4),
            (:PSL2odd, 7),
            (:PSL2even, 5),
            (:Suzuki, 8),
            (:GU3, 2),
            (:SU3, 2),
            (Symbol("P:Q"), [5, 4]),
            (:ExtraspecialPlusOdd, 27),
           ];
  for (series, para) in inputs
    t = character_table(series, para)
    @test length(character_parameters(t)) == length(t)
    @test length(class_parameters(t)) == length(t)
  end
end

@testset "symmetrizations" begin
    t = character_table("S5");
    irr = [chi for chi in t];
    @test [chi[1] for chi in symmetrizations(irr, 2)] == [0, 1, 0, 1, 15, 21, 6, 10, 6, 10, 10, 15, 10, 15]
    @test [chi[1] for chi in symmetric_parts(irr, 2)] == [1, 1, 21, 10, 10, 15, 15]
    @test [chi[1] for chi in anti_symmetric_parts(irr, 2)] == [0, 0, 15, 6, 6, 10, 10]
    @test [chi[1] for chi in orthogonal_components(irr, 2)] == [15, 20, 6, 9, 6, 9, 10, 14, 10, 14]
    @test [chi[1] for chi in symplectic_components(irr, 2)] == [14, 21, 5, 10, 5, 10, 9, 15, 9, 15]
    @test [exterior_power(chi, 3)[1] for chi in t] == [0, 0, 20, 4, 4, 10, 10]
    @test [symmetric_power(chi, 3)[1] for chi in t] == [1, 1, 56, 20, 20, 35, 35]

    empty = Oscar.GAPGroupClassFunction[]
    @test symmetrizations(empty, 2) == empty
    @test symmetric_parts(empty, 2) == empty
    @test anti_symmetric_parts(empty, 2) == empty
    @test orthogonal_components(empty, 2) == empty
    @test symplectic_components(empty, 2) == empty
end
