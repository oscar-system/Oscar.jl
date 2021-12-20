@testset "show and print character tables" begin
  io = IOBuffer();
  t_a4 = character_table(alternating_group(4))
  t_a5 = character_table("A5")

  # `print` shows an abbrev. form
  print(io, t_a4)
  @test String(take!(io)) == "character_table(Alt( [ 1 .. 4 ] ))"

  # default `show`
  show(io, t_a4)
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
  show(IOContext(io, :with_legend => true), t_a4)
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

  # show some separating lines, in the screen format ...
  show(IOContext(io, :separators_col => [0,5],
                     :separators_row => [0,5]), t_a5)
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
  show(IOContext(io, :separators_col => [0],
                     :separators_row => [0],
                     :portions_col => [2,3]), t_a5)
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
  show(IOContext(io, :separators_col => [0],
                     :separators_row => [0],
                     :portions_row => [2,3]), t_a5)
  @test String(take!(io)) ==
  """A5

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

  scp = scalar_product(t[1], t[1])
  @test scp == 1
  @test scp isa fmpz
  scp = scalar_product(t[1], t[2])
  @test scp == 0
  @test scp isa fmpz

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
