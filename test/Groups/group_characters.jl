@testset "show and print character tables" begin
  io = IOBuffer();
  t_a4 = character_table(alternating_group(4))
  t_a5 = character_table("A5")

  # `print` shows an abbrev. form
  print(io, t_a4)
  @test String(take!(io)) == "character_table(Alt( [ 1 .. 4 ] ))"

  # default `show`
  show(io, t_a4)
  @test String(take!(io)) == "Alt( [ 1 .. 4 ] )\n\n 2  2  2       .       .\n 3  1  .       1       1\n                        \n   1a 2a      3a      3b\n2P 1a 1a      3b      3a\n3P 1a 2a      1a      1a\n                        \nχ₁  1  1       1       1\nχ₂  1  1 -ζ₃ - 1      ζ₃\nχ₃  1  1      ζ₃ -ζ₃ - 1\nχ₄  3 -1       .       .\n"

  # LaTeX format
  show(io, MIME("text/html"), t_a4)
  @test String(take!(io)) == "\$Alt( [ 1 .. 4 ] )\n\n\\begin{array}{rrrrr}\n2 & 2 & 2 & . & . \\\\\n3 & 1 & . & 1 & 1 \\\\\n &  &  &  &  \\\\\n & 1a & 2a & 3a & 3b \\\\\n2P & 1a & 1a & 3b & 3a \\\\\n3P & 1a & 2a & 1a & 1a \\\\\n &  &  &  &  \\\\\n\\chi_{1} & 1 & 1 & 1 & 1 \\\\\n\\chi_{2} & 1 & 1 & -\\zeta_{3} - 1 & \\zeta_{3} \\\\\n\\chi_{3} & 1 & 1 & \\zeta_{3} & -\\zeta_{3} - 1 \\\\\n\\chi_{4} & 3 & -1 & . & . \\\\\n\\end{array}\n\$"

  # show a legend of irrationalities instead of self-explanatory values,
  # in the screen format ...
  show(IOContext(io, :with_legend => true), t_a4)
  @test String(take!(io)) == "Alt( [ 1 .. 4 ] )\n\n 2  2  2  .  .\n 3  1  .  1  1\n              \n   1a 2a 3a 3b\n2P 1a 1a 3b 3a\n3P 1a 2a 1a 1a\n              \nχ₁  1  1  1  1\nχ₂  1  1  A  A̅\nχ₃  1  1  A̅  A\nχ₄  3 -1  .  .\n\nA = -ζ₃ - 1\nA̅ = ζ₃\n"

  # ... and in LaTeX format
  show(IOContext(io, :with_legend => true), MIME("text/html"), t_a4)
  @test String(take!(io)) == "\$Alt( [ 1 .. 4 ] )\n\n\\begin{array}{rrrrr}\n2 & 2 & 2 & . & . \\\\\n3 & 1 & . & 1 & 1 \\\\\n &  &  &  &  \\\\\n & 1a & 2a & 3a & 3b \\\\\n2P & 1a & 1a & 3b & 3a \\\\\n3P & 1a & 2a & 1a & 1a \\\\\n &  &  &  &  \\\\\n\\chi_{1} & 1 & 1 & 1 & 1 \\\\\n\\chi_{2} & 1 & 1 & A & \\overline{A} \\\\\n\\chi_{3} & 1 & 1 & \\overline{A} & A \\\\\n\\chi_{4} & 3 & -1 & . & . \\\\\n\\end{array}\n\nA = -\\zeta_{3} - 1\n\\overline{A} = \\zeta_{3}\n\$"

  # show the screen format for a table with real and non-real irrationalities
  show(IOContext(io, :with_legend => true), character_table("L2(11)"))
  @test String(take!(io)) == "L2(11)\n\n  2  2  2  1  .  .  1   .   .\n  3  1  1  1  .  .  1   .   .\n  5  1  .  .  1  1  .   .   .\n 11  1  .  .  .  .  .   1   1\n                             \n    1a 2a 3a 5a 5b 6a 11a 11b\n 2P 1a 1a 3a 5b 5a 3a 11b 11a\n 3P 1a 2a 1a 5b 5a 2a 11a 11b\n 5P 1a 2a 3a 1a 1a 6a 11a 11b\n11P 1a 2a 3a 5a 5b 6a  1a  1a\n                             \n χ₁  1  1  1  1  1  1   1   1\n χ₂  5  1 -1  .  .  1   C   C̅\n χ₃  5  1 -1  .  .  1   C̅   C\n χ₄ 10 -2  1  .  .  1  -1  -1\n χ₅ 10  2  1  .  . -1  -1  -1\n χ₆ 11 -1 -1  1  1 -1   .   .\n χ₇ 12  .  .  A  B  .   1   1\n χ₈ 12  .  .  B  A  .   1   1\n\nA = -ζ₅³ - ζ₅² - 1\nB = ζ₅³ + ζ₅²\nC = ζ₁₁⁹ + ζ₁₁⁵ + ζ₁₁⁴ + ζ₁₁³ + ζ₁₁\nC̅ = -ζ₁₁⁹ - ζ₁₁⁵ - ζ₁₁⁴ - ζ₁₁³ - ζ₁₁ - 1\n"

  # show some separating lines, in the screen format ...
  show(IOContext(io, :separators_col => [0,5],
                     :separators_row => [0,5]), t_a5)
  @test String(take!(io)) == "A5\n\n 2│ 2  2  .             .             .│\n 3│ 1  .  1             .             .│\n 5│ 1  .  .             1             1│\n  │                                    │\n  │1a 2a 3a            5a            5b│\n2P│1a 1a 3a            5b            5a│\n3P│1a 2a 1a            5b            5a│\n5P│1a 2a 3a            1a            1a│\n  │                                    │\n──┼────────────────────────────────────┼\nχ₁│ 1  1  1             1             1│\nχ₂│ 3 -1  . ζ₅³ + ζ₅² + 1    -ζ₅³ - ζ₅²│\nχ₃│ 3 -1  .    -ζ₅³ - ζ₅² ζ₅³ + ζ₅² + 1│\nχ₄│ 4  .  1            -1            -1│\nχ₅│ 5  1 -1             .             .│\n──┼────────────────────────────────────┼\n"

  # ... and in LaTeX format
  show(IOContext(io, :separators_col => [0,5],
                     :separators_row => [0,5]),
                     MIME("text/html"), t_a5)
  @test String(take!(io)) == "\$A5\n\n\\begin{array}{r|rrrrr|}\n2 & 2 & 2 & . & . & . \\\\\n3 & 1 & . & 1 & . & . \\\\\n5 & 1 & . & . & 1 & 1 \\\\\n &  &  &  &  &  \\\\\n & 1a & 2a & 3a & 5a & 5b \\\\\n2P & 1a & 1a & 3a & 5b & 5a \\\\\n3P & 1a & 2a & 1a & 5b & 5a \\\\\n5P & 1a & 2a & 3a & 1a & 1a \\\\\n &  &  &  &  &  \\\\\n\\hline\n\\chi_{1} & 1 & 1 & 1 & 1 & 1 \\\\\n\\chi_{2} & 3 & -1 & . & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 & -\\zeta_{5}^{3} - \\zeta_{5}^{2} \\\\\n\\chi_{3} & 3 & -1 & . & -\\zeta_{5}^{3} - \\zeta_{5}^{2} & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 \\\\\n\\chi_{4} & 4 & . & 1 & -1 & -1 \\\\\n\\chi_{5} & 5 & 1 & -1 & . & . \\\\\n\\hline\n\\end{array}\n\$"

  # distribute the table into column portions, in the screen format ...
  show(IOContext(io, :separators_col => [0],
                     :separators_row => [0],
                     :portions_col => [2,3]), t_a5)
  @test String(take!(io)) == "A5\n\n 2│ 2  2\n 3│ 1  .\n 5│ 1  .\n  │     \n  │1a 2a\n2P│1a 1a\n3P│1a 2a\n5P│1a 2a\n  │     \n──┼─────\nχ₁│ 1  1\nχ₂│ 3 -1\nχ₃│ 3 -1\nχ₄│ 4  .\nχ₅│ 5  1\n\n 2│ .             .             .\n 3│ 1             .             .\n 5│ .             1             1\n  │                              \n  │3a            5a            5b\n2P│3a            5b            5a\n3P│1a            5b            5a\n5P│3a            1a            1a\n  │                              \n──┼──────────────────────────────\nχ₁│ 1             1             1\nχ₂│ . ζ₅³ + ζ₅² + 1    -ζ₅³ - ζ₅²\nχ₃│ .    -ζ₅³ - ζ₅² ζ₅³ + ζ₅² + 1\nχ₄│ 1            -1            -1\nχ₅│-1             .             .\n"

  # ... and in LaTeX format
  show(IOContext(io, :separators_col => [0],
                     :separators_row => [0],
                     :portions_col => [2,3]), MIME("text/html"), t_a5)
  @test String(take!(io)) == "\$A5\n\n\\begin{array}{r|rr}\n2 & 2 & 2 \\\\\n3 & 1 & . \\\\\n5 & 1 & . \\\\\n &  &  \\\\\n & 1a & 2a \\\\\n2P & 1a & 1a \\\\\n3P & 1a & 2a \\\\\n5P & 1a & 2a \\\\\n &  &  \\\\\n\\hline\n\\chi_{1} & 1 & 1 \\\\\n\\chi_{2} & 3 & -1 \\\\\n\\chi_{3} & 3 & -1 \\\\\n\\chi_{4} & 4 & . \\\\\n\\chi_{5} & 5 & 1 \\\\\n\\end{array}\n\n\\begin{array}{r|rrr}\n2 & . & . & . \\\\\n3 & 1 & . & . \\\\\n5 & . & 1 & 1 \\\\\n &  &  &  \\\\\n & 3a & 5a & 5b \\\\\n2P & 3a & 5b & 5a \\\\\n3P & 1a & 5b & 5a \\\\\n5P & 3a & 1a & 1a \\\\\n &  &  &  \\\\\n\\hline\n\\chi_{1} & 1 & 1 & 1 \\\\\n\\chi_{2} & . & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 & -\\zeta_{5}^{3} - \\zeta_{5}^{2} \\\\\n\\chi_{3} & . & -\\zeta_{5}^{3} - \\zeta_{5}^{2} & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 \\\\\n\\chi_{4} & 1 & -1 & -1 \\\\\n\\chi_{5} & -1 & . & . \\\\\n\\end{array}\n\$"

  # distribute the table into row portions,
  # in the screen format (perhaps not relevant) ...
  show(IOContext(io, :separators_col => [0],
                     :separators_row => [0],
                     :portions_row => [2,3]), t_a5)
  @test String(take!(io)) == "A5\n\n 2│ 2  2  .             .             .\n 3│ 1  .  1             .             .\n 5│ 1  .  .             1             1\n  │                                    \n  │1a 2a 3a            5a            5b\n2P│1a 1a 3a            5b            5a\n3P│1a 2a 1a            5b            5a\n5P│1a 2a 3a            1a            1a\n  │                                    \n──┼────────────────────────────────────\nχ₁│ 1  1  1             1             1\nχ₂│ 3 -1  . ζ₅³ + ζ₅² + 1    -ζ₅³ - ζ₅²\n\n 2│ 2  2  .             .             .\n 3│ 1  .  1             .             .\n 5│ 1  .  .             1             1\n  │                                    \n  │1a 2a 3a            5a            5b\n2P│1a 1a 3a            5b            5a\n3P│1a 2a 1a            5b            5a\n5P│1a 2a 3a            1a            1a\n  │                                    \n──┼────────────────────────────────────\nχ₃│ 3 -1  .    -ζ₅³ - ζ₅² ζ₅³ + ζ₅² + 1\nχ₄│ 4  .  1            -1            -1\nχ₅│ 5  1 -1             .             .\n"

  # ... and in LaTeX format (may be interesting)
  show(IOContext(io, :separators_col => [0],
                     :separators_row => [0],
                     :portions_row => [2,3]), MIME("text/html"), t_a5)
  @test String(take!(io)) == "\$A5\n\n\\begin{array}{r|rrrrr}\n2 & 2 & 2 & . & . & . \\\\\n3 & 1 & . & 1 & . & . \\\\\n5 & 1 & . & . & 1 & 1 \\\\\n &  &  &  &  &  \\\\\n & 1a & 2a & 3a & 5a & 5b \\\\\n2P & 1a & 1a & 3a & 5b & 5a \\\\\n3P & 1a & 2a & 1a & 5b & 5a \\\\\n5P & 1a & 2a & 3a & 1a & 1a \\\\\n &  &  &  &  &  \\\\\n\\hline\n\\chi_{1} & 1 & 1 & 1 & 1 & 1 \\\\\n\\chi_{2} & 3 & -1 & . & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 & -\\zeta_{5}^{3} - \\zeta_{5}^{2} \\\\\n\\end{array}\n\n\\begin{array}{r|rrrrr}\n2 & 2 & 2 & . & . & . \\\\\n3 & 1 & . & 1 & . & . \\\\\n5 & 1 & . & . & 1 & 1 \\\\\n &  &  &  &  &  \\\\\n & 1a & 2a & 3a & 5a & 5b \\\\\n2P & 1a & 1a & 3a & 5b & 5a \\\\\n3P & 1a & 2a & 1a & 5b & 5a \\\\\n5P & 1a & 2a & 3a & 1a & 1a \\\\\n &  &  &  &  &  \\\\\n\\hline\n\\chi_{3} & 3 & -1 & . & -\\zeta_{5}^{3} - \\zeta_{5}^{2} & \\zeta_{5}^{3} + \\zeta_{5}^{2} + 1 \\\\\n\\chi_{4} & 4 & . & 1 & -1 & -1 \\\\\n\\chi_{5} & 5 & 1 & -1 & . & . \\\\\n\\end{array}\n\$"
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
