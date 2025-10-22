# Make sure that the tests start without cached character tables.
# (Running the tests will store information that changes some test outputs,
# thus running the tests twice needs these calls.)
GAP.Globals.UnloadCharacterTableData()
empty!(Oscar.character_tables_by_id)

Oscar.@_AuxDocTest "show and print character tables", (fix = false),
raw"""
```jldoctest group_characters.test
julia> using Oscar

julia> t_a4 = Oscar._sort_for_stable_tests(character_table(alternating_group(4)));

julia> t_a5 = Oscar._sort_for_stable_tests(character_table("A5"));

julia> t_a4_2 = mod(t_a4, 2);

julia> t_a5_2 = mod(t_a5, 2);
```

`print` shows an abbrev. form

```jldoctest group_characters.test
julia> print(t_a4)
character table of Alt(4)

julia> print(t_a5)
character table of A5

julia> print(t_a4_2)
2-modular Brauer table of Alt(4)

julia> print(t_a5_2)
2-modular Brauer table of A5
```

`show` uses the abbrev. form for nested objects

```jldoctest group_characters.test
julia> show([t_a4])
Oscar.GAPGroupCharacterTable[character table of Alt(4)]

julia> show([t_a5])
Oscar.GAPGroupCharacterTable[character table of A5]

julia> show([t_a4_2])
Oscar.GAPGroupCharacterTable[2-modular Brauer table of Alt(4)]

julia> show([t_a5_2])
Oscar.GAPGroupCharacterTable[2-modular Brauer table of A5]
```

terse printing
```jldoctest group_characters.test
julia> print(AbstractAlgebra.terse(stdout), t_a4)
character table of a group

julia> print(AbstractAlgebra.terse(stdout), t_a5)
character table of a group

julia> print(AbstractAlgebra.terse(stdout), t_a4_2)
2-modular Brauer table of a group

julia> print(AbstractAlgebra.terse(stdout), t_a5_2)
2-modular Brauer table of a group
```

default `show` with unicode

```jldoctest group_characters.test
julia> Oscar.with_unicode() do
         show(stdout, MIME("text/plain"), t_a4)
       end
Character table of alternating group of degree 4

 2  2  2       .       .
 3  1  .       1       1

   1a 2a      3a      3b
2P 1a 1a      3b      3a
3P 1a 2a      1a      1a

χ₁  1  1       1       1
χ₂  1  1 -ζ₃ - 1      ζ₃
χ₃  1  1      ζ₃ -ζ₃ - 1
χ₄  3 -1       .       .
```

default `show` without unicode

```jldoctest group_characters.test
julia> show(stdout, MIME("text/plain"), t_a4)
Character table of alternating group of degree 4

  2  2  2        .        .
  3  1  .        1        1

    1a 2a       3a       3b
 2P 1a 1a       3b       3a
 3P 1a 2a       1a       1a

X_1  1  1        1        1
X_2  1  1 -z_3 - 1      z_3
X_3  1  1      z_3 -z_3 - 1
X_4  3 -1        .        .
```

LaTeX format

```jldoctest group_characters.test
julia> show(stdout, MIME("text/latex"), t_a4)
Character table of alternating group of degree 4

$\begin{array}{rrrrr}
2 & 2 & 2 & . & . \\
3 & 1 & . & 1 & 1 \\
 &  &  &  &  \\
 & 1a & 2a & 3a & 3b \\
2P & 1a & 1a & 3b & 3a \\
3P & 1a & 2a & 1a & 1a \\
 &  &  &  &  \\
\chi_{1} & 1 & 1 & 1 & 1 \\
\chi_{2} & 1 & 1 & -\zeta_{3} - 1 & \zeta_{3} \\
\chi_{3} & 1 & 1 & \zeta_{3} & -\zeta_{3} - 1 \\
\chi_{4} & 3 & -1 & . & . \\
\end{array}
$

```

show a legend of irrationalities instead of self-explanatory values,
in the screen format ...

```jldoctest group_characters.test
julia> Oscar.with_unicode() do
         show(IOContext(stdout, :with_legend => true), MIME("text/plain"), t_a4)
       end
Character table of alternating group of degree 4

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
```

```jldoctest group_characters.test
julia> show(IOContext(stdout, :with_legend => true), MIME("text/plain"), t_a4)
Character table of alternating group of degree 4

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
```

... and in LaTeX format
```jldoctest group_characters.test
julia> show(IOContext(stdout, :with_legend => true), MIME("text/latex"), t_a4)
Character table of alternating group of degree 4

$\begin{array}{rrrrr}
2 & 2 & 2 & . & . \\
3 & 1 & . & 1 & 1 \\
 &  &  &  &  \\
 & 1a & 2a & 3a & 3b \\
2P & 1a & 1a & 3b & 3a \\
3P & 1a & 2a & 1a & 1a \\
 &  &  &  &  \\
\chi_{1} & 1 & 1 & 1 & 1 \\
\chi_{2} & 1 & 1 & A & \overline{A} \\
\chi_{3} & 1 & 1 & \overline{A} & A \\
\chi_{4} & 3 & -1 & . & . \\
\end{array}

\begin{array}{l}
A = -\zeta_{3} - 1 \\
\overline{A} = \zeta_{3} \\
\end{array}
$
```

Show the screen format for a table with real and non-real irrationalities.
```jldoctest group_characters.test
julia> Oscar.with_unicode() do
         show(IOContext(stdout, :with_legend => true),
              MIME("text/plain"), character_table("L2(11)"))
       end
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
```

```jldoctest group_characters.test
julia> show(IOContext(stdout, :with_legend => true),
            MIME("text/plain"), character_table("L2(11)"))
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
```

Show some separating lines, in the screen format ...
```jldoctest group_characters.test
julia> Oscar.with_unicode() do
         show(IOContext(stdout, :separators_col => [0,5],
                                :separators_row => [0,5]),
              MIME("text/plain"), t_a5)
       end
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
```

```jldoctest group_characters.test
julia> show(IOContext(stdout, :separators_col => [0,5],
                          :separators_row => [0,5]),
            MIME("text/plain"), t_a5)
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
```

... and in LaTeX format

```jldoctest group_characters.test
julia> show(IOContext(stdout, :separators_col => [0,5],
                          :separators_row => [0,5]),
            MIME("text/latex"), t_a5)
A5

$\begin{array}{r|rrrrr|}
2 & 2 & 2 & . & . & . \\
3 & 1 & . & 1 & . & . \\
5 & 1 & . & . & 1 & 1 \\
 &  &  &  &  &  \\
 & 1a & 2a & 3a & 5a & 5b \\
2P & 1a & 1a & 3a & 5b & 5a \\
3P & 1a & 2a & 1a & 5b & 5a \\
5P & 1a & 2a & 3a & 1a & 1a \\
 &  &  &  &  &  \\
\hline
\chi_{1} & 1 & 1 & 1 & 1 & 1 \\
\chi_{2} & 3 & -1 & . & \zeta_{5}^{3} + \zeta_{5}^{2} + 1 & -\zeta_{5}^{3} - \zeta_{5}^{2} \\
\chi_{3} & 3 & -1 & . & -\zeta_{5}^{3} - \zeta_{5}^{2} & \zeta_{5}^{3} + \zeta_{5}^{2} + 1 \\
\chi_{4} & 4 & . & 1 & -1 & -1 \\
\chi_{5} & 5 & 1 & -1 & . & . \\
\hline
\end{array}
$
```

distribute the table into column portions, in the screen format ...

```jldoctest group_characters.test
julia> Oscar.with_unicode() do
         show(IOContext(stdout, :separators_col => [0],
                            :separators_row => [0],
                            :portions_col => [2,3]),
              MIME("text/plain"), t_a5)
       end
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
```


```jldoctest group_characters.test
julia> show(IOContext(stdout, :separators_col => [0],
                     :separators_row => [0],
                     :portions_col => [2,3]),
       MIME("text/plain"), t_a5)
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
```

... and in LaTeX format

```jldoctest group_characters.test
julia> show(IOContext(stdout, :separators_col => [0],
                     :separators_row => [0],
                     :portions_col => [2,3]),
       MIME("text/latex"), t_a5)
A5

$\begin{array}{r|rr}
2 & 2 & 2 \\
3 & 1 & . \\
5 & 1 & . \\
 &  &  \\
 & 1a & 2a \\
2P & 1a & 1a \\
3P & 1a & 2a \\
5P & 1a & 2a \\
 &  &  \\
\hline
\chi_{1} & 1 & 1 \\
\chi_{2} & 3 & -1 \\
\chi_{3} & 3 & -1 \\
\chi_{4} & 4 & . \\
\chi_{5} & 5 & 1 \\
\end{array}

\begin{array}{r|rrr}
2 & . & . & . \\
3 & 1 & . & . \\
5 & . & 1 & 1 \\
 &  &  &  \\
 & 3a & 5a & 5b \\
2P & 3a & 5b & 5a \\
3P & 1a & 5b & 5a \\
5P & 3a & 1a & 1a \\
 &  &  &  \\
\hline
\chi_{1} & 1 & 1 & 1 \\
\chi_{2} & . & \zeta_{5}^{3} + \zeta_{5}^{2} + 1 & -\zeta_{5}^{3} - \zeta_{5}^{2} \\
\chi_{3} & . & -\zeta_{5}^{3} - \zeta_{5}^{2} & \zeta_{5}^{3} + \zeta_{5}^{2} + 1 \\
\chi_{4} & 1 & -1 & -1 \\
\chi_{5} & -1 & . & . \\
\end{array}
$
```

distribute the table into row portions,
in the screen format (perhaps not relevant) ...

```jldoctest group_characters.test
julia> Oscar.with_unicode() do
         show(IOContext(stdout, :separators_col => [0],
                            :separators_row => [0],
                            :portions_row => [2,3]),
              MIME("text/plain"), t_a5)
       end
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
```


```jldoctest group_characters.test
julia> show(IOContext(stdout, :separators_col => [0],
                          :separators_row => [0],
                          :portions_row => [2,3]),
            MIME("text/plain"), t_a5)
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
```

... and in LaTeX format (may be interesting)

```jldoctest group_characters.test
julia> show(IOContext(stdout, :separators_col => [0],
                          :separators_row => [0],
                          :portions_row => [2,3]),
            MIME("text/latex"), t_a5)
A5

$\begin{array}{r|rrrrr}
2 & 2 & 2 & . & . & . \\
3 & 1 & . & 1 & . & . \\
5 & 1 & . & . & 1 & 1 \\
 &  &  &  &  &  \\
 & 1a & 2a & 3a & 5a & 5b \\
2P & 1a & 1a & 3a & 5b & 5a \\
3P & 1a & 2a & 1a & 5b & 5a \\
5P & 1a & 2a & 3a & 1a & 1a \\
 &  &  &  &  &  \\
\hline
\chi_{1} & 1 & 1 & 1 & 1 & 1 \\
\chi_{2} & 3 & -1 & . & \zeta_{5}^{3} + \zeta_{5}^{2} + 1 & -\zeta_{5}^{3} - \zeta_{5}^{2} \\
\end{array}

\begin{array}{r|rrrrr}
2 & 2 & 2 & . & . & . \\
3 & 1 & . & 1 & . & . \\
5 & 1 & . & . & 1 & 1 \\
 &  &  &  &  &  \\
 & 1a & 2a & 3a & 5a & 5b \\
2P & 1a & 1a & 3a & 5b & 5a \\
3P & 1a & 2a & 1a & 5b & 5a \\
5P & 1a & 2a & 3a & 1a & 1a \\
 &  &  &  &  &  \\
\hline
\chi_{3} & 3 & -1 & . & -\zeta_{5}^{3} - \zeta_{5}^{2} & \zeta_{5}^{3} + \zeta_{5}^{2} + 1 \\
\chi_{4} & 4 & . & 1 & -1 & -1 \\
\chi_{5} & 5 & 1 & -1 & . & . \\
\end{array}
$
```

show indicators in the screen format ...

```jldoctest group_characters.test
julia> Oscar.with_unicode() do
         show(IOContext(stdout, :indicator => [2]), MIME("text/plain"), t_a4)
       end
Character table of alternating group of degree 4

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
```

... and in LaTeX format

```jldoctest group_characters.test
julia> show(IOContext(stdout, :indicator => [2]), MIME("text/latex"), t_a4)
Character table of alternating group of degree 4

$\begin{array}{rrrrrr}
 & 2 & 2 & 2 & . & . \\
 & 3 & 1 & . & 1 & 1 \\
 &  &  &  &  &  \\
 &  & 1a & 2a & 3a & 3b \\
 & 2P & 1a & 1a & 3b & 3a \\
 & 3P & 1a & 2a & 1a & 1a \\
 & 2 &  &  &  &  \\
\chi_{1} & + & 1 & 1 & 1 & 1 \\
\chi_{2} & o & 1 & 1 & -\zeta_{3} - 1 & \zeta_{3} \\
\chi_{3} & o & 1 & 1 & \zeta_{3} & -\zeta_{3} - 1 \\
\chi_{4} & + & 3 & -1 & . & . \\
\end{array}
$
```

ordinary table:
show character field degrees in the screen format ...

```jldoctest group_characters.test
julia> Oscar.with_unicode() do
         show(IOContext(stdout, :character_field => true), MIME("text/plain"), t_a4)
       end
Character table of alternating group of degree 4

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
```

... and in LaTeX format

```jldoctest group_characters.test
julia> show(IOContext(stdout, :character_field => true), MIME("text/latex"), t_a4)
Character table of alternating group of degree 4

$\begin{array}{rrrrrr}
 & 2 & 2 & 2 & . & . \\
 & 3 & 1 & . & 1 & 1 \\
 &  &  &  &  &  \\
 &  & 1a & 2a & 3a & 3b \\
 & 2P & 1a & 1a & 3b & 3a \\
 & 3P & 1a & 2a & 1a & 1a \\
 & d &  &  &  &  \\
\chi_{1} & 1 & 1 & 1 & 1 & 1 \\
\chi_{2} & 2 & 1 & 1 & -\zeta_{3} - 1 & \zeta_{3} \\
\chi_{3} & 2 & 1 & 1 & \zeta_{3} & -\zeta_{3} - 1 \\
\chi_{4} & 1 & 3 & -1 & . & . \\
\end{array}
$
```

Brauer table:
show character field degrees in the screen format ...

```jldoctest group_characters.test
julia> Oscar.with_unicode() do
         show(IOContext(stdout, :character_field => true), MIME("text/plain"), mod(t_a4, 2))
       end
2-modular Brauer table of alternating group of degree 4

    2  2       .       .
    3  1       1       1

      1a      3a      3b
   2P 1a      3b      3a
   3P 1a      1a      1a
    d
χ₁  1  1       1       1
χ₂  2  1 -ζ₃ - 1      ζ₃
χ₃  2  1      ζ₃ -ζ₃ - 1
```

... and in LaTeX format

```jldoctest group_characters.test
julia> show(IOContext(stdout, :character_field => true), MIME("text/latex"), mod(t_a4, 2))
2-modular Brauer table of alternating group of degree 4

$\begin{array}{rrrrr}
 & 2 & 2 & . & . \\
 & 3 & 1 & 1 & 1 \\
 &  &  &  &  \\
 &  & 1a & 3a & 3b \\
 & 2P & 1a & 3b & 3a \\
 & 3P & 1a & 1a & 1a \\
 & d &  &  &  \\
\chi_{1} & 1 & 1 & 1 & 1 \\
\chi_{2} & 2 & 1 & -\zeta_{3} - 1 & \zeta_{3} \\
\chi_{3} & 2 & 1 & \zeta_{3} & -\zeta_{3} - 1 \\
\end{array}
$
```

Test the case where a group has a custom name where the first character
should not be turned into lowercase
```jldoctest group_characters.test
julia> character_table(SL(2,2))
Character table of SL(2,2)

  2  1  1  .
  3  1  .  1

    1a 2a 3a

X_1  1 -1  1
X_2  2  . -1
X_3  1  1  1
```
"""


@testset "create character tables" begin
  @testset "library of character tables" begin
    @test character_table("J5") === nothing
    @test_throws ErrorException character_table("A5mod2", 3)
  end

  @test character_table(alternating_group(5), 2) === nothing
end

@testset "create class functions" begin
  g = symmetric_group(3)
  tbl = character_table(g)
  n = number_of_conjugacy_classes(tbl)
  triv = trivial_character(tbl)
  for X in [tbl, g]
    @test triv == Oscar.class_function(X, [1 for i in 1:n])
    @test triv == Oscar.class_function(X, GapObj(triv))
    @test triv == Oscar.class_function(X, GapObj([1 for i in 1:n]))
    @test_throws ArgumentError Oscar.class_function(X, GapObj([true]))
    @test_throws ErrorException Oscar.class_function(X, GapObj([1 for i in 1:2*n]))
  end
end

@testset "access fields in character tables" begin
  # table without group
  t = character_table("A5")
  @test GapObj(t) === t.GAPTable
  @test characteristic(t) == t.characteristic
  @test_throws UndefRefError t.group
  @test_throws UndefRefError t.isomorphism

  # table with `GAPGroup` group
  g = symmetric_group(4)
  t = character_table(g)
  @test GapObj(t) === t.GAPTable
  @test characteristic(t) == t.characteristic
  @test group(t) === t.group === g
  @test Oscar.isomorphism_to_GAP_group(t) === t.isomorphism

  # table with `FinGenAbGroup` group
  g = abelian_group([2, 4])
  t = character_table(g)
  @test GapObj(t) === t.GAPTable
  @test characteristic(t) == t.characteristic
  @test group(t) === t.group === g
  @test Oscar.isomorphism_to_GAP_group(t) === t.isomorphism
end

@testset "attributes of character tables" begin
  ordtbl = character_table("A5")
  modtbl = mod(ordtbl, 2)

  @test modtbl === rem(ordtbl, 2)
  @test modtbl === ordtbl % 2

  @test characteristic(ordtbl) == 0
  @test characteristic(modtbl) == 2
  @test character_parameters(ordtbl) == [[1, 1, 1, 1, 1], [[3, 1, 1], '+'], [[3, 1, 1], '-'], [2, 1, 1, 1], [2, 2, 1]]
  @test character_parameters(modtbl) === nothing
  @test class_lengths(ordtbl) == [1, 15, 20, 12, 12]
  @test class_lengths(modtbl) == [1, 20, 12, 12]
  @test class_names(ordtbl) == ["1a", "2a", "3a", "5a", "5b"]
  @test class_names(modtbl) == ["1a", "3a", "5a", "5b"]
  @test class_parameters(ordtbl) == [[1, 1, 1, 1, 1], [2, 2, 1], [3, 1, 1], [[5], '+'], [[5], '-']]
  @test class_parameters(modtbl) == [[1, 1, 1, 1, 1], [3, 1, 1], [[5], '+'], [[5], '-']]
  @test decomposition_matrix(modtbl) == matrix(ZZ, [1 0 0 0; 1 0 1 0; 1 1 0 0; 0 0 0 1; 1 1 1 0])
  @test identifier(ordtbl) == "A5"
  @test identifier(modtbl) == "A5mod2"
  @test !is_duplicate_table(ordtbl)
  @test maxes(ordtbl) == ["a4", "D10", "S3"]
  @test maxes(modtbl) === nothing
  @test "D10" in names_of_fusion_sources(ordtbl)
  @test order(ordtbl) == 60
  @test order(modtbl) == 60
  @test ordinary_table(modtbl) === ordtbl
  @test_throws ArgumentError ordinary_table(ordtbl)
  @test orders_centralizers(ordtbl) == [60, 4, 3, 5, 5]
  @test orders_centralizers(modtbl) == [60, 3, 5, 5]
  @test orders_class_representatives(ordtbl) == [1, 2, 3, 5, 5]
  @test orders_class_representatives(modtbl) == [1, 3, 5, 5]
  @test trivial_character(ordtbl)[1] == 1
  @test trivial_character(modtbl)[1] == 1
  @test [power_map(ordtbl, 3, i) for i in 1:5] == [1, 2, 1, 5, 4]
end

@testset "p-blocks" begin
  tbl = character_table("A5")
  p = 2
  bl = block_distribution(tbl, p)
  @test bl == Dict(:block => [1, 1, 1, 2, 1], :defect => [2, 0])
  @test bl == block_distribution(tbl, ZZ(p))
  @test_throws ArgumentError block_distribution(mod(tbl, p), p)

  G = alternating_group(5)
  tbl = character_table(G)
  bl = block_distribution(tbl, p)
  defects = bl[:defect]
  @test [order(defect_group(tbl, p, b)[1]) for b in 1:length(defects)] ==
        [p^d for d in defects]
end

@testset "characters" begin
  g = symmetric_group(4)
  t = character_table(g)
  @test nrows(t) == 5
  @test ncols(t) == 5
  @test_throws ErrorException t[6]
  tr = trivial_character(t)
  @test parent(tr) === t
  @test tr in t
  @test tr != t[1]
  @test tr == t[5]
  @test tr == t[end]
  @test tr == trivial_character(g)
  @test !is_faithful(tr)
  @test 2 * tr == tr + tr
  @test ZZ(2) * tr == tr + tr
  @test tr^2 == tr
  @test tr^ZZ(2) == tr
  re = regular_character(g)
  @test coordinates(re) == degree.(t)
  re = regular_character(t)
  @test re == permutation_character(g, trivial_subgroup(g)[1])
  @test coordinates(re) == degree.(t)
  lin = @inferred linear_characters(g)
  @test length(lin) == 2
  lin = @inferred linear_characters(t)
  @test length(lin) == 2
  chi = t[2]
  @test chi isa Oscar.GAPGroupClassFunction
  @test chi[4] == t[2,4]
  @test [chi[i] for i in 1:5] == values(chi)
  @test [chi[nam] for nam in class_names(t)] == values(chi)
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
  @test sort!([order(kernel(chi)[1]) for chi in t]) == [1, 1, 4, 12, 24]
  @test sort!([order(center(chi)[1]) for chi in t]) == [1, 1, 4, 24, 24]
  @test allunique(t)

  @test all(chi -> chi * chi == tensor_product(chi, chi), t)

  scp = scalar_product(t[1], t[1])
  @test scp == 1
  @test scp isa QQFieldElem
  for T in [ZZRingElem, QQFieldElem, Int64, QQAbFieldElem]
    scpT = scalar_product(T, t[1],t[1])
    @test scpT == scp
    @test scpT isa T
  end
  scp = scalar_product(t[1], t[2])
  @test scp == 0
  @test scp isa QQFieldElem
  res = coordinates(chi)
  @test res == [0, 0, 1, 0, 0]
  @test res isa Vector{QQFieldElem}
  @test coordinates(Int, chi) isa Vector{Int}

  orders = orders_class_representatives(t)
  pos = findfirst(==(4), orders)
  chi = filter(x -> x[1] == 2, collect(t))[1]
  ev = multiplicities_eigenvalues(chi, pos)
  @test ev == [0, 1, 0, 1]
  @test ev isa Vector{Int}
  ev = multiplicities_eigenvalues(ZZRingElem, chi, pos)
  @test ev == [0, 1, 0, 1]
  @test ev isa Vector{ZZRingElem}

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
  @test mod(t, 2) === nothing

  g = general_linear_group(2, 3)
  h = derived_subgroup(g)[1]
  p = 2
  t = mod(character_table(h), p)
  chi = t[2]
  x = gen(g, 1)
  @test ! (x in h)
  psi = chi^x
  @test all(y -> mod(order(y), p) == 0 || chi(x*y*x^-1) == psi(y), collect(h))

  t = character_table("A5")
  chi = t[4]
  @test values(central_character(chi)) == [1, 0, 5, -3, -3]
end

@testset "Galois conjugacy of characters" begin
  g = alternating_group(4)
  t = character_table(g)
  @test all(x -> conj(x) == QQAbAutomorphism(5)(x), t)
  @test all(x -> x == QQAbAutomorphism(4)(x), t)

  t = character_table("A5")
  sums = [galois_orbit_sum(x) for x in t]
  rep = galois_representative_and_multiplicity(sums[2])
  @test_throws ArgumentError galois_representative_and_multiplicity(0*t[1])
  @test rep[1] == t[2] && rep[3] == 1
  degrees = [degree_of_character_field(x) for x in t]
  @test degrees == [degree(character_field(x)[1]) for x in t]
  @test degrees == [1, 2, 2, 1, 1]
  @test all(i -> sums[i][1] == t[i][1] * degrees[i], 1:length(sums))
  @test all(x -> degree_of_character_field(x) == 1, sums)
  m = mod(t, 2)
  sums = [galois_orbit_sum(x) for x in m]
  degrees = [degree_of_character_field(x) for x in m]
  @test degrees == [degree(character_field(x)[1]) for x in m]
  @test degrees == [1, 2, 2, 1]
  @test all(i -> sums[i][1] == m[i][1] * degrees[i], 1:length(sums))
  @test all(x -> degree_of_character_field(x) == 1, sums)

  # irreducibles not all rational but all are defined over the prime field
  m = character_table("L3(2)", 2)
  @test all(x -> x == galois_orbit_sum(x), m)
end

@testset "induction and restriction of characters" begin
  g = symmetric_group(4)
  h = derived_subgroup(g)[1]
  t = character_table(g)

  # `induced_cyclic`: no group is needed
  indcyc = @inferred induced_cyclic(t)
  @test sort!([degree(chi) for chi in indcyc]) == [6, 8, 12, 12, 24]
  @test all(x -> scalar_product(trivial_character(t), x) == 1, indcyc)
  @test indcyc == @inferred induced_cyclic(t, 1:nrows(t))

  # `induce` for character tables with groups
  ind = [chi^t for chi in character_table(h)]
  @test sort!([degree(chi) for chi in ind]) == [2, 2, 2, 6]
  @test [scalar_product(trivial_character(t), x) for x in ind] == [1, 0, 0, 0]

  # `induce` for character tables without groups (may fail)
  ta5 = character_table("A5")
  ta6 = character_table("A6")
  ta7 = character_table("A7")
  ind = [chi^ta7 for chi in ta6]
  @test_throws ArgumentError trivial_character(ta5)^ta7

  # `restrict` for character tables with groups
  rest = [restrict(chi, character_table(h)) for chi in t]
  @test [degree(chi) for chi in rest] == [degree(chi) for chi in t]

  # `restrict` for character tables without groups (may fail)
  rest = [restrict(chi, ta6) for chi in ta7]
  @test_throws ArgumentError restrict(ta7[2], ta5)
end

@testset "natural characters" begin
  G = symmetric_group(4)
  chi = Oscar.natural_character(G)
  @test degree(chi) == 4
  psi = chi
  @test scalar_product(chi, psi) == 2
  @test scalar_product(chi, chi) == 2
  @test scalar_product(chi, trivial_character(G)) == 1

  K, a = cyclotomic_field(5, "a")
  L, b = cyclotomic_field(3, "b")
  F, sqrt3 = quadratic_field(3)

  inputs = [
    [ matrix(ZZ, [0 1 0; -1 0 0; 0 0 -1]) ],
    [ matrix(QQ, [0 1 0; -1 0 0; 0 0 -1]) ],
    [ matrix(K, [a 0; 0 a]) ],
    [ matrix(L, 2, 2, [b, 0, -b - 1, 1]), matrix(L, 2, 2, [1, b + 1, 0, b]) ],
    [ matrix(F, [1 0 0 ; 0 -1 0 ; 0 0 -1]),
      matrix(F, [1//2 -sqrt3//2 0 ; sqrt3//2 1//2 0 ; 0 0 1]) ],
  ]

  @testset "... over ring $(base_ring(mats[1]))" for mats in inputs
    G = matrix_group(mats)
    chi = Oscar.natural_character(G)
    @test degree(chi) == degree(G)
    @test chi == natural_character(hom(G, G, gens(G)))
  end

  G = symmetric_group(3)
  @test values(natural_character(hom(G, G, gens(G)))) == [3, 1, 0]

  K, a = cyclotomic_field(8, "a")
  o = one(K)
  z = zero(K)
  G = general_linear_group(2, 3)
  @test values(natural_character(hom(G, G, gens(G)))) ==
        [QQAbFieldElem(x, 8) for x in [2*o, -2*o, z, -a^3-a, a^3+a, z]]

  G = small_group(4, 1)  # pc group
  @test_throws MethodError natural_character(G)
  @test_throws ArgumentError natural_character(hom(G, G, gens(G)))
end

@testset "class fusions" begin
  subtbl = character_table("A5")
  tbl = character_table("A6")
  fus1 = possible_class_fusions(subtbl, tbl)
  fus2 = possible_class_fusions(subtbl, tbl, decompose = false)
  @test fus1 == fus2
  fus3 = possible_class_fusions(subtbl, tbl, fusionmap = fus1[1])
  @test length(fus3) == 1 && fus3[1] == fus1[1]
  @test approximate_class_fusion(subtbl, tbl) == [1, 2, [3, 4], [6, 7], [6, 7]]
  @test approximate_class_fusion(tbl, subtbl) == []
end

@testset "normal subgroups" begin
  tbl = character_table("2.A5.2")
  @test class_positions_of_center(tbl) == [1, 2]
  @test class_positions_of_derived_subgroup(tbl) == [1, 2, 3, 4, 5, 6, 7]
  @test class_positions_of_pcore(tbl, 2) == [1, 2]
  @test class_positions_of_pcore(tbl, 3) == [1]
  @test class_positions_of_solvable_residuum(tbl) == [1, 2, 3, 4, 5, 6, 7]
end

@testset "character tables of factor groups" begin
  g = symmetric_group(4)
  gtbl = character_table(g)
  n = class_positions_of_pcore(gtbl, 2)
  ftbl, fus = quo(gtbl, n)
  f = group(ftbl)
  @test conjugacy_classes(f) == conjugacy_classes(ftbl)
  @test fus == known_class_fusion(gtbl, ftbl)[2]
  @test all(x -> restrict(x, gtbl) in gtbl, ftbl)
end

@testset "character fields of ordinary characters" begin
  tbl = character_table("A5")
  degrees = [degree_of_character_field(chi) for chi in tbl]
  @test degrees == [degree(character_field(chi)[1]) for chi in tbl]
  @test degrees == [1, 2, 2, 1, 1]
  @test character_field(tbl[2])[1] === character_field(tbl[3])[1] # caching
  @test characteristic(character_field(tbl[1])[1]) == 0
  deg = degree_of_character_field(collect(tbl))
  @test deg == degree(character_field(collect(tbl))[1])
  @test deg == 2

  for id in [ "C5", "A5" ]   # cyclotomic and non-cyclotomic number fields
    for chi in character_table(id)
      F1, phi = character_field(chi)
      F2, _ = number_field(QQ, chi)
      F3, _ = QQ[chi]
      @test degree(F1) == degree(F2)
      @test degree(F1) == degree(F3)
      @test conductor(chi) == conductor(phi)
      @test conductor(Int, chi) isa Int
      for i in 1:length(chi)
        x = chi[i]
        xF = preimage(phi, x)
        @test parent(xF) == F1
        @test phi(xF) == x
      end
    end
  end

  # embeddings of abelian number fields into the abelian closure
  # - for cyclotomic character fields
  t = character_table("A4")
  chi = t[2]
  F, emb = character_field(chi)
  @test Hecke.is_cyclotomic_type(F)[1]
  for elm in [gen(F), one(F)]
    img = emb(elm)
    @test preimage(emb, img) == elm
    @test has_preimage_with_preimage(emb, img) == (true, elm)
    z5 = gen(parent(img))(5)
    @test has_preimage_with_preimage(emb, z5)[1] == false
    @test_throws ErrorException preimage(emb, z5)
  end

  # - for non-cyclotomic character fields
  t = character_table("A5")
  chi = t[2]
  F, emb = character_field(chi)
  @test ! Hecke.is_cyclotomic_type(F)[1]
  for elm in [gen(F), one(F)]
    img = emb(elm)
    @test preimage(emb, img) == elm
    @test has_preimage_with_preimage(emb, img) == (true, elm)
    z5 = gen(parent(img))(5)
    @test has_preimage_with_preimage(emb, z5)[1] == false
    @test_throws ErrorException preimage(emb, z5)
  end
end

@testset "character fields of Brauer characters" begin
  ordtbl = character_table("A5")
  modtbl = mod(ordtbl, 2)
  @test character_field(modtbl[2])[1] === character_field(modtbl[3])[1] # caching
  @test [order_field_of_definition(chi) for chi in modtbl] == [2, 4, 4, 2]
  @test [order(character_field(chi)[1]) for chi in modtbl] == [2, 4, 4, 2]
  @test order_field_of_definition(Int, modtbl[1]) isa Int
  @test_throws ArgumentError order_field_of_definition(ordtbl[1])
  deg = degree_of_character_field(collect(modtbl))
  @test deg == degree(character_field(collect(modtbl))[1])
  @test deg == 2
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
  @test sort!(map(schur_index, collect(t))) == append!(repeat([1],15), repeat([2],4))

  # For a character table without group, we can fail.
  t = character_table("S6")
  for chi in t
    try
      schur_index(chi)
    catch(e)
      msg = sprint(showerror, e)
      @test msg == "ArgumentError: cannot determine the Schur index with the currently used criteria"
    end
  end

  # The function is defined only for ordinary characters.
  t = mod(t, 2)
  @test_throws ArgumentError schur_index(t[1])
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

  t = character_table_complex_reflection_group(4, 1, 3)
  @test length(character_parameters(t)) == length(t)
  @test length(class_parameters(t)) == length(t)
  @test get_attribute(t, :type) == (4, 1, 3)
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

@testset "character functions for FinGenAbGroup" begin
  @testset for para in [ Int[], [2, 3, 4], [2, 4] ]
    G1 = abelian_group(FinGenAbGroup, para)
    iso = isomorphism(PcGroup, G1)
    G2 = codomain(iso)
    n = Int(order(G1))

    # create character tables
    tbl1 = character_table(G1)
    tbl2 = character_table(G2)
    @test length(tbl1) == length(tbl2) == n

    # characters
    tr1 = trivial_character(tbl1)
    tr2 = trivial_character(tbl2)
    @test values(tr1) == values(tr2)
    @test tr1 in tbl1
    chi = tbl1[n]
    @test chi[n] == tbl1[n,n]
    @test [chi[i] for i in 1:n] == values(chi)
    @test [2*chi[i] for i in 1:n] == values(chi + chi)
    @test [chi[i]^2 for i in 1:n] == values(chi * chi)
    @test [chi[i]^2 for i in 1:n] == values(chi^2)
    @test [-chi[i] for i in 1:n] == values(-chi)
    @test [0*chi[i] for i in 1:n] == values(chi-chi)
    @test [1 for i in 1:n] == values(one(chi))
    @test [0 for i in 1:n] == values(zero(chi))
    @test degree(chi) == chi[1]
    @test length(chi) == n
    @test 1 in chi
    @test [chi(representative(c)) for c in conjugacy_classes(tbl1)] == values(chi)
    @test all(x -> degree(x) == 1, linear_characters(G1))
    @test all(is_irreducible, tbl1)
    @test all(chi -> order(center(chi)[1]) == n, tbl1)
    @test all(chi -> is_subgroup(kernel(chi)[1], G1)[1], tbl1)

    @test all(chi -> chi * chi == tensor_product(chi, chi), tbl1)

    scp = scalar_product(tbl1[1], tbl1[1])
    @test scp == 1
    @test scp isa QQFieldElem
    for T in [ZZRingElem, QQFieldElem, Int64, QQAbFieldElem]
      scpT = scalar_product(T, tbl1[1],tbl1[1])
      @test scpT == scp
      @test scpT isa T
    end
    if n > 1
      scp = scalar_product(tbl1[1], tbl1[2])
      @test scp == 0
      @test scp isa QQFieldElem
    end

    # induced characters
    H = sylow_subgroup(G1, 3)[1]
    subtbl = character_table(H)
    ind = [chi^tbl1 for chi in subtbl]
    @test all(chi -> degree(chi) == index(G1, H), ind)
    @test ind == [chi^G1 for chi in subtbl]

    # restricted characters
    rest = [restrict(psi, subtbl) for psi in ind]
    @test all(chi -> degree(chi) == index(G1, H), rest)
    @test rest == [restrict(psi, H) for psi in ind]

    # conjugate characters
    H, _ = pcore(G1, 2)
    th = character_table(H)
    chi = th[1]
    for x in gens(G1)
      psi = chi^x
      @test psi == chi
    end

    # a corresponding Brauer table
    tmod2 = mod(tbl1, 2)
    @test tmod2 isa Oscar.GAPGroupCharacterTable
    n = ncols(tmod2)
    chi = tmod2[1]
    @test degree(chi) == chi[1]
    @test length(chi) == n
    @test chi[n] == tmod2[1,n]
    tr = trivial_character(tmod2)
    @test ! (-1 in tr)

    mat = decomposition_matrix(tmod2)
    @test nrows(mat) == nrows(tbl1)
    @test ncols(mat) == nrows(tmod2)

    # character fields
    for chi in tbl1
      F, phi = character_field(chi)
      for x in chi
        xF = preimage(phi, x)
        @test parent(xF) == F
        @test phi(xF) == x
      end
    end
  end
end

@testset "read off group properties from character tables" begin
  t = character_table("A5")
  @test ! is_abelian(t)
  @test is_almost_simple(t)
  @test ! is_cyclic(t)
  @test ! is_elementary_abelian(t)
  @test ! is_nilpotent(t)
  @test is_perfect(t)
  @test is_quasisimple(t)
  @test is_simple(t)
  @test ! is_solvable(t)
  @test ! is_sporadic_simple(t)
  @test ! is_supersolvable(t)
end

@testset "action  on class functions" begin
  g = symmetric_group(4)
  h = pcore(g,2)[1]
  t = character_table(h)
  chi = t[2]
  @test order(stabilizer(g, chi)[1]) == 8
  chi = t[3]  # now the action on classes is already stored
  @test order(stabilizer(g, chi)[1]) == 8

  h = sylow_subgroup(g,2)[1]
  t = character_table(h)
  chi = t[2]
  @test_throws ArgumentError stabilizer(g, chi)
end
