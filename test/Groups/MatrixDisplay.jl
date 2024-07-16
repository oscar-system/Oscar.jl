Oscar.@_AuxDocTest "show and print formatted matrices", (fix = false),
raw"""
create the test matrix

```jldoctest labeled_matrix_formatted.test
julia> using Oscar

julia> m = 2; n = 2;  mat = Array{String}(undef, m, n);

julia> for i in 1:m for j in 1:n mat[i,j] = string(i)*"/"*string(j); end; end;
```

text format

no labels

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout), mat)
1/1 1/2
2/1 2/2
```

with header

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :header => ["header", ""]), mat)
header

1/1 1/2
2/1 2/2
```

with footer

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :footer => ["footer"]), mat)
1/1 1/2
2/1 2/2
footer
```

with row labels as vector

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m]), mat)
1: 1/1 1/2
2: 2/1 2/2
```

with row labels as matrix

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_row => reshape([string(i)*":" for i in 1:m], m, 1)), mat)
1: 1/1 1/2
2: 2/1 2/2
```

with too wide row labels

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_row => reshape(["_"^30*string(i)*":" for i in 1:m], m, 1),
         :displaysize => (52, 30)), mat)
(row label part is too wide for the screen)
```

with column labels as vector

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_col => [string(j) for j in 1:n]), mat)
  1   2
1/1 1/2
2/1 2/2
```

with column labels as matrix

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_col => reshape([string(j) for j in 1:n], 1, n)), mat)
  1   2
1/1 1/2
2/1 2/2
```

with row and column labels

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m],
         :labels_col => [string(j) for j in 1:n]), mat)
     1   2
1: 1/1 1/2
2: 2/1 2/2
```

with row and column labels and corner as vector

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m],
         :labels_col => [string(j) for j in 1:n],
         :corner => ["(i,j)"]), mat)
(i,j)   1   2
   1: 1/1 1/2
   2: 2/1 2/2
```

with row and column labels and corner as matrix

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m],
         :labels_col => [string(j) for j in 1:n],
         :corner => reshape(["(i,j)"], 1, 1)), mat)
(i,j)   1   2
   1: 1/1 1/2
   2: 2/1 2/2
```

with a too wide column

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m],
         :labels_col => ["_"^30*string(j) for j in 1:n],
         :corner => reshape(["(i,j)"], 1, 1),
         :displaysize => (52, 30)), mat)
(i,j) ______________________________1
   1:                             1/1
   2:                             2/1

(i,j) ______________________________2
   1:                             1/2
   2:                             2/2
```

with row separators but without column labels

```jldoctest labeled_matrix_formatted.test
julia> Oscar.with_unicode() do
         labeled_matrix_formatted(IOContext(stdout,
           :separators_row => [0, 1, 2]), mat)
       end
───────
1/1 1/2
───────
2/1 2/2
───────
julia> labeled_matrix_formatted(IOContext(stdout,
         :separators_row => [0, 1, 2]), mat)
-------
1/1 1/2
-------
2/1 2/2
-------
```

with row separators and column labels

```jldoctest labeled_matrix_formatted.test
julia> Oscar.with_unicode() do
         labeled_matrix_formatted(IOContext(stdout,
         :labels_col => [string(j) for j in 1:n],
         :separators_row => [0,1,2]), mat)
       end
  1   2
───────
1/1 1/2
───────
2/1 2/2
───────
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_col => [string(j) for j in 1:n],
         :separators_row => [0,1,2]), mat)
  1   2
-------
1/1 1/2
-------
2/1 2/2
-------
```

with column separators but without row labels

```jldoctest labeled_matrix_formatted.test
julia> Oscar.with_unicode() do
         labeled_matrix_formatted(IOContext(stdout,
         :separators_col => [0, 1, 2]), mat)
       end
│1/1│1/2│
│2/1│2/2│
julia> labeled_matrix_formatted(IOContext(stdout,
         :separators_col => [0, 1, 2]), mat)
|1/1|1/2|
|2/1|2/2|
```

with column separators and row labels

```jldoctest labeled_matrix_formatted.test
julia> Oscar.with_unicode() do
         labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m],
         :separators_col => [0, 1, 2]), mat)
       end
1:│1/1│1/2│
2:│2/1│2/2│
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m],
         :separators_col => [0, 1, 2]), mat)
1:|1/1|1/2|
2:|2/1|2/2|
```

with row and column labels and separators

```jldoctest labeled_matrix_formatted.test
julia> Oscar.with_unicode() do
         labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m],
         :labels_col => [string(j) for j in 1:n],
         :separators_row => [0, 1, 2],
         :separators_col => [0, 1, 2]), mat)
       end
  │  1│  2│
──┼───┼───┼
1:│1/1│1/2│
──┼───┼───┼
2:│2/1│2/2│
──┼───┼───┼
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m],
         :labels_col => [string(j) for j in 1:n],
         :separators_row => [0, 1, 2],
         :separators_col => [0, 1, 2]), mat)
  |  1|  2|
--+---+---+
1:|1/1|1/2|
--+---+---+
2:|2/1|2/2|
--+---+---+
```

with row and column portions

```jldoctest labeled_matrix_formatted.test
julia> Oscar.with_unicode() do
         labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m],
         :labels_col => [string(j) for j in 1:n],
         :separators_row => [0],
         :separators_col => [0],
         :portions_row => [1,1],
         :portions_col => [1,1]), mat)
       end
  │  1
──┼───
1:│1/1

  │  1
──┼───
2:│2/1

  │  2
──┼───
1:│1/2

  │  2
──┼───
2:│2/2
julia> labeled_matrix_formatted(IOContext(stdout,
         :labels_row => [string(i)*":" for i in 1:m],
         :labels_col => [string(j) for j in 1:n],
         :separators_row => [0],
         :separators_col => [0],
         :portions_row => [1,1],
         :portions_col => [1,1]), mat)
  |  1
--+---
1:|1/1

  |  1
--+---
2:|2/1

  |  2
--+---
1:|1/2

  |  2
--+---
2:|2/2
```

LaTeX format

```jldoctest labeled_matrix_formatted.test
julia> m = 2; n = 2;  mat = Array{String}(undef, m, n);

julia> for i in 1:m for j in 1:n mat[i,j] = string(i)*"/"*string(j); end; end;
```

no labels

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true), mat)
$\begin{array}{rr}
1/1 & 1/2 \\
2/1 & 2/2 \\
\end{array}
$
```

with header

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :header => ["header", ""]), mat)
header

$\begin{array}{rr}
1/1 & 1/2 \\
2/1 & 2/2 \\
\end{array}
$
```

with footer

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :footer => ["footer"]), mat)
$\begin{array}{rr}
1/1 & 1/2 \\
2/1 & 2/2 \\
\end{array}

\begin{array}{l}
footer \\
\end{array}
$
```

with row labels as vector

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :labels_row => [string(i)*":" for i in 1:m]), mat)
$\begin{array}{rrr}
1: & 1/1 & 1/2 \\
2: & 2/1 & 2/2 \\
\end{array}
$
```

with row labels as matrix

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :labels_row => reshape([string(i)*":" for i in 1:m], m, 1)), mat)
$\begin{array}{rrr}
1: & 1/1 & 1/2 \\
2: & 2/1 & 2/2 \\
\end{array}
$
```

with column labels as vector

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
        :labels_col => [string(j) for j in 1:n]), mat)
$\begin{array}{rr}
1 & 2 \\
1/1 & 1/2 \\
2/1 & 2/2 \\
\end{array}
$
```

with column labels as matrix

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :labels_col => reshape([string(j) for j in 1:n], 1, n)), mat)
$\begin{array}{rr}
1 & 2 \\
1/1 & 1/2 \\
2/1 & 2/2 \\
\end{array}
$
```

with row and column labels

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :labels_row => [string(i)*":" for i in 1:m],
         :labels_col => [string(j) for j in 1:n]), mat)
$\begin{array}{rrr}
 & 1 & 2 \\
1: & 1/1 & 1/2 \\
2: & 2/1 & 2/2 \\
\end{array}
$
```

with row separators but without column labels

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :separators_row => [0, 1, 2]), mat)
$\begin{array}{rr}
\hline
1/1 & 1/2 \\
\hline
2/1 & 2/2 \\
\hline
\end{array}
$
```

with row separators and column labels

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :labels_col => [string(j) for j in 1:n],
         :separators_row => [0,1,2]), mat)
$\begin{array}{rr}
1 & 2 \\
\hline
1/1 & 1/2 \\
\hline
2/1 & 2/2 \\
\hline
\end{array}
$
```

with column separators but without row labels

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :separators_col => [0, 1, 2]), mat)
$\begin{array}{|r|r|}
1/1 & 1/2 \\
2/1 & 2/2 \\
\end{array}
$
```

with column separators and row labels

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :labels_row => [string(i)*":" for i in 1:m],
         :separators_col => [0, 1, 2]), mat)
$\begin{array}{r|r|r|}
1: & 1/1 & 1/2 \\
2: & 2/1 & 2/2 \\
\end{array}
$
```

with row and column labels and separators

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :labels_row => [string(i)*":" for i in 1:m],
         :labels_col => [string(j) for j in 1:n],
         :separators_row => [0, 1, 2],
         :separators_col => [0, 1, 2]), mat)
$\begin{array}{r|r|r|}
 & 1 & 2 \\
\hline
1: & 1/1 & 1/2 \\
\hline
2: & 2/1 & 2/2 \\
\hline
\end{array}
$
```

with row and column portions

```jldoctest labeled_matrix_formatted.test
julia> labeled_matrix_formatted(IOContext(stdout, :TeX => true,
         :labels_row => [string(i) for i in 1:m],
         :labels_col => [string(j) for j in 1:n],
         :separators_row => [0],
         :separators_col => [0],
         :portions_row => [1,1],
         :portions_col => [1,1]), mat)
$\begin{array}{r|r}
 & 1 \\
\hline
1 & 1/1 \\
\end{array}

\begin{array}{r|r}
 & 1 \\
\hline
2 & 2/1 \\
\end{array}

\begin{array}{r|r}
 & 2 \\
\hline
1 & 1/2 \\
\end{array}

\begin{array}{r|r}
 & 2 \\
\hline
2 & 2/2 \\
\end{array}
$
```
"""
