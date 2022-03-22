#############################################################################
##
##  Provide functionality for displaying labelled matrices,
##  that is, two-dimensional arrays of strings
##  together with row labels, columns labels, header, and footer
##
##  Typical examples are character tables of groups,
##  where the character names and class names appear as row and column labels,
##  respectively, the name of the table is shown in the header,
##  and a legend of irrationalities info in the footer.
##
##  There is support for the output of \LaTeX format strings
##  as well as for column aligned text in a Julia session.
##
##  The code was inspired by the ideas underlying
##  [GAP's Browse package](http://www.math.rwth-aachen.de/~Browse/),
##  and by Jean Michel's implementation of character table display in his
##  [Julia package Gapjm.jl](https://github.com/jmichel7/Gapjm.jl).
##
##  Some open items (to be addressed if there is interest):
##  Support left/centered/right alignment of columns and of the whole table.

export labelled_matrix_formatted


# `lpad` is wrong for example for strings containing overlined characters,
# such as "A̅".
# (This had been observed already by Jean Michel.
# See https://github.com/JuliaLang/julia/pull/39044,
# a fix seems to become available in Julia 1.7.)
mylpad(s::String,n::Int) = " "^(n-textwidth(s))*s


# for conversions from \LaTeX format to unicode format,
# deals with subscripts and superscripts in braces `{}`,
# character values, character names
const superscript = Dict(zip("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
const subscript = Dict(zip("-0123456789", "₋₀₁₂₃₄₅₆₇₈₉"))

function replace_TeX(str::String)
  if Oscar.is_unicode_allowed()
    str = replace(str, "\\chi" => "χ")
    str = replace(str, "\\phi" => "φ")
    str = replace(str, "\\varphi" => "φ")
    str = replace(str, "\\zeta" => "ζ")
    str = replace(str, r"\^\{[-0123456789]*\}" => (s -> map(x->superscript[x], s[3:end-1])))
    str = replace(str, r"_\{[-0123456789]*\}" => (s -> map(x->subscript[x], s[3:end-1])))
    str = replace(str, r"\\overline\{(?<nam>[A-Z]*)\}" => s"\g<nam>\u0305")
  else
    str = replace(str, "\\chi" => "X")
    str = replace(str, "\\phi" => "Y")
    str = replace(str, "\\varphi" => "Y")
    str = replace(str, "\\zeta" => "z")
    str = replace(str, r"\^\{[-0123456789]*\}" => (s -> "^"*s[3:end-1]))
    str = replace(str, r"_\{[-0123456789]*\}" => (s -> "_"*s[3:end-1]))
    str = replace(str, r"\\overline\{(?<nam>[A-Z]*)\}" => s"/\g<nam>")
  end
  return str
end


"""
    labelled_matrix_formatted(io::IO, mat::Matrix{String})

Write a formatted version of `mat` to `io`.
The following attributes of `io` are supported.

- `:TeX`:
  Produce a LaTeX format `array` if the value is `true`,
  otherwise produce column aligned text format (if `Oscar.is_unicode_allowed()`
  returns `true` then with unicode characters and sub-/superscripts).

- `:subset_row`:
  array of positions of those rows of `mat` that shall be shown,

- `:subset_col`:
  array of positions of those columns of `mat` that shall be shown,

- `:header`:
  array of strings to be shown above `mat`,

- `:corner`:
  matrix (or vector) of strings to be shown on the left of (and aligned to)
  the column labels and above (and aligned to) the row labels

- `:labels_row`:
  matrix (or vector) of strings to be shown on the left of and aligned to
  the rows of `mat`,
  (the `i`-th row corresponds to the `i`-th row of `mat`,
  also if `:subset_row` restricts the shown rows),

- `:labels_col`:
  matrix (or vector) of strings to be shown above and aligned to
  the columns of `mat`
  (the `j`-th column corresponds to the `j`-th column of `mat`,
  also if `:subset_col` restricts the shown columns),

- `:footer`:
  array of strings to be shown below `mat`,

- `:separators_row`:
  array of numbers `i` such that a horizontal line shall be drawn
  below the `i`-th shown row
  (enter `0` for a line above the first row),

- `:separators_col`:
  array of numbers `j` such that a vertical line shall be drawn
  instead of a blank behind the `j`-th shown column
  (enter `0` for a line in front of the first column),

- `:portions_row`:
  array of numbers of rows after which a new labelled table shall be started;
  the default is to have just one portion.

- `:portions_col`:
  array of numbers of column after which a new labelled table shall be started;
  the default is to have just one portion in the `:TeX` case,
  and to create portions according to the screen width otherwise,

## Examples
```jldoctest
julia> m = 3; n = 4;  mat = Array{String}(undef, m, n);

julia> for i in 1:m for j in 1:n mat[i,j] = string( (i,j) ); end; end

julia> mat
3×4 Matrix{String}:
 "(1, 1)"  "(1, 2)"  "(1, 3)"  "(1, 4)"
 "(2, 1)"  "(2, 2)"  "(2, 3)"  "(2, 4)"
 "(3, 1)"  "(3, 2)"  "(3, 3)"  "(3, 4)"

julia> io = IOBuffer();   # simply using `stdout` does not work in doctests

julia> ioc = IOContext(io,
              :header => ["", "a labelled matrix", ""],
              :labels_row => [string(i) for i in 1:m],
              :labels_col => [string(j) for j in 1:n],
              :separators_row => [0],
              :separators_col => [0],
              :footer => ["", "with footer"],
             );

julia> Oscar.with_unicode() do
         labelled_matrix_formatted(ioc, mat)
       end;

julia> print(String(take!(io)))

a labelled matrix

 │     1      2      3      4
─┼───────────────────────────
1│(1, 1) (1, 2) (1, 3) (1, 4)
2│(2, 1) (2, 2) (2, 3) (2, 4)
3│(3, 1) (3, 2) (3, 3) (3, 4)

with footer

julia> ioc = IOContext(io,
              :labels_row => [string(i) for i in 1:m],
              :labels_col => [string(j) for j in 1:n],
              :separators_row => [0],
              :separators_col => [0],
              :portions_row => [2,1],
              :portions_col => [2,2],
             );

julia> Oscar.with_unicode() do
         labelled_matrix_formatted(ioc, mat)
       end;

julia> print(String(take!(io)))
 │     1      2
─┼─────────────
1│(1, 1) (1, 2)
2│(2, 1) (2, 2)

 │     1      2
─┼─────────────
3│(3, 1) (3, 2)

 │     3      4
─┼─────────────
1│(1, 3) (1, 4)
2│(2, 3) (2, 4)

 │     3      4
─┼─────────────
3│(3, 3) (3, 4)

```
"""
function labelled_matrix_formatted(io::IO, mat::Matrix{String})
    TeX = get(io, :TeX, false)

    if Oscar.is_unicode_allowed()
      pipe = "\u2502"
      horz_bar = "\u2500"
      cross = "\u253C"
    else
      pipe = "|"
      horz_bar = "-"
      cross = "+"
    end

    subset_row = get(io, :subset_row, 1:size(mat, 1))
    subset_col = get(io, :subset_col, 1:size(mat, 2))
    m2 = size(subset_row, 1)
    n2 = size(subset_col, 1)

    # Get the three surrounding matrices.
    # If some of them were given as vectors then turn them into matrices.
    labels_col = get(io, :labels_col, fill("", 0, n2))
    labels_row = get(io, :labels_row, fill("", m2, 0))

    if length(size(labels_row)) == 1
      # one column of row labels
      labels_row = reshape(labels_row, length(labels_row), 1)
    end
    if length(size(labels_col)) == 1
      # one row of column labels
      labels_col = reshape(labels_col, 1, length(labels_col))
    end

    m1 = size(labels_col, 1)
    n1 = size(labels_row, 2)
    corner = get(io, :corner, fill("", m1, n1))
    if length(size(corner)) == 1
      # one column of labels above the row labels
      corner = reshape(corner, length(corner), 1)
    end

    @assert size(corner, 1) == m1 "row numbers in corner and column labels differ"
    @assert size(corner, 2) == n1 "column numbers in corner and row labels differ"
    @assert size(labels_row, 1) == m2 "row numbers in row labels and matrix differ"
    @assert size(labels_col, 2) == n2 "column numbers in column labels and matrix differ"

    m = m1 + m2

    header = get(io, :header, [])
    footer = get(io, :footer, [])

    # Concatenate corner part and row labels (restricted to the relevant rows).
    leftpart = vcat(corner, labels_row[subset_row, 1:size(labels_row, 2)])

    # Concatenate column labels (restricted to the relevant columns)
    # and `mat` (restricted to the relevant rows and columns).
    rightpart = vcat(labels_col[1:size(labels_col, 1), subset_col],
                     mat[subset_row, subset_col])

    if ! TeX
      # If no `TeX` output is required then
      # replace markup in the matrices, ...
      leftpart = map(replace_TeX, leftpart)
      rightpart = map(replace_TeX, rightpart)

      header = map(replace_TeX, header)
      footer = map(replace_TeX, footer)

      # ... compute columns widths of the corner entries and row labels ...
      leftwidth = map(textwidth, leftpart)
      leftrows = [leftwidth[i,:] for i in 1:m]
      widths_leftpart = map(max, leftrows...)

      # ... and of the column labels and matrix entries.
      rightwidth = map(textwidth, rightpart)
      rightrows = [rightwidth[i,:] for i in 1:m]
      widths_rightpart = map(max, rightrows...)

      # Compute the width of the row labels part.
      leftwidthsum = sum(widths_leftpart) + length(widths_leftpart) - 1
      if leftwidthsum == -1
        leftwidthsum = 0
      end
    end

    # Print the header (a vector of strings).
    for line in header
      write(io, line, "\n")
    end

    separators_col = get(io, :separators_col, [])
    separators_row = map(x -> x+m1, get(io, :separators_row, []))

    # Determine the column portions.
    portions_col = get(io, :portions_col, [])
    if portions_col == []
      # Set a default.
      if TeX
        # Do not split into column portions.
        portions_col = [n2]
      else
        # Distribute the matrix into column portions,
        # according to the screen width.
        screen_width = displaysize(io)[2]
        if leftwidthsum >= screen_width
          write(io, "(row label part is too wide for the screen)\n")
          return
        end

        num = 0
        width = leftwidthsum
        for w in widths_rightpart
          width = width + w + 1
          if width >= screen_width
            if num == 0
              # This column will be too wide
              # but we have to take at least one.
              push!(portions_col, 1)
              width = leftwidthsum
            else
              push!(portions_col, num)
              num = 1
              width = leftwidthsum + w + 1
            end
          else
            num = num+1
          end
        end
        if num != 0
          push!(portions_col, num)
        end
      end
    end

    # We do not distribute the table into row portions
    # if this is not wanted.
    portions_row = get(io, :portions_row, [m2])

    # Run over the column portions.
    cpos = 0
    for ci in 1:length(portions_col)
      cnum = portions_col[ci]
      crange = (cpos + 1):(cpos + cnum)
      cpos = cpos + cnum

      if TeX
        cpattern = 0 in separators_col ? "|" : ""
        for j in crange
          cpattern = cpattern * (j in separators_col ? "r|" : "r")
        end
      else
        if 0 in separators_col
          cprefix = pipe
        elseif n1 == 0
          cprefix = ""
        else
          cprefix = " "
        end
        cpattern = String[]
        for j in crange
          push!(cpattern, j in separators_col ? pipe : " ")
        end
        if ! (crange[end] in separators_col)
          cpattern[end] = ""
        end
      end

      # For each column portion, run over the row portions.
      rpos = 0
      for ri in 1:length(portions_row)
        rnum = portions_row[ri]
        # Print first the column labels and then the row portion.
        rrange = vcat(1:m1, (rpos + m1 + 1):(rpos + m1 + rnum))
        rpos = rpos + rnum
        if TeX
          write(io, "\\begin{array}{" * repeat("r", n1) * cpattern * "}", "\n")
        end

        # Compute the horizontal separator line.
        if TeX
          hline = "\\hline"
        else
          hline = repeat(horz_bar, leftwidthsum)
          if 0 in separators_col
            hline = hline * cross
          elseif n1 != 0
            hline = hline * horz_bar
          end
          for jj in 1:length(crange)
            j = crange[jj]
            hline = hline * repeat(horz_bar, widths_rightpart[j])
            if jj < length(crange) || j in separators_col
              if j in separators_col
                hline = hline * cross
              else
                hline = hline * horz_bar
              end
            end
          end
        end

        # Print the matrices.
        for i in rrange
          if 0 in separators_row && ri == 1 && m1 == 0 && i == rrange[1]
            write(io, hline, "\n")
          end
          if TeX
            write(io, join(leftpart[i,:], " & "))
            if n1 > 0
              write(io, " & ")
            end
            write(io, join(rightpart[i,:][crange], " & "), " \\\\\n")
          else
            write(io, join(map(mylpad, leftpart[i,:], widths_leftpart), " "),
                      cprefix,
                      join(join.(collect(zip(map(mylpad, rightpart[i,:][crange], widths_rightpart[crange]), cpattern)))), "\n")
          end
          if i in separators_row
            write(io, hline, "\n")
          end
        end
        if TeX
          write(io, "\\end{array}\n")
        end
        if ci != length(portions_col) || ri != length(portions_row)
          write(io, "\n")
        end
      end
    end

    # footer (vector of strings)
    for line in footer
      write(io, line, "\n")
    end
end

