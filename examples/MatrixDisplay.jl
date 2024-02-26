#############################################################################
##
##  functionality for displaying labeled matrices,
##  that is, 2-dimensional arrays of strings
##  together with row labels, columns labels, header, and footer
##
##  Typical examples are character tables of groups.
##
##  There is support for the output of \LaTeX format strings
##  as well as for column aligned text in a Julia session.
##
##  The code was inspired by the ideas underlying
##  [GAP's Browse package](http://www.math.rwth-aachen.de/~Browse/),
##  and by Jean Michel's implementation of character table display in his
##  [Julia package Gapjm.jl](https://github.com/jmichel7/Gapjm.jl).
##
#T TODO:
#T - Support left/centered/right alignment of columns and of the whole table.
#T - Concerning irrational values in character tables:
#T   in `matrix_of_strings`,
#T   identify also unique Galois conjugates (mark with superscript *);
#T   show both values in the footer.
#T   (For elements in quadratic fields, show also an expression in terms of
#T   square roots.)


# `lpad` is wrong for example for strings containing overlined characters,
# such as "A̅".
# (This had been observed already by Jean Michel.)
mylpad(s::String,n::Int) = " "^(n-textwidth(s))*s

# for conversions from \LaTeX format to unicode format,
# deals with subscripts and superscripts in braces `{}`,
# character values, character names
const sup = Dict(zip("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
const sub = Dict(zip("-0123456789", "₋₀₁₂₃₄₅₆₇₈₉"))

function replace_TeX(str::String)
  str = replace(str, "\\chi" => "χ")
  str = replace(str, "\\phi" => "φ")
  str = replace(str, "\\varphi" => "φ")
  str = replace(str, "\\zeta" => "ζ")
  str = replace(str, r"\^\{[-0123456789]*\}" => (s -> map(x->sup[x], s[3:end-1])))
  str = replace(str, r"_\{[-0123456789]*\}" => (s -> map(x->sub[x], s[3:end-1])))
  str = replace(str, r"\\overline\{(?<nam>[A-Z]*)\}" => s"\g<nam>\u0305")
  return str
end

"""
    labeled_matrix_formatted(io::IO, mat::Matrix{String})

Output a formatted version of `mat` to `io`.
The following attributes of `io` are supported.

- `:TeX`:
  Produce a LaTeX format `array` if the value is `true`,
  otherwise produce column aligned text format with unicode characters and
  sub-/superscripts.

- `:subset_row`:
  array of positions of those rows of `mat` that shall be shown,

- `:subset_col`:
  array of positions of those columns of `mat` that shall be shown,

- `:header`:
  array of strings to be shown above `mat`,

- `:corner`:
  matrix (or vector) of strings to be shown on the left of (and aligned to)
  the the column labels and above (and aligned to) the row labels

- `:labels_col`:
  matrix of strings to be shown above and aligned to the columns of `mat`,

- `:labels_row`:
  matrix (or vector) of strings to be shown on the left of (and aligned to)
  the rows of `mat`,

- `:footer`:
  array of strings to be shown below `mat`,

- `:column_separators`:
  array of column numbers after which a vertical line shall be drawn
  instead of a blank (enter `0` for a line in front of the first column),

- `:row_separators`:
  array of row numbers after which a horizontal line shall be drawn
  (enter `0` for a line above the first row),

- `:column_portions`:
  array of column numbers after which a new labeled table shall be started;
  the default is to have just one portion in the `:TeX` case,
  and to create portions according to the screen width otherwise,

- `:row_portions`:
  array of row numbers after which a new labeled table shall be started;
  the default is to have just one portion.

## Examples:
...
"""
function labeled_matrix_formatted(io::IO, mat::Matrix{String})
    TeX = get(io, :TeX, false)

    subset_row = get(io, :subset_row, 1:size(mat, 1))
    subset_col = get(io, :subset_col, 1:size(mat, 2))
    m2 = size(subset_row, 1)
    n2 = size(subset_col, 1)

    # Get the three surrounding matrices.
    labels_col = get(io, :labels_col, fill("", 0, n2))
    labels_row = get(io, :labels_row, fill("", m2, 0))
    m1 = size(labels_col, 1)
    n1 = size(labels_row, 2)
    corner = get(io, :corner, fill("", m1, n1))
#FIXME:
# If `subset_row`, `subset_col` are given, are `labels_row`, `labels_col`
# to be read relative to the subsets or to `mat`?

    @assert size(corner, 1) == m1 "row numbers in corner and column labels differ"
    @assert size(corner, 2) == n1 "column numbers in corner and row labels differ"
    @assert size(labels_row, 1) == m2 "row numbers in row labels and matrix differ"
    @assert size(labels_col, 2) == n2 "column numbers in column labels and matrix differ"

    m = m1 + m2

    header = get(io, :header, [])
    footer = get(io, :footer, [])

    # Concatenate corner part and row labels.
    leftpart = vcat(corner, labels_row)

    # Concatenate column labels and `mat`.
    rightpart = vcat(labels_col, mat[subset_row, subset_col])

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
    end

    # Print the header (a vector of strings).
    for line in header
      println(io, line)
    end

    column_separators = get(io, :column_separators, [])
    row_separators = map(x -> x+m1, get(io, :row_separators, []))

    # Determine the column portions.
    column_portions = get(io, :column_portions, [])
    if column_portions == []
      # Set a default.
      if TeX
        # Do not split into column portions.
        column_portions = [n2]
      else
        # Distribute the matrix into column portions,
        # according to the screen width.
        screen_width = displaysize(io)[2]
        if leftwidthsum >= screen_width
          println(io, "(row label part is too wide for the screen)")
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
              push!(column_portions, 1)
              num = 0
              width = leftwidthsum
            else
              push!(column_portions, num)
              num = 1
              width = leftwidthsum + w + 1
            end
          else
            num = num+1
          end
        end
        push!(column_portions, num)
      end
    end

    # We do not distribute the table into row portions
    # if this is not wanted.
    row_portions = get(io, :row_portions, [m2])

    # Run over the column portions.
    cpos = 0
    lastportion = (column_portions[end], row_portions[end])
    for cnum in column_portions
      crange = (cpos + 1):(cpos + cnum)
      cpos = cpos + cnum

      if TeX
        cpattern = 0 in column_separators ? "|" : ""
        for j in crange
          cpattern = cpattern * (j in column_separators ? "r|" : "r")
        end
      else
        cprefix = 0 in column_separators ? "\u2502" : " "
        cpattern = String[]
        for j in crange
          push!(cpattern, j in column_separators ? "\u2502" : " ")
        end
      end

      # For each column portion, run over the row portions.
      rpos = 0
      for rnum in row_portions
        # Print first the column labels and then the row portion.
        rrange = vcat(1:m1, (rpos + m1 + 1):(rpos + m1 + rnum))
        rpos = rpos + rnum
        if TeX
          println(io, "\\begin{array}{" * repeat("r", n1) * cpattern * "}")
        end

        # Compute the horizontal separator line.
        if TeX
          hline = "\\hline"
        else
          hline = repeat("\u2500", leftwidthsum)
          if 0 in column_separators
            hline = hline * "\u253C"
          else
            hline = hline * "\u2500"
          end
          for j in crange
            hline = hline * repeat("\u2500", widths_rightpart[j])
            if j in column_separators
              hline = hline * "\u253C"
            else
              hline = hline * "\u2500"
            end
          end
        end

        # Print the matrices.
        for i in rrange
          if TeX
            println(io, join(leftpart[i,:], " & "), " & ", join(rightpart[i,:][crange], " & "), " \\\\")
            if i in row_separators
              println(io, "\\hline")
            end
          else
            println(io, join(map(mylpad, leftpart[i,:], widths_leftpart), " "),
                        cprefix,
                        join(join.(collect(zip(map(mylpad, rightpart[i,:][crange], widths_rightpart[crange]), cpattern)))))
            if i in row_separators
              println(io, hline)
            end
          end
        end
        if TeX
          println(io, "\\end{array}")
        end
        if (cnum, rnum) != lastportion
          println(io, "")
        end
      end
    end

    # footer (vector of strings)
    for line in footer
      println(io, line)
    end
end

