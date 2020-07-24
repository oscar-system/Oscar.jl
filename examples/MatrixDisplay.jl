#############################################################################
##
##  functionality for displaying labelled matrices,
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
#T - Support separating horizontal/vertical lines (for both LaTeX and screen
#T   output).
#T - Distribute the table into column parts if required:
#T   if the `:column_parts` attribute of `io` prescribes an array of
#T   arrays of column indices (most likely a partition of the set of all
#T   column indices) then produce a column block for each array of indices,
#T   each together with the row labels.
#T   If `:column_parts` is not given and `:TeX` is `true` then do not
#T   distribute the table into column blocks.
#T   If `:column_parts` is not given and `:TeX` is not `true` then
#T   use the screensize to break the table into column portions of maximal
#T   width.
#T - Support left/centered/right alignment of columns and of the whole table.
#T - In the case of character tables, support a variant that shows
#T   irrational values as symbols and shows their values in the footer of
#T   the table.

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
  return str
end

"""
    function labelled_matrix_formatted(io::IO, mat::Matrix{String})

Output a formatted version of `mat` to `io`.
Produce a LaTeX format `array` if the `:TeX` attribute of `io` is `true`,
otherwise produce column aligned text format with unicode characters and
sub-/superscripts.
"""
function labelled_matrix_formatted(io::IO, mat::Matrix{String})
    TeX = get(io, :TeX, false)
    m2 = size(mat, 1)
    n2 = size(mat, 2)

    # Get the three surrounding matrices.
    labels_col = get(io, :labels_col, fill("", 0, n2))
    labels_row = get(io, :labels_row, fill("", m2, 0))
    m1 = size(labels_col, 1)
    n1 = size(labels_row, 2)
    corner = get(io, :corner, fill("", m1, n1))

    @assert size(corner, 1) == m1 "row numbers in corner and column labels differ"
    @assert size(corner, 2) == n1 "column numbers in corner and row labels differ"
    @assert size(labels_row, 1) == m2 "row numbers in row labels and matrix differ"
    @assert size(labels_col, 2) == n2 "column numbers in column labels and matrix differ"

    m = m1 + m2

    # Concatenate corner part and row labels.
    leftpart = vcat(corner, labels_row)

    # Concatenate column labels and `mat`.
    rightpart = vcat(labels_col, mat)

    if ! TeX
      # If no `TeX` output is required then replace markup in the matrices.
      leftpart = map(replace_TeX, leftpart)
      rightpart = map(replace_TeX, rightpart)

      # If no `TeX` output is required then compute columns widths
      # of the corner entries and row labels ...
      leftwidth = map(textwidth, leftpart)
      leftrows = [leftwidth[i,:] for i in 1:m]
      widths_leftpart = map(max, leftrows...)

      # ... and of the column labels and matrix entries
      rightwidth = map(textwidth, rightpart)
      rightrows = [rightwidth[i,:] for i in 1:m]
      widths_rightpart = map(max, rightrows...)
    end

    # Print the header (a vector of strings).
    header = get(io, :header, [])
    for line in header
      println(io, line)
    end
    if TeX
      println(io, "\\begin{array}{" * repeat("r", n1+n2) * "}")
    end

    # Print the matrices.
    for i in 1:m
      if TeX
        println(io, join(leftpart[i,:], " & "), " & ", join(rightpart[i,:], " & "), " \\\\")
      else
        println(io, join(map(lpad, leftpart[i,:], widths_leftpart), " "),
                    " ",
                    join(map(lpad, rightpart[i,:], widths_rightpart), " "))
      end
    end
    if TeX
      println(io, "\\end{array}")
    end

    # footer (vector of strings)
    footer = get(io, :footer, [])
    for line in footer
      println(io, line)
    end
end

