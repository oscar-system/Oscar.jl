@testset "labelled_matrix_formatted, text format" begin
  m = 2; n = 2;  mat = Array{String}(undef, m, n);
  for i in 1:m for j in 1:n mat[i,j] = string(i)*"/"*string(j); end; end

  old_allow = Oscar.is_unicode_allowed()

  io = IOBuffer();

  # no labels
  ioc = IOContext(io,
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "1/1 1/2\n2/1 2/2\n"

  # with header
  ioc = IOContext(io,
          :header => ["header", ""],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "header\n\n1/1 1/2\n2/1 2/2\n"

  # with footer
  ioc = IOContext(io,
          :footer => ["footer"],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "1/1 1/2\n2/1 2/2\nfooter\n"

  # with row labels as vector
  ioc = IOContext(io,
          :labels_row => [string(i)*":" for i in 1:m],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "1: 1/1 1/2\n2: 2/1 2/2\n"

  # with row labels as matrix
  ioc = IOContext(io,
          :labels_row => reshape([string(i)*":" for i in 1:m], m, 1),
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "1: 1/1 1/2\n2: 2/1 2/2\n"

  # with too wide row labels
  ioc = IOContext(io,
          :labels_row => reshape(["_"^30*string(i)*":" for i in 1:m], m, 1),
          :displaysize => (52, 30),
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "(row label part is too wide for the screen)\n"

  # with column labels as vector
  ioc = IOContext(io,
          :labels_col => [string(j) for j in 1:n],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "  1   2\n1/1 1/2\n2/1 2/2\n"

  # with column labels as matrix
  ioc = IOContext(io,
          :labels_col => reshape([string(j) for j in 1:n], 1, n),
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "  1   2\n1/1 1/2\n2/1 2/2\n"

  # with row and column labels
  ioc = IOContext(io,
          :labels_row => [string(i)*":" for i in 1:m],
          :labels_col => [string(j) for j in 1:n],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "     1   2\n1: 1/1 1/2\n2: 2/1 2/2\n"

  # with row and column labels and corner as vector
  ioc = IOContext(io,
          :labels_row => [string(i)*":" for i in 1:m],
          :labels_col => [string(j) for j in 1:n],
          :corner => ["(i,j)"],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "(i,j)   1   2\n   1: 1/1 1/2\n   2: 2/1 2/2\n"

  # with row and column labels and corner as matrix
  ioc = IOContext(io,
          :labels_row => [string(i)*":" for i in 1:m],
          :labels_col => [string(j) for j in 1:n],
          :corner => reshape(["(i,j)"], 1, 1),
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "(i,j)   1   2\n   1: 1/1 1/2\n   2: 2/1 2/2\n"

  # with a too wide column
  ioc = IOContext(io,
          :labels_row => [string(i)*":" for i in 1:m],
          :labels_col => ["_"^30*string(j) for j in 1:n],
          :corner => reshape(["(i,j)"], 1, 1),
          :displaysize => (52, 30),
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "(i,j) ______________________________1\n   1:                             1/1\n   2:                             2/1\n\n(i,j) ______________________________2\n   1:                             1/2\n   2:                             2/2\n"

  # with row separators but without column labels
  ioc = IOContext(io,
          :separators_row => [0, 1, 2],
        );
  Oscar.allow_unicode(true)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "───────\n1/1 1/2\n───────\n2/1 2/2\n───────\n"
  Oscar.allow_unicode(false)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "-------\n1/1 1/2\n-------\n2/1 2/2\n-------\n"

  # with row separators and column labels
  ioc = IOContext(io,
          :labels_col => [string(j) for j in 1:n],
          :separators_row => [0,1,2],
        );
  Oscar.allow_unicode(true)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "  1   2\n───────\n1/1 1/2\n───────\n2/1 2/2\n───────\n"
  Oscar.allow_unicode(false)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "  1   2\n-------\n1/1 1/2\n-------\n2/1 2/2\n-------\n"

  # with column separators but without row labels
  ioc = IOContext(io,
          :separators_col => [0, 1, 2],
        );
  Oscar.allow_unicode(true)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "│1/1│1/2│\n│2/1│2/2│\n"
  Oscar.allow_unicode(false)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "|1/1|1/2|\n|2/1|2/2|\n"

  # with column separators and row labels
  ioc = IOContext(io,
          :labels_row => [string(i)*":" for i in 1:m],
          :separators_col => [0, 1, 2],
        );
  Oscar.allow_unicode(true)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "1:│1/1│1/2│\n2:│2/1│2/2│\n"
  Oscar.allow_unicode(false)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "1:|1/1|1/2|\n2:|2/1|2/2|\n"

  # with row and column labels and separators
  ioc = IOContext(io,
          :labels_row => [string(i)*":" for i in 1:m],
          :labels_col => [string(j) for j in 1:n],
          :separators_row => [0, 1, 2],
          :separators_col => [0, 1, 2],
        );
  Oscar.allow_unicode(true)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "  │  1│  2│\n──┼───┼───┼\n1:│1/1│1/2│\n──┼───┼───┼\n2:│2/1│2/2│\n──┼───┼───┼\n"
  Oscar.allow_unicode(false)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "  |  1|  2|\n--+---+---+\n1:|1/1|1/2|\n--+---+---+\n2:|2/1|2/2|\n--+---+---+\n"

  # with row and column portions
  ioc = IOContext(io,
          :labels_row => [string(i)*":" for i in 1:m],
          :labels_col => [string(j) for j in 1:n],
          :separators_row => [0],
          :separators_col => [0],
          :portions_row => [1,1],
          :portions_col => [1,1],
        );
  Oscar.allow_unicode(true)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "  │  1\n──┼───\n1:│1/1\n\n  │  1\n──┼───\n2:│2/1\n\n  │  2\n──┼───\n1:│1/2\n\n  │  2\n──┼───\n2:│2/2\n"
  Oscar.allow_unicode(false)
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "  |  1\n--+---\n1:|1/1\n\n  |  1\n--+---\n2:|2/1\n\n  |  2\n--+---\n1:|1/2\n\n  |  2\n--+---\n2:|2/2\n"

  Oscar.allow_unicode(old_allow)
end

@testset "labelled_matrix_formatted, LaTeX format" begin
  m = 2; n = 2;  mat = Array{String}(undef, m, n);
  for i in 1:m for j in 1:n mat[i,j] = string(i)*"/"*string(j); end; end

  io = IOBuffer();

  # no labels
  ioc = IOContext(io, :TeX => true,
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{rr}\n1/1 & 1/2 \\\\\n2/1 & 2/2 \\\\\n\\end{array}\n"

  # with header
  ioc = IOContext(io, :TeX => true,
          :header => ["header", ""],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "header\n\n\\begin{array}{rr}\n1/1 & 1/2 \\\\\n2/1 & 2/2 \\\\\n\\end{array}\n"

  # with footer
  ioc = IOContext(io, :TeX => true,
          :footer => ["footer"],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{rr}\n1/1 & 1/2 \\\\\n2/1 & 2/2 \\\\\n\\end{array}\nfooter\n"

  # with row labels as vector
  ioc = IOContext(io, :TeX => true,
          :labels_row => [string(i)*":" for i in 1:m],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{rrr}\n1: & 1/1 & 1/2 \\\\\n2: & 2/1 & 2/2 \\\\\n\\end{array}\n"

  # with row labels as matrix
  ioc = IOContext(io, :TeX => true,
          :labels_row => reshape([string(i)*":" for i in 1:m], m, 1),
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{rrr}\n1: & 1/1 & 1/2 \\\\\n2: & 2/1 & 2/2 \\\\\n\\end{array}\n"

  # with column labels as vector
  ioc = IOContext(io, :TeX => true,
          :labels_col => [string(j) for j in 1:n],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{rr}\n1 & 2 \\\\\n1/1 & 1/2 \\\\\n2/1 & 2/2 \\\\\n\\end{array}\n"

  # with column labels as matrix
  ioc = IOContext(io, :TeX => true,
          :labels_col => reshape([string(j) for j in 1:n], 1, n),
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{rr}\n1 & 2 \\\\\n1/1 & 1/2 \\\\\n2/1 & 2/2 \\\\\n\\end{array}\n"

  # with row and column labels
  ioc = IOContext(io, :TeX => true,
          :labels_row => [string(i)*":" for i in 1:m],
          :labels_col => [string(j) for j in 1:n],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{rrr}\n & 1 & 2 \\\\\n1: & 1/1 & 1/2 \\\\\n2: & 2/1 & 2/2 \\\\\n\\end{array}\n"

  # with row separators but without column labels
  ioc = IOContext(io, :TeX => true,
          :separators_row => [0, 1, 2],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{rr}\n\\hline\n1/1 & 1/2 \\\\\n\\hline\n2/1 & 2/2 \\\\\n\\hline\n\\end{array}\n"

  # with row separators and column labels
  ioc = IOContext(io, :TeX => true,
          :labels_col => [string(j) for j in 1:n],
          :separators_row => [0,1,2],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{rr}\n1 & 2 \\\\\n\\hline\n1/1 & 1/2 \\\\\n\\hline\n2/1 & 2/2 \\\\\n\\hline\n\\end{array}\n"

  # with column separators but without row labels
  ioc = IOContext(io, :TeX => true,
          :separators_col => [0, 1, 2],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{|r|r|}\n1/1 & 1/2 \\\\\n2/1 & 2/2 \\\\\n\\end{array}\n"

  # with column separators and row labels
  ioc = IOContext(io, :TeX => true,
          :labels_row => [string(i)*":" for i in 1:m],
          :separators_col => [0, 1, 2],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{r|r|r|}\n1: & 1/1 & 1/2 \\\\\n2: & 2/1 & 2/2 \\\\\n\\end{array}\n"

  # with row and column labels and separators
  ioc = IOContext(io, :TeX => true,
          :labels_row => [string(i)*":" for i in 1:m],
          :labels_col => [string(j) for j in 1:n],
          :separators_row => [0, 1, 2],
          :separators_col => [0, 1, 2],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{r|r|r|}\n & 1 & 2 \\\\\n\\hline\n1: & 1/1 & 1/2 \\\\\n\\hline\n2: & 2/1 & 2/2 \\\\\n\\hline\n\\end{array}\n"

  # with row and column portions
  ioc = IOContext(io, :TeX => true,
          :labels_row => [string(i) for i in 1:m],
          :labels_col => [string(j) for j in 1:n],
          :separators_row => [0],
          :separators_col => [0],
          :portions_row => [1,1],
          :portions_col => [1,1],
        );
  labelled_matrix_formatted(ioc, mat)
  @test String(take!(io)) == "\\begin{array}{r|r}\n & 1 \\\\\n\\hline\n1 & 1/1 \\\\\n\\end{array}\n\n\\begin{array}{r|r}\n & 1 \\\\\n\\hline\n2 & 2/1 \\\\\n\\end{array}\n\n\\begin{array}{r|r}\n & 2 \\\\\n\\hline\n1 & 1/2 \\\\\n\\end{array}\n\n\\begin{array}{r|r}\n & 2 \\\\\n\\hline\n2 & 2/2 \\\\\n\\end{array}\n"
end
