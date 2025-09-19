using Pkg
Pkg.activate(; temp=true)
Pkg.add(name="JuliaFormatter", version="1")

using JuliaFormatter

# Disclaimer:
# The goal should be to work with 
# https://github.com/julia-actions/julia-format
# . However, currently too many files have broken formatting and since there
# are many ongoing pull requests, we will need to extend proper formatting to
# the entire codebase in a step by step fashion.
#
# In case you format some code, also add the commit hash to
# .git-blame-ignore-revs for `git blame` to ignore these commits.

oscardir = joinpath(@__DIR__, "..")

files_and_dirs_to_be_formatted = [
  "src/InvariantTheory",
  "src/PolyhedralGeometry",
  "test/PolyhedralGeometry",
  "src/LieTheory",
  "test/LieTheory",
  "src/aliases.jl",
  "src/AlgebraicGeometry/ToricVarieties",
  "test/AlgebraicGeometry/ToricVarieties",
  "experimental/BasisLieHighestWeight",
  "experimental/ExperimentalTemplate",
  "experimental/ExteriorAlgebra",
  "experimental/LieAlgebras",
  "experimental/LinearQuotients",
  "experimental/FTheoryTools",
]
# a few files for one reason or another may need to be skipped
ignore = [
  "experimental/InvariantTheory/src/InvariantTheory.jl",
]

for fname in files_and_dirs_to_be_formatted
  while !format(fname; ignore)
    # invoke format until it reports no more changes, to work around
    # <https://github.com/domluna/JuliaFormatter.jl/issues/897>
  end
end
