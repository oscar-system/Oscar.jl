using Pkg
Pkg.activate(; temp=true)
Pkg.add("Test")
Pkg.add(name="JuliaFormatter", version="1")
using Test
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
  "experimental/BasisLieHighestWeight",
  "experimental/ExperimentalTemplate",
  "experimental/ExteriorAlgebra",
  "experimental/LieAlgebras",
  "experimental/LinearQuotients",
]
# a few files for one reason or another may need to be skipped
ignore = [
   "experimental/InvariantTheory/src/InvariantTheory.jl",
]

result = true
for fname in files_and_dirs_to_be_formatted
  global result
  result &= format(fname; ignore)
end

run(`git status`)

exit(!result)
