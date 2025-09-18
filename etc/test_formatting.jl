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

file = @__FILE__
oscardir = replace(file, "etc/test_formatting.jl" => "")

result = 0

@testset "Formatting" begin
  function _gather_source_files(path::String)
    queue = [path]
    result = String[]
    while length(queue) > 0
      current_folder = pop!(queue)
      next_batch = readdir(current_folder; join=true)
      for elem in next_batch
        if isfile(elem) && endswith(elem, ".jl")
          push!(result, elem)
        elseif isdir(elem)
          push!(queue, elem)
        end
      end
    end
    return result
  end

  # Since right now only very few files are formatted, we also use a whitelist
  # approach.
  enabled = [
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
  skip = [
    "experimental/InvariantTheory/src/InvariantTheory.jl",# the path matches the whitelist entry "src/InvariantTheory"
  ]

  # Collect all code files.
  entire = [
    _gather_source_files(joinpath(oscardir, "etc/"))
    _gather_source_files(joinpath(oscardir, "src/"))
    _gather_source_files(joinpath(oscardir, "experimental/"))
  ]

  failed = String[]

  for file in entire
    is_enabled = !isnothing(findfirst(occursin(file), enabled))
    should_skip = !isnothing(findfirst(occursin(file), skip))
    if is_enabled && !should_skip
      # Do not actually format file, only check whether format is ok.
      res = @test format(file; overwrite=false)
      if res isa Test.Fail
        result = 1
        filename = replace(file, oscardir => "")
        push!(failed, filename)
      else
        filename = replace(file, oscardir => "")
      end
    else
      @test format(file; overwrite=false) skip = true
    end
  end

  # Produce some nice output at the end so users do not have to scroll through
  # entire output.
  if length(failed) > 0
    println("""
The files
$failed
were not properly formatted. You can fix the formatting by running
```
using JuliaFormatter""")
    for fn in failed
      println("format(\"$fn\")")
    end
    println("```\n")
  end
end

exit(result)
