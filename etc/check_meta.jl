# this can be run from the command line with:
# > julia -e 'include("etc/check_meta.jl")'
#
# note: the @main function requires julia 1.11 or newer

function (@main)(_)
  issues = String[]
  oscardir = normpath(@__DIR__, "..")
  for (dir, _, files) in walkdir(oscardir)
    contains(dir, "docs/src") || continue
    for file in files
      endswith(file, ".md") || continue
      path = joinpath(dir, file)
      if success(`grep "@meta" $(path)`)
        push!(issues, path)
      end
    end
  end
  isempty(issues) && return
  println("The below markdown files do have a @meta block, but probably don't need it:") # remove this CI check once we actually need @meta blocks
  for file in issues
    println(replace(file, oscardir => ""))
  end
  exit(-2)
end
