# this can be run from the command line with:
# > julia -e 'include("etc/check_meta.jl")'
#
# note: the @main function requires julia 1.11 or newer

function (@main)(_)
  issues = String[]
  meta = """```@meta
  CurrentModule = Oscar
  CollapsedDocStrings = true
  DocTestSetup = Oscar.doctestsetup()
  ```
  """
  meta_n_lines = count('\n', meta)
  oscardir = normpath(@__DIR__, "..")
  for (dir, _, files) in walkdir(oscardir)
    contains(dir, "docs/src") || continue
    for file in files
      endswith(file, ".md") || continue
      path = joinpath(dir, file)
      if read(`head -n$(meta_n_lines) $(path)`, String) != meta
        push!(issues, path)
      end
    end
  end
  isempty(issues) && return
  println("The standard @meta block is:")
  println(meta)
  println("The below markdown files do not have a standard @meta block:")
  for file in issues
    println(replace(file, oscardir => ""))
  end
  exit(-2)
end
