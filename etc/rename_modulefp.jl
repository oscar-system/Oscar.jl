# this can be run from the command line with:
# > julia --project=. etc/rename_modulefp.jl

function replace_in_file(filepath::String, old::String, new::String)
  pattern = Regex("\\b$old\\b")
  content = read(filepath, String)

  if occursin(pattern, content)
    new_content = replace(content, pattern => new)
    open(filepath, "w") do io
      write(io, new_content)
    end
  end
end

function gather_source_files(path::String)
  result = String[]

  # Check if the user has git.
  if !isnothing(Sys.which("git"))
    output = read(
      `git ls-files "$path/*.jl" "$path/**/*.jl" "$path/*.md" "$path/**/*.md"`, String
    )
    result = string.(split(output))
  else
    queue = [path]
    while length(queue) > 0
      current_folder = pop!(queue)
      next_batch = readdir(current_folder; join=true)
      for elem in next_batch
        if isfile(elem) && (endswith(elem, ".jl") || endswith(elem, ".md"))
          push!(result, elem)
        elseif isdir(elem)
          push!(queue, elem)
        end
      end
    end
  end

  return result
end

replacements = [
  "ModuleFP" => "OFPModule",
  "ModuleFPHom" => "OFPModuleHom",
  "ModuleFPElem" => "OFPModuleElem",
  "AdmissibleModuleFPRing" => "AdmissibleOFPModuleRing",
  "AdmissibleModuleFPRingElem" => "AdmissibleOFPModuleRingElem",
  "ModuleFP_dec" => "OFPModule_dec",
  "ModuleFPElem_dec" => "OFPModuleElem_dec",
  "ModuleFPHomDummy" => "OFPModuleHomDummy",
]

file = @__FILE__
oscardir = replace(file, "etc/rename_modulefp.jl" => "")

files = [
  gather_source_files(joinpath(oscardir, "src"))
  gather_source_files(joinpath(oscardir, "docs"))
  gather_source_files(joinpath(oscardir, "test"))
  gather_source_files(joinpath(oscardir, "experimental"))
]

# Rename all occurrences of old typename to new typename.
for file in files
  for (old, new) in replacements
    replace_in_file(file, old, new)
  end
end
