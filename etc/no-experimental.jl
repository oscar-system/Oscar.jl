function move_all(source, target)
  isdir(source) || error("'$(source)' directory does not exist")
  mkpath(target)
  for obj in readdir(source)
    mv(joinpath(source, obj), joinpath(target, obj))
  end
  rm(source)
end

function move_all_but_list(source, target, list)
  isdir(source) || error("'$(source)' directory does not exist")
  issubset(list, readdir(source)) || error(
    "Not all files in list are in source directory. You might want to update the list"
  )
  mkpath(target)
  for obj in readdir(source)
    if !(obj in list)
      mv(joinpath(source, obj), joinpath(target, obj))
    end
  end
end

function main(args)
  length(args) == 1 || error("No command given. Possible commands are 'init' and 'cleanup'")
  command = only(args)
  if command == "init"
    move_all_but_list(
      "experimental",
      "experimental2",
      [
        "Experimental.jl",              # needs to stay here forever
        "BasisLieHighestWeight",        # TODO: remove from this list
        "DoubleAndHyperComplexes",      # TODO: remove from this list
        "EnumerativeCombinatorics",     # TODO: remove from this list
        "ExperimentalTemplate",         # TODO: remove from this list
        "ExteriorAlgebra",              # TODO: remove from this list
        "FTheoryTools",                 # TODO: remove from this list
        "GaloisGrp",                    # TODO: remove from this list
        "GITFans",                      # TODO: remove from this list
        "GModule",                      # TODO: remove from this list
        "IntersectionTheory",           # TODO: remove from this list
        "InvariantTheory",              # TODO: remove from this list
        "LieAlgebras",                  # TODO: remove from this list
        "LinearQuotients",              # TODO: remove from this list
        "MatrixGroups",                 # TODO: remove from this list
        "MatroidRealizationSpaces",     # TODO: remove from this list
        "ModStd",                       # TODO: remove from this list
        "OrthogonalDiscriminants",      # TODO: remove from this list
        "QuadFormAndIsom",              # TODO: remove from this list
        "Rings",                        # TODO: remove from this list
        "Schemes",                      # TODO: remove from this list
        "SetPartitions",                # TODO: remove from this list
        "StandardFiniteFields",         # TODO: remove from this list
        "SymmetricIntersections",       # TODO: remove from this list
      ],
    )
  elseif command == "cleanup"
    move_all("experimental2", "experimental")
  else
    error("Error, unknown command: $command. Possible commands are 'init' and 'cleanup'")
  end
end

main(ARGS)
