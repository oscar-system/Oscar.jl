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
        "DoubleAndHyperComplexes",      # TODO: remove from this list. `MethodError: no method matching simplify(::SubquoModule{MPolyQuoRingElem{MPolyDecRingElem{FqMPolyRingElem, AbstractAlgebra.Generic.MPoly{FqMPolyRingElem}}}})`
        "ExteriorAlgebra",              # TODO: remove from this list. `undefined binding 'exterior_algebra' in `@docs` block in src/NoncommutativeAlgebra/PBWAlgebras/quotients.md:40-42` and `Error During Test at /home/runner/work/Oscar.jl/Oscar.jl/test/Rings/PBWAlgebraQuo.jl:41`
        "GaloisGrp",                    # TODO: remove from this list. `no docs found for 'fixed_field(C::Oscar.GaloisGrp.GaloisCtx, s::Vector{PermGroup})' in `@docs` block in src/NumberTheory/galois.md:275-278`
        "GModule",                      # TODO: remove from this list. many doctest failures and `MethodError: no method matching numerator(::QQPolyRingElem)`
        "InvariantTheory",              # TODO: remove from this list. `undefined binding 'linearly_reductive_group' in `@docs` block in src/InvariantTheory/reductive_groups.md:53-55` and more docs errors`
        "ModStd",                       # TODO: remove from this list. `MethodError: no method matching monomial(::QQMPolyRing, ::Vector{Int64})` and many similar errors
        "Rings",                        # TODO: remove from this list. used by Rings/binomial_ideals.jl, see https://github.com/oscar-system/Oscar.jl/blob/13282dfd07b6aee58e433a45353f48261cda787b/src/Oscar.jl#L268
        "Schemes",                      # TODO: untangle src/AlgebraicGeometry/Schemes/ and experimental/Schemes/
      ],
    )
  elseif command == "cleanup"
    move_all("experimental2", "experimental")
  else
    error("Error, unknown command: $command. Possible commands are 'init' and 'cleanup'")
  end
end

main(ARGS)
