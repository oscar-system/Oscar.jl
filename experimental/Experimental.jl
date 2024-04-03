# We read experimental and filter out all packages that follow our desired
# scheme. Remember those packages to avoid doing this all over again for docs
# and test.
# We don't want to interfere with existing stuff in experimental though.
const expdir = joinpath(@__DIR__, "../experimental")
const oldexppkgs = [
  "ExteriorAlgebra",
  "GModule",
  "ModStd",
  "Rings",
  "Schemes",
  "FTheoryTools" # Must be loaded after the schemes.
]
# DEVELOPER OPTION:
# If an experimental package A depends on another experimental package B, one
# can add `"B", "A"` to `orderedpkgs` below to ensure that B is loaded before A.
# DO NOT USE THIS UNLESS YOU KNOW THE CONSEQUENCES.
# For more background, see https://github.com/oscar-system/Oscar.jl/issues/2300.
const orderedpkgs = [
  "LieAlgebras",
  "BasisLieHighestWeight",   # needs code from LieAlgebras
]
exppkgs = filter(x->isdir(joinpath(expdir, x)) && !(x in oldexppkgs) && !(x in orderedpkgs), readdir(expdir))
append!(exppkgs, orderedpkgs)

# force trigger recompile when folder changes
include_dependency(".")

# setup for the NoExperimental CI job
# For local testing, run `ln -s NoExperimental_whitelist_.jl experimental/NoExperimental_whitelist.jl` to initialize
# and `rm experimental/NoExperimental_whitelist.jl` for cleanup
isfile(joinpath(expdir, "NoExperimental_whitelist_.jl")) || error("experimental/NoExperimental_whitelist_.jl is missing")
if islink(joinpath(expdir, "NoExperimental_whitelist.jl"))
  include(joinpath(expdir, "NoExperimental_whitelist.jl"))
  filter!(in(whitelist), exppkgs)
  filter!(in(whitelist), oldexppkgs)
end

# Error if something is incomplete in experimental
for pkg in exppkgs
  if !isfile(joinpath(expdir, pkg, "src", "$pkg.jl"))
    error("experimental/$pkg is incomplete: $pkg/src/$pkg.jl missing. See the documentation at https://docs.oscar-system.org/dev/Experimental/intro/ for details.")
  end
  path = joinpath(expdir, pkg, "test")
  if !isdir(path) || length(filter(endswith(".jl"), readdir(path))) == 0
    error("experimental/$pkg is incomplete: $pkg/test/ missing or empty. See the documentation at https://docs.oscar-system.org/dev/Experimental/intro/ for details.")
  end
  # Load the package
  include(joinpath(expdir, pkg, "src", "$pkg.jl"))
end


# Force some structure for `oldexppkgs`
for pkg in oldexppkgs
  if !isfile(joinpath(expdir, pkg, "$pkg.jl"))
    error("experimental/$pkg is incomplete: $pkg/$pkg.jl missing. Please fix this or remove $pkg from `oldexppkgs`.")
  end
  # Load the package
  include(joinpath(expdir, pkg, "$pkg.jl"))
end
