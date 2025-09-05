const expdir = joinpath(@__DIR__, "../experimental")

# DEVELOPER OPTION:
# If an experimental package A depends on another experimental package B, one
# can add `"B", "A"` to `orderedpkgs` below to ensure that B is loaded before A.
# DO NOT USE THIS UNLESS YOU KNOW THE CONSEQUENCES.
# For more background, see https://github.com/oscar-system/Oscar.jl/issues/2300.
const orderedpkgs = [
  "LieAlgebras",
  "BasisLieHighestWeight",   # needs code from LieAlgebras
  "AlgebraicShifting",       # Needs code from Lie Algebras (`isomorphism(PermGroup, ::WeylGroup)` specifically)
  "SetPartitions",
  "PartitionedPermutations", # needs code from SetPartitions
  "Schemes",
  "FTheoryTools",            # must be loaded after Schemes and LieAlgebras
  "IntersectionTheory",      # must be loaded after Schemes
  "Parallel",
  "AlgebraicStatistics" # must be loaded after Parallel
]
const exppkgs = filter(x->isdir(joinpath(expdir, x)) && !(x in orderedpkgs), readdir(expdir))
append!(exppkgs, orderedpkgs)

# force trigger recompile when folder changes
include_dependency(".")

# setup for the NoExperimental CI job
# For local testing, run `ln -s NoExperimental_whitelist_.jl experimental/NoExperimental_whitelist.jl` to initialize
# and `rm experimental/NoExperimental_whitelist.jl` for cleanup
isfile(joinpath(expdir, "NoExperimental_whitelist_.jl")) || error("experimental/NoExperimental_whitelist_.jl is missing")
if islink(joinpath(expdir, "NoExperimental_whitelist.jl"))
  include(joinpath(expdir, "NoExperimental_whitelist.jl"))
  issubset(whitelist, exppkgs) || error("experimental/NoExperimental_whitelist.jl contains unknown packages")
  filter!(in(whitelist), exppkgs)
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

# We modify the documentation of the experimental part to attach a warning to
# every exported function that this function is part of experimental.
# Furthermore we give a link for users to read up on what this entails.
#
# Note that there are functions in the docs of experimental that are not
# exported. These then also do not get the warning attached.
using Markdown
const warnexp = Markdown.parse(raw"""
!!! warning "Experimental"
    This function is part of the experimental code in Oscar. Please read
    [here](https://docs.oscar-system.org/v1/Experimental/intro/) for more
    details.
""")
for name in names(Oscar)
  if isdefined(Oscar, name)
    md = Base.Docs.doc(getfield(Oscar, name))
    # Loop over all definitions of a function
    for entry in md.meta[:results]
      # Test whether function was defined in experimental
      if startswith(entry.data[:path], joinpath(Oscar.oscardir, "experimental"))
        append!(entry.object.content[1].content, warnexp.content)
      end
    end
  end
end
