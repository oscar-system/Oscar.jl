# We read experimental and filter out all packages that follow our desired
# scheme. Remember those packages to avoid doing this all over again for docs
# and test.
# We don't want to interfere with existing stuff in experimental though.
const expdir = joinpath(@__DIR__, "../experimental")
const oldexppkgs = [
  "ExteriorAlgebra",
  "GModule",
  "MatrixGroups",
  "ModStd",
  "Rings",
  "Schemes",
]
# DEVELOPER OPTION:
# If an experimental package A depends on another experimental package B, one
# can add `"B", "A"` to `orderedpkgs` below to ensure that B is loaded before A.
# DO NOT USE THIS UNLESS YOU KNOW THE CONSEQUENCES.
# For more background, see https://github.com/oscar-system/Oscar.jl/issues/2300.
const orderedpkgs = [
  "LieAlgebras",
  "BasisLieHighestWeight",   # nees code from LieAlgebras
]
exppkgs = filter(x->isdir(joinpath(expdir, x)) && !(x in oldexppkgs) && !(x in orderedpkgs), readdir(expdir))
append!(exppkgs, orderedpkgs)

# Error if something is incomplete in experimental
for pkg in exppkgs
  if !isfile(joinpath(expdir, pkg, "src", "$pkg.jl"))
    error("experimental/$pkg is incomplete: $pkg/src/$pkg.jl missing. See the documentation at https://docs.oscar-system.org/dev/Experimental/intro/ for details.")
  end
  path = joinpath(expdir, pkg, "test")
  if !isdir(path) || length(filter(endswith(".jl"), readdir(path))) == 0
    error("experimental/$pkg is incomplete: $pkg/test/ missing or empty. See the documentation at https://docs.oscar-system.org/dev/Experimental/intro/ for details.")
  end
end

# force trigger recompile when folder changes
include_dependency(".")

for pkg in Oscar.exppkgs
  include("$pkg/src/$pkg.jl")
end

include("Rings.jl")
include("ModStd.jl")
include("GModule.jl")

include("MatrixGroups/matrix.jl")

include("Schemes/Types.jl")
include("Schemes/CoveredScheme.jl")
include("Schemes/FunctionFields.jl")
include("Schemes/ProjectiveModules.jl")
include("Schemes/SpaceGerms.jl")
include("Schemes/Sheaves.jl")
include("Schemes/IdealSheaves.jl")
include("Schemes/AlgebraicCycles.jl")
include("Schemes/WeilDivisor.jl")
include("Schemes/CoveredProjectiveSchemes.jl")

include("Schemes/SimplifiedAffineScheme.jl")
include("Schemes/CoherentSheaves.jl")
include("Schemes/LazyGluing.jl")
include("Schemes/CartierDivisor.jl")
include("Schemes/Auxiliary.jl")
include("Schemes/BlowupMorphism.jl")
include("Schemes/duValSing.jl")
include("Schemes/elliptic_surface.jl")
include("Schemes/MorphismFromRationalFunctions.jl")

include("Schemes/ToricIdealSheaves/auxiliary.jl")
include("Schemes/ToricIdealSheaves/constructors.jl")

include("Schemes/ToricDivisors/constructors.jl")
include("Schemes/ToricDivisors/attributes.jl")

include("Schemes/NormalToricVarieties/attributes.jl")

include("Schemes/ToricBlowups/types.jl")
include("Schemes/ToricBlowups/constructors.jl")
include("Schemes/ToricBlowups/attributes.jl")
include("Schemes/ToricBlowups/methods.jl")

include("ExteriorAlgebra/ExteriorAlgebra.jl")

include("Schemes/DerivedPushforward.jl")

