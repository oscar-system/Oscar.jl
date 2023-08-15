# We read experimental and filter out all packages that follow our desired
# scheme. Remember those packages to avoid doing this all over again for docs
# and test.
# We don't want to interfere with existing stuff in experimental though.
const expdir = joinpath(@__DIR__, "../experimental")
const oldexppkgs = [
  "ExteriorAlgebra",
  "GModule",
  "Matrix",
  "ModStd",
  "Rings",
  "Schemes"
]
# DEVELOPER OPTION:
# The following lines ensure that ToricSchemes is loaded before FTheoryTools.
# DO NOT USE THIS UNLESS YOU KNOW THE CONSEQUENCES.
# For more background, see https://github.com/oscar-system/Oscar.jl/issues/2300.
const orderedpkgs = [
  "ToricSchemes",
  "FTheoryTools",
  "JuLie",
  "IntersectionTheory",
  "OrthogonalDiscriminants",  # needs code from JuLie
]
exppkgs = filter(x->isdir(joinpath(expdir, x)) && !(x in oldexppkgs) && !(x in orderedpkgs), readdir(expdir))
append!(exppkgs, orderedpkgs)

# Error if something is incomplete in experimental
for pkg in exppkgs
  if !isfile(joinpath(expdir, pkg, "src", "$pkg.jl"))
    error("experimental/$pkg is incomplete: $pkg/src/$pkg.jl missing. See the documentation at https://docs.oscar-system.org/dev/Experimental/intro/ for details.")
  end
  if !isfile(joinpath(expdir, pkg, "test", "runtests.jl"))
    error("experimental/$pkg is incomplete: $pkg/test/runtests.jl missing. See the documentation at https://docs.oscar-system.org/dev/Experimental/intro/ for details.")
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

include("Schemes/Types.jl")
include("Schemes/SpecialTypes.jl")
include("Schemes/CoveredScheme.jl")
include("Schemes/FunctionFields.jl")
include("Schemes/ProjectiveModules.jl")
include("Schemes/SpaceGerms.jl")
include("Schemes/Sheaves.jl")
include("Schemes/IdealSheaves.jl")
include("Schemes/AlgebraicCycles.jl")
include("Schemes/WeilDivisor.jl")
include("Schemes/CoveredProjectiveSchemes.jl")

include("Matrix/matrix.jl")
include("Schemes/SimplifiedSpec.jl")
include("Schemes/CoherentSheaves.jl")
include("Schemes/LazyGlueing.jl")
include("Schemes/CartierDivisor.jl")
include("Schemes/Auxiliary.jl")
include("Schemes/BlowupMorphism.jl")
include("Schemes/duValSing.jl")
include("Schemes/elliptic_surface.jl")
include("Schemes/RationalMap.jl")

include("ExteriorAlgebra/ExteriorAlgebra.jl")
