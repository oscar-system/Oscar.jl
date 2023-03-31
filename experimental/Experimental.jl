for pkg in Oscar.exppkgs
  include("$pkg/src/$pkg.jl")
end

include("GaloisGrp.jl")
include("Rings.jl")
include("ModStd.jl")
include("GITFans.jl")
include("GModule.jl")
include("MPolyRingSparse.jl")
include("SymmetricIntersections.jl")
include("LinearQuotients.jl")

include("JuLie.jl")

include("Schemes/Types.jl")
include("Schemes/SpecialTypes.jl")
include("Schemes/CoveredScheme.jl")
include("Schemes/FunctionFields.jl")
include("Schemes/ProjectiveModules.jl")
include("Schemes/SpaceGerms.jl")
include("Schemes/singular_locus.jl")
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
include("Schemes/ToricSchemes/include.jl")

include("ExteriorAlgebra/ExteriorAlgebra.jl")

