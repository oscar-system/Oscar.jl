using Oscar
using Test

# Used in both Rings/slpolys-test.jl and StraightLinePrograms/runtests.jl
const SLP = Oscar.StraightLinePrograms
include("printing.jl")

include("PolyhedralGeometry/runtests.jl")
include("Combinatorics/runtests.jl")

include("GAP/runtests.jl")
include("Groups/runtests.jl")

include("Rings/runtests.jl")

include("NumberTheory/nmbthy-test.jl")

if Oscar.is_dev
  include("Examples/GITFans-test.jl")
end

include("Experimental/PlaneCurve-test.jl")
include("Experimental/galois-test.jl")
include("Experimental/gmodule-test.jl")
include("Experimental/ModStdQt-test.jl")
include("Experimental/ModStdNF-test.jl")
include("Experimental/MPolyRingSparse-test.jl")
include("Experimental/MatrixGroups-test.jl")

include("Modules/UngradedModules.jl")
include("Modules/ModulesGraded.jl")
include("Modules/module-localizations.jl")

include("InvariantTheory/runtests.jl")

include("ToricVarieties/runtests.jl")

include("Schemes/AffineSchemes.jl")
include("Schemes/SpecOpen.jl")
include("Schemes/Glueing.jl")
include("Schemes/ProjectiveSchemes.jl")
include("Schemes/CoveredScheme.jl")

include("TropicalGeometry/runtests.jl")
include("Serialization/runtests.jl")

include("StraightLinePrograms/runtests.jl")
