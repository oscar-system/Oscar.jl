using Oscar
using Test

import Oscar.Nemo.AbstractAlgebra
include(joinpath(pathof(AbstractAlgebra), "..", "..", "test", "Rings-conformance-tests.jl"))

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
  include("Experimental/GITFans-test.jl")
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
include("Schemes/FunctionFields.jl")
include("Modules/ProjectiveModules.jl")
include("Schemes/singular_locus.jl")
include("Schemes/SpaceGerms.jl")
include("Schemes/SpecialTypes.jl")

include("TropicalGeometry/runtests.jl")
include("Serialization/runtests.jl")

include("StraightLinePrograms/runtests.jl")
