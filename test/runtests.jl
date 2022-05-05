using Oscar
using Test

include("printing.jl")

include("PolyhedralGeometry/runtests.jl")
include("Combinatorics/runtests.jl")

include("GAP/runtests.jl")

include("Rings/integer-test.jl")
include("Rings/rational-test.jl")
include("Rings/mpoly-test.jl")
include("Rings/orderings-test.jl")
include("Rings/affine-alghom-test.jl")
include("Rings/mpoly-graded-test.jl")
include("Rings/mpoly-local-test.jl")
include("Rings/mpoly-localizations.jl")
include("Rings/mpolyquo-localizations.jl")
include("Rings/integer-localizations.jl")
include("Rings/nmod-localizations.jl")
include("Rings/mpoly-nested-test.jl")
include("Rings/MPolyQuo_test.jl")
include("Rings/groebner-test.jl")
include("Rings/msolve-test.jl")
include("Rings/FractionalIdeal-test.jl")
include("Rings/mpoly_affine_algebras_test.jl")
include("Rings/slpolys-test.jl")
include("NumberTheory/nmbthy-test.jl")
include("Groups/runtests.jl")
include("Rings/NumberField.jl")
include("Rings/FunctionField-test.jl")
include("Rings/AbelianClosure.jl")

include("Rings/MPolyAnyMap/MPolyRing.jl")
include("Rings/MPolyAnyMap/MPolyQuo.jl")
include("Rings/MPolyAnyMap/AffineAlgebras.jl")

if Oscar.is_dev
  include("Examples/GITFans-test.jl")
end

include("Rings/binomial-ideals-test.jl")

include("Experimental/PlaneCurve-test.jl")
include("Experimental/galois-test.jl")
include("Experimental/gmodule-test.jl")
include("Experimental/ModStdQt-test.jl")
include("Experimental/ModStdNF-test.jl")
include("Experimental/Permutations-test.jl")

include("Modules/UngradedModules.jl")
include("Modules/ModulesGraded.jl")

include("InvariantTheory/runtests.jl")

include("ToricVarieties/runtests.jl")

include("Schemes/AffineSchemes.jl")
include("Schemes/SpecOpen.jl")
include("Schemes/Glueing.jl")
include("Schemes/ProjectiveSchemes.jl")

include("TropicalGeometry/runtests.jl")
