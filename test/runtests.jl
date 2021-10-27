using Oscar
using Test

include("Polytopes/runtests.jl")

include("GAP/runtests.jl")

include("Rings/integer-test.jl")
include("Rings/rational-test.jl")
include("Rings/mpoly-test.jl")
include("Rings/affine-alghom-test.jl")
include("Rings/mpoly-graded-test.jl")
include("Rings/mpoly-local-test.jl")
include("Rings/mpoly-localizations.jl")
include("Rings/integer-localizations.jl")
include("Rings/nmod-localizations.jl")
include("Rings/mpoly-nested-test.jl")
include("Rings/MPolyQuo_test.jl")
include("Rings/msolve-test.jl")
include("Rings/FractionalIdeal-test.jl")
include("Rings/mpoly_affine_algebras_test.jl")
include("Rings/slpolys-test.jl")
include("NumberTheory/nmbthy-test.jl")
include("Groups/runtests.jl")
include("Rings/NumberField.jl")
include("Rings/FunctionField-test.jl")
include("Rings/AbelianClosure.jl")

if Oscar.is_dev
  include("Examples/PlaneCurve-test.jl")
end

include("Rings/binomial-ideals-test.jl")

include("Examples/galois-test.jl")
include("Examples/ModStdQt-test.jl")
include("Examples/ModStdNF-test.jl")

include("Modules/UngradedModules.jl")

include("InvariantTheory/invariant_rings-test.jl")
include("InvariantTheory/primary_invariants-test.jl")

include("ToricVarieties/runtests.jl")
