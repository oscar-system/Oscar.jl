using Oscar
using Test

include("GAP/runtests.jl")

include("Rings/integer-test.jl")
include("Rings/rational-test.jl")
include("Rings/mpoly-test.jl")
include("Rings/mpoly-graded-test.jl")
include("Rings/mpoly-local-test.jl")
include("Rings/MPolyQuo.jl")
include("Rings/slpolys-test.jl")
include("Polymake/nmbthy-test.jl")
include("Groups/runtests.jl")

if Oscar.is_dev
  include("Examples/PlaneCurve-test.jl")
  include("Examples/galois-test.jl")
  include("Examples/ModStdQt-test.jl")
  include("Examples/ModStdNF-test.jl")
end
