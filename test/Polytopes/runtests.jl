const pm = Polymake

using Oscar

# This enables a debug flag that will produce some more output
# during polymake's wrapper compilation (at run-time).
# We want to keep this active for CI to make it easier to debug future
# compilation errors, especially when they only appear during CI.
if (haskey(ENV, "GITHUB_ACTIONS"))
    Polymake.shell_execute(raw"$Verbose::cpp=3;")
end

include("types.jl")
include("iterators.jl")
include("Cone.jl")
include("Group.jl")
include("Polyhedron.jl")
include("PolyhedralFan.jl")
include("SubdivisionOfPoints.jl")
include("extended.jl")
include("Graph.jl")
