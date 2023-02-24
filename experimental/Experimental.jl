for pkg in Oscar.exppkgs
  include("$pkg/src/$pkg.jl")
end

include("GaloisGrp.jl")
include("Rings.jl")
include("ModStd.jl")
include("PlaneCurve.jl")
include("IntersectionTheory.jl")
