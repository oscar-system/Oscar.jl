for pkg in Oscar.exppkgs
  include("$pkg/test/runtests.jl")
end

# legacy files
include("GModule/test/runtests.jl")
include("ModStd/test/runtests.jl")
include("MatrixGroups/test/runtests.jl")
include("ExteriorAlgebra/test/runtests.jl")
