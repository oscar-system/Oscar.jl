module SymInt
using Oscar, Markdown

macro req(args...)
  @assert length(args) == 2
  quote
    if !($(esc(args[1])))
      throw(ArgumentError($(esc(args[2]))))
    end
  end
end

include("Types.jl")
include("Elevators.jl")
include("Representations.jl")
include("SymmetricGrassmannians.jl")
include("HomogeneousPolynomialsActions.jl")
include("K3Models.jl")

end
