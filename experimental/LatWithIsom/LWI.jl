module LWI
using Oscar, Markdown
import Hecke

macro req(args...)
  @assert length(args) == 2
  quote
    if !($(esc(args[1])))
      throw(ArgumentError($(esc(args[2]))))
    end
  end
end

include("Types.jl")
include("HMM.jl")
include("LatWithIsom.jl")
include("Enumeration.jl")
include("Embeddings.jl")

end
