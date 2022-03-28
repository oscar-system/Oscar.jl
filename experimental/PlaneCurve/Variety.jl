module VarietyModule

using Oscar

export Variety

mutable struct Variety
  I::Oscar.MPolyIdeal

  function Variety(I::Oscar.MPolyIdeal{<:MPolyElem{<:FieldElem}})
    return new(I)
  end
end

function Base.show(io::IO, V::Variety)
  print(io, "Variety defined by the ", V.I)
end


end
using .VarietyModule
