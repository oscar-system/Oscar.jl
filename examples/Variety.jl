module VarietyModule

using Oscar

export Variety

mutable struct Variety
  I::Oscar.MPolyIdeal

  function Variety(I::Oscar.MPolyIdeal{<:MPolyElem{<:FieldElem}})
    r = new()
    r.I = I
    return r
  end
end

function Base.show(io::IO, V::Variety)
  println(io, "Variety defined by ", V.I)
end


end
using .VarietyModule
