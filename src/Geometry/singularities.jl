const CRing = Union{MPolyRing, MPolyQuo{<:Oscar.MPolyElem}, MPolyRing_dec, MPolyQuo{<:Oscar.MPolyElem_dec}}
const CRingElem = Union{MPolyElem, MPolyQuoElem{<:Oscar.MPolyElem}, MPolyElem_dec, MPolyQuoElem{<:Oscar.MPolyElem_dec}}

struct Germ{T <: CRingElem} 
  repres::T
  milnor_number::Int
  #milnor_algebra
end

function Base.getproperty(g::Germ, s::Symbol)
  if s == :milnor_number
    #setfield!(g,s, 10)
  elseif s == :milnor_algebra
    #setfield!(g,:milnor_algebra, _)
  else
    return getfield(g, s)
  end
end