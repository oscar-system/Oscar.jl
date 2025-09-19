const FinAbGroupElemDict{T} = Dict{FinGenAbGroupElem, T} where T

function Base.getindex(d::FinAbGroupElemDict, arg::Tuple{Vararg{Int}})
  isempty(d) && throw(KeyError(arg))
  G = parent(first(d).first)
  key = G(collect(arg))
  return d[key]
end
