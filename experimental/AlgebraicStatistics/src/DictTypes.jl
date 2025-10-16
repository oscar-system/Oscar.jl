################################################################################
# important dict type, leaving here for now
const GraphGenDict = Dict{Union{Int, Edge}, T} where T <: MPolyRingElem
const GenDict = Dict{S, T} where {S, T <: MPolyRingElem}

function Base.getindex(D::GraphGenDict, i::Int, j::Int)
  return D[Edge(i, j)]
end

const GraphTransDict = Dict{Tuple{VarName, Edge}, T} where T <: MPolyRingElem

function Base.getindex(D::GraphTransDict, s::VarName, i::Int, j::Int)
  return D[(s, Edge(i, j))]
end

const FinAbGroupElemDict{T} = Dict{FinGenAbGroupElem, T} where T

function Base.getindex(d::FinAbGroupElemDict, arg::Tuple{Vararg{Int}})
  isempty(d) && throw(KeyError(arg))
  G = parent(first(d).first)
  key = G(collect(arg))
  return d[key]
end
