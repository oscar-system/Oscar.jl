################################################################################
# important dict type, leaving here for now
const GraphDict{T} = Dict{<:Union{Int, Edge}, T}
const GenDict{S, T} = Dict{S, T} where {S, T <: MPolyRingElem}

function Base.getindex(D::GraphDict, i::Int, j::Int)
  return D[Edge(i, j)]
end

function Base.getindex(D::GenDict{Tuple{Set{Int}, T}}, s::Vector{Int}, t::T) where T <: Tuple{Vararg{Int}}
  return D[Set{Int}(s), t]
end

function Base.getindex(D::GenDict{Tuple{Int, Int, Set{Int}}}, i::Int, j::Int, s::Vector{Int})
  return D[i, j, Set{Int}(s)]
end

const GraphTransDict{T} = Dict{Tuple{VarName, Edge}, T} where T <: MPolyRingElem

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
