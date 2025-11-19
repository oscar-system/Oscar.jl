################################################################################
# GraphDict
struct GraphDict{T}
  d::Dict{<:Union{Int, Edge}, T}
end

Base.getindex(D::GraphDict, i::Int) = D.d[i]
Base.getindex(D::GraphDict, e::Edge) = D.d[e]

function Base.getindex(D::GraphDict, i::Int, j::Int)
  return D.d[Edge(i, j)]
end

################################################################################
# GenDict
struct GenDict{S} 
  d::Dict{S, <: MPolyRingElem}
end

Base.getindex(D::GenDict{S}, arg::S) where S = D.d[arg]

function Base.getindex(D::GenDict{Tuple{Set{Int}, Tuple}}, s::Vector{Int}, t::Tuple)
  return D.d[Set{Int}(s), t]
end

function Base.getindex(D::GenDict{Tuple{Int, Int, Set{Int}}}, i::Int, j::Int, s::Vector{Int})
  return D.d[i, j, Set{Int}(s)]
end

################################################################################
# GenTransDict
struct GraphTransDict{T}
  d::Dict{Tuple{VarName, Edge}, T}
end

Base.getindex(D::GraphTransDict, arg::Tuple{VarName, Edge}) = D.[arg]

function Base.getindex(D::GraphTransDict, s::VarName, i::Int, j::Int)
  return D.d[(s, Edge(i, j))]
end

################################################################################
# FinGenAbGroupElemDict
struct FinAbGroupElemDict{T}
  d::Dict{FinGenAbGroupElem, T}
end

Base.getindex(D::FinAbGroupElemDict, arg::FinGenAbGroupElem) = D.d[arg]

function Base.getindex(D::FinAbGroupElemDict, arg::Tuple{Vararg{Int}})
  isempty(D.d) && throw(KeyError(arg))
  G = parent(first(D.d).first)
  key = G(collect(arg))
  return d[key]
end

################################################################################
# dict functions
const SpecialDictUnion = Union{GraphDict, GenDict, GraphTransDict, FinAbGroupElemDict}

keys(D::T) where T <: SpecialDictUnion = keys(D.d)
values(D::T) where T <: SpecialDictUnion = values(D.d)
isempty(D::T) where T <: SpecialDictUnion = isempty(D.d)
show(D::T) where T <: SpecialDictUnion = show(D.d)


