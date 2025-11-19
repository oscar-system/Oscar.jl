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

Base.getindex(D::GraphTransDict, arg::Tuple{VarName, Edge}) = D.d[arg]

Base.getindex(D::GraphTransDict, s::VarName, e::Edge) = D.d[(s, e)]

Base.getindex(D::GraphTransDict, s::VarName, i::Int, j::Int) = D.d[(s, Edge(i, j))]

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

Base.keys(D::T) where T <: SpecialDictUnion = keys(D.d)
Base.values(D::T) where T <: SpecialDictUnion = values(D.d)
Base.isempty(D::T) where T <: SpecialDictUnion = isempty(D.d)
Base.length(D::T) where T <: SpecialDictUnion = length(D.d)
Base.iterate(D::T) where T <: SpecialDictUnion = iterate(D.d)
Base.iterate(D::T, i::Int) where T <: SpecialDictUnion = iterate(D.d, i)

function Base.show(io::IO,  m::MIME"text/plain", D::T) where T <: SpecialDictUnion
  io = pretty(io)
  println(io, "$T with underlying Dict")
  print(io, Indent())
  show(io,  m, D.d)
  print(io, Dedent())
end


