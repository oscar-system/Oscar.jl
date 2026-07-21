mutable struct TransversalChain
  C::Vector{Tuple{PermGroup, Vector{PermGroupElem}}}
end

data(T::TransversalChain) = T.C

function Base.size(T::TransversalChain)
  return size(data(T))
end

function Base.length(T::TransversalChain)
  return length(data(T))
end

Base.firstindex(T::TransversalChain) = 1
Base.lastindex(T::TransversalChain) = length(T)

function Base.getindex(T::TransversalChain, i::IntegerUnion)
  return getindex(data(T), Int(i))
end


function TransversalChain(U::PermGroup)
  U0 = U
  C = Tuple{PermGroup, Vector{PermGroupElem}}[]
  while length(U0) != 1
    s = first(moved_points(U0))
    U = stabilizer(U0, s)[1]
    T = right_transversal(U0, U)
    push!(C, (U0, T))
    U0 = U
  end
  push!(C, (U0, [one(U0)]))

  return TransversalChain(C)
end

get_transversal_chain(U::PermGroup) = TransversalChain(U)




mutable struct LadderStep
  is_up_step::Bool

  Aprev::PermGroup
  A::PermGroup

  T::Union{Vector{PermGroupElem}, Oscar.SubgroupTransversal{PermGroup, PermGroup, PermGroupElem}}
  Tmap::Map{PermGroup, PermGroup}

  DSt::Dict{PermGroupElem, TransversalChain}


  Im::Dict{PermGroupElem, Tuple{PermGroupElem, PermGroupElem}}

  F

  function LadderStep(Aprev::PermGroup, A::PermGroup)
    if is_subset(Aprev, A)
      u = true
    elseif is_subset(A, Aprev)
      u = false
    else
      # throw an error
    end
    return new(u, Aprev, A)
  end

  function LadderStep(A::PermGroup, Aprev::PermGroup, T::Vector, Tmap::Map)
    if is_subset(Aprev, A)
      u = true
    elseif is_subset(A, Aprev)
      u = false
    else
      # throw an error
    end
    return new(u, Aprev, A, T, Tmap)
  end
end

function Base.getproperty(S::LadderStep, s::Symbol)
  if s === :F
    return getfield(S, :F)::SubArray{LadderStep, 1, SubgroupLadder, Tuple{UnitRange{Int64}}, false}
  else
    return getfield(S, s)
  end
end



function _get_previous_steps(S::LadderStep)
  return S.F
end

# function _set_previous_steps!(S::LadderStep, )


# SubgroupLadder stores Vector{LadderStep}
mutable struct SubgroupLadder <: AbstractVector{LadderStep}
  S::Vector{LadderStep}

  function SubgroupLadder(S::Vector{LadderStep})
    L = new(S)
    for (i, s) in enumerate(L)
      s.F = view(L, 1:i)
    end
    return L
  end
end

function SubgroupLadder(S::Vector{PermGroup})
  LS = LadderStep[]
  isempty(S) && return new(LS)
  push!(LS, LadderStep(first(S), first(S)))
  for i in 2:length(S)
    push!(L, LadderStep(S[i-1], S[i]))
  end

  return SubgroupLadder(LS)
end

data(L::SubgroupLadder) = L.S

################################################################################
#
#  Array-like functionality
#
################################################################################

function Base.size(L::SubgroupLadder)
  return size(data(L))
end

function Base.length(L::SubgroupLadder)
  return length(data(L))
end

Base.firstindex(L::SubgroupLadder) = 1
Base.lastindex(L::SubgroupLadder) = length(L)

function Base.getindex(L::SubgroupLadder, i::IntegerUnion)
  return getindex(data(L), Int(i))
end

# I don't think we want to support this
# function Base.setindex!(L::SubgroupLadder, s::LadderStep, i::IntegerUnion)
#   return setindex!(data(L), s, Int(i))
# end

function Base.copy(L::SubgroupLadder)
  return SubgroupLadder(copy(data(L)))
end

# function Base.push!(L::SubgroupLadder, s::LadderStep)
#   # require or assert that this is actually a step
#   @req last(L).A == s.Aprev "Incompatible LadderStep"
#   push!(L.S, s)
#   return L
# end
