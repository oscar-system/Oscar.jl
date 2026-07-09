mutable struct LadderStep
  is_up_step::Bool

  Aprev::PermGroup
  A::PermGroup


  T::Union{Vector{PermGroupElem}, Oscar.SubgroupTransversal}
  Tmap::Map

  # Should these also be a Dict??????
  D::Vector
  St::Vector

  # TODO I/m should be implemented as a Dict
  I::Vector
  m::Vector

  F

  function LadderStep(Aprev::PermGroup, A::PermGroup)
    S = new()
    if is_subset(Aprev, A)
      S.is_up_step = true
    elseif is_subset(Aprev, A)
      S.is_up_step = false
    else
      # throw an error
    end
    S.A = A
    S.Aprev = Aprev
    # construct T / Tmap or do this only when needed?
    return S
  end

  function LadderStep(Aprev::PermGroup, A::PermGroup, T::Vector)
    S = new()
    if is_subset(Aprev, A)
      S.is_up_step = true
    elseif is_subset(Aprev, A)
      S.is_up_step = false
    else
      # throw an error
    end
    S.A = A
    S.Aprev = Aprev
    S.T = T
    # construct a tmap
    return S
  end

  function LadderStep(A::PermGroup, Aprev::PermGroup, T::Vector, Tmap::Map)
    rec = new()
    if is_subset(Aprev, A)
      rec.is_up_step = true
    elseif is_subset(Aprev, A)
      rec.is_up_step = false
    else
      # throw an error
    end
    rec.A = A
    rec.Aprev = Aprev
    rec.T = T
    rec.Tmap = Tmap
    return rec
  end

end


# SubgroupLadder stores Vector{LadderStep}
# with some special tools for adding steps.

mutable struct SubgroupLadder <: AbstractVector{LadderStep}
  S::Vector{LadderStep}

  function SubgroupLadder(S::Vector{PermGroup})
    rec = new()

    rec.S = [ LadderStep(S[i+1], S[i])]
  end
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

function Base.push!(L::SubgroupLadder, s::LadderStep)
  # require or assert that this is actually a step
  @req last(L).A == s.Aprev "Incompatible LadderStep"
  push!(L.S, s)
  return L
end
