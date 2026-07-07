struct LadderStep
  is_up_step::Bool
  A::PermGroup
  Aprev::PermGroup
  T::Vector
  Tmap::Map
  D::Vector
  St::Vector
  I::Vector
  m::Map

  function LadderStep(A::PermGroup, Aprev::PermGroup)
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
    # construct T / Tmap or do this only when needed?
    return rec
  end

  function LadderStep(A::PermGroup, Aprev::PermGroup, T::Vector)
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
    # construct a tmap
    return rec
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

struct SubgroupLadder <: AbstractVector{LadderStep}
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
