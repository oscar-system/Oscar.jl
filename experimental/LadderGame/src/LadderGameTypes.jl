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


struct SubgroupLadder
  L::Vector{LadderStep}

  function SubgroupLadder(S::Vector{PermGroup})
    rec = new()

    rec.L = [ LadderStep(S[i+1], S[i])]
  end
end
