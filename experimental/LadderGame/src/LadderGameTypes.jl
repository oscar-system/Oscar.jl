struct LadderData
  A::PermGroup
  Aprev::PermGroup
  T::Vector
  Tmap::Map
  D::Vector
  St::Vector
  I::Vector
  m::Map

  function LadderData(A::PermGroup, Aprev::PermGroup)
    rec = new()
    rec.A = A
    rec.Aprev = Aprev
    return rec
  end
end
