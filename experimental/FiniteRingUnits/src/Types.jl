struct EffectivePresentation{S}
  A::S
  G#=::FPGroup=#
  forward
  backward

  # G an FPGroup (most of the time)
  function EffectivePresentation(A::S, G, forward, backward) where {S}
    z = new{S}(A, G, forward, backward)
     return z
  end
end

struct RingMultMap{S, T} <: Map{S, T, Any, RingMultMap}
  R::S
  A::T
  f # map from R to A
  g # map from A to R

  function RingMultMap(R::S, A::T, f, g) where {S, T}
    z = new{S, T}(R, A, f, g)
  end
end
