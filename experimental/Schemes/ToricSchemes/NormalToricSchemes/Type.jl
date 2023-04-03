@attributes mutable struct ToricCoveredScheme{BRT} <: AbsCoveredScheme{BRT}
  ntv::NormalToricVariety
  X::CoveredScheme

  function ToricCoveredScheme(ntv::NormalToricVariety)
    return new{typeof(QQ)}(ntv)
  end
end

function Base.show(io::IO, X::ToricCoveredScheme)
  print(io, "Scheme of a toric variety with fan spanned by $(rays(fan(X)))")
end

