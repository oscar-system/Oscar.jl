###############################
### 1: Types for toric schemes
###############################

@attributes mutable struct ToricCoveredScheme{BRT} <: AbsCoveredScheme{BRT}
  ntv::NormalToricVariety
  ToricCoveredScheme(ntv::NormalToricVariety) = new{typeof(coefficient_ring(ntv))}(ntv)
end

@attributes mutable struct ToricSpec{BRT, RT} <: AbsSpec{BRT, RT}
  antv::AffineNormalToricVariety
  var_names::Vector{Symbol}
  ToricSpec(antv::AffineNormalToricVariety; var_names::Vector{Symbol} = Vector{Symbol}()) = new{typeof(coefficient_ring(antv)), MPolyQuoRing}(antv, var_names)
end



###############################
### 2: Display
###############################

function Base.show(io::IO, X::ToricSpec)
  print(io, "Spec of an affine toric variety with cone spanned by $(rays(polyhedral_fan(X)))")
end

function Base.show(io::IO, X::ToricCoveredScheme)
  print(io, "Scheme of a toric variety with fan spanned by $(rays(polyhedral_fan(X)))")
end
