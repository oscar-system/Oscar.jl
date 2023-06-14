import .AbstractAlgebra.Generic: LaurentMPolyWrapRing, LaurentMPolyWrap
import .AbstractAlgebra: LaurentMPolyRing, LaurentMPolyRingElem

@attributes mutable struct LaurentMPolyAnyMap{D, C} <: Map{D, C, Hecke.Hecke.Map, LaurentMPolyAnyMap}
  R::D
  S::C
  image_of_gens
  data # map from _polyringquo(R) -> C
  
  function LaurentMPolyAnyMap(R::D, S::C, image_of_gens) where {D, C}
    @assert all(x -> parent(x) === S, image_of_gens)
    return new{D, C}(R, S, image_of_gens)
  end
end

mutable struct LaurentMPolyIdeal{T}
  R
  gens::Vector{T}
  data # ideal of of _polyringquo(R)

  function LaurentMPolyIdeal(R::LaurentMPolyRing, gens::Vector)
    @assert all(x -> parent(x) === R, gens)
    return new{eltype(gens)}(R, gens)
  end
end

mutable struct _LaurentMPolyBackend{D, C, M}
  R::D
  Q::C
  inv::M
  _gens_cache

  function _LaurentMPolyBackend(R::D, Q::C) where {D, C}
    _inv = hom(Q, R, vcat(gens(R), inv.(gens(R))))
    return new{D, C, typeof(_inv)}(R, Q, _inv)
  end
end


