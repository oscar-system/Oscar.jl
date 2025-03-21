### Production of the chains
struct KoszulChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  v::FreeModElem

  function KoszulChainFactory(v::FreeModElem)
    return new{typeof(parent(v))}(v)
  end
end

function (fac::KoszulChainFactory)(c::AbsHyperComplex, i::Tuple)
  r = rank(parent(fac.v))
  if first(i) == -1 || first(i) == r + 1
    M = c[0]
    R = base_ring(M)
    return (is_graded(M) ? graded_free_module(R, elem_type(grading_group(R))[]) : FreeMod(R, 0))
  end
  return exterior_power(parent(fac.v), r-first(i))[1]
end

function can_compute(fac::KoszulChainFactory, c::AbsHyperComplex, i::Tuple)
  r = rank(parent(fac.v))
  return -1 <= first(i) <= r + 1
end

### Production of the morphisms 
struct KoszulMorphismFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  v::FreeModElem

  function KoszulMorphismFactory(v::FreeModElem)
    return new{morphism_type(parent(v), parent(v))}(v)
  end
end

function (fac::KoszulMorphismFactory)(c::AbsHyperComplex, p::Int, i::Tuple)
  r = rank(parent(fac.v))
  dom = c[i]
  cod = c[first(i) - 1]
  if first(i) == 0 || first(i) == r + 1
    return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)]; check=false)
  end
  return wedge_multiplication_map(dom, cod, fac.v)
end

function can_compute(fac::KoszulMorphismFactory, c::AbsHyperComplex, p::Int, i::Tuple)
  r = rank(parent(fac.v))
  return 0 <= first(i) <= r + 1
end

### The concrete struct
@attributes mutable struct KoszulComplex{ChainType, MorphismType} <: AbsSimpleComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}
  v::FreeModElem

  function KoszulComplex(v::FreeModElem)
    chain_fac = KoszulChainFactory(v)
    map_fac = KoszulMorphismFactory(v)

    # Assuming d is the dimension of the new complex
    internal_complex = HyperComplex(1, chain_fac, map_fac, [:chain], lower_bounds=[0], upper_bounds=[rank(parent(v))])
    # Assuming that ChainType and MorphismType are provided by the input
    return new{typeof(parent(v)), morphism_type(parent(v), parent(v))}(internal_complex, v)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::KoszulComplex) = c.internal_complex

### User facing constructor
function koszul_complex(::Type{T}, v::FreeModElem) where {T<:AbsSimpleComplex}
  return KoszulComplex(v)
end

function koszul_complex(::Type{T}, v::Vector{<:RingElem}) where {T<:AbsSimpleComplex}
  isempty(v) && error("can not build from empty list")
  R = parent(first(v))
  @assert all(u->parent(u)===R, v) "elements must belong to the same ring"
  r = length(v)
  F = (is_graded(R) ? graded_free_module(R, degree.(v)) : FreeMod(R, r))
  w = sum(u*e for (u, e) in zip(v, gens(F)); init = zero(F))
  return koszul_complex(T, w)
end

function koszul_complex(::Type{T}, v::RingElem...) where {T<:AbsSimpleComplex}
  return koszul_complex(T, collect(v))
end

function defining_element(K::KoszulComplex)
  return K.v
end
#=
function truncated_cech_complex(M::ModuleFP, x::Vector{T}, d::Int) where {T<:RingElem}
  R = base_ring(M)
  @assert all(u->parent(u)===R, x) "elements do not belong to the correct ring"
  @assert d>=0 "pole order must be non-negative"

  r = length(x)
  F = (is_graded(M) ? graded_free_module(R, degree.(x)) : FreeMod(R, r))
  v = sum(u^d*e for (u, e) in zip(x, gens(F)); init = zero(F))

  kosz = koszul_complex(AbsSimpleComplex, v)
  return hom(kosz, M)
end
=#
