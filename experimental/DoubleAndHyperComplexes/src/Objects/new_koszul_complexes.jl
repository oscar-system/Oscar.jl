########################################################################
# Homogeneous Koszul complexes
#
# Given a graded ring R and homogeneous elements a₁,…,aₙ ∈ R of degrees
# dᵢ = deg aᵢ one wishes to build the Koszul complex 
#
# R¹ ← R[-d₁] ⊕ … ⊕ R[-dᵢ] ← R[-d₁-d₂] ⊕ R[-d₁-d₃] ⊕ … 
#
# which is homogeneous with all differentials of relative degree zero. 
# 
# The construction via exterior multiplication with ∑ᵢaᵢ⋅eᵢ ∈ Rⁿ does 
# not provide this. So we introduce an extra method here. 
########################################################################

### Production of the chains
struct HomogKoszulComplexChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  S::Ring
  seq::Vector{<:RingElem}
  degs::Vector{FinGenAbGroupElem}

  function HomogKoszulComplexChainFactory(S::Ring, seq::Vector{<:RingElem}; check::Bool=true)
    @check all(is_homogeneous(a) for a in seq) "elements must be homogeneous"
    return new{FreeMod{elem_type(S)}}(S, seq, [degree(a; check) for a in seq])
  end
end

function (fac::HomogKoszulComplexChainFactory)(self::AbsHyperComplex, I::Tuple)
  i = first(I)
  G = grading_group(fac.S)
  is_zero(i) && return graded_free_module(fac.S, [zero(G)])
  n = length(fac.seq)
  degs = elem_type(G)[-sum(fac.degs[indices(omi)]; init=zero(G)) for omi in OrderedMultiIndexSet(i, n)]
  result = graded_free_module(fac.S, degs)
end

function can_compute(fac::HomogKoszulComplexChainFactory, self::AbsHyperComplex, i::Tuple)
  return 0 <= first(i) <= length(fac.seq)
end

### Production of the morphisms 
struct HomogKoszulComplexMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  function HomogKoszulComplexMapFactory(S::Ring)
    return new{FreeModuleHom}()
  end
end

function (fac::HomogKoszulComplexMapFactory)(self::AbsHyperComplex, p::Int, I::Tuple)
  i = first(I)
  cfac = chain_factory(self)
  n = length(cfac.seq)
  dom = self[I]
  cod = self[i+1]
  img_gens = elem_type(cod)[]
  inds = [OrderedMultiIndex([k], n) for k in 1:n]
  for omi in OrderedMultiIndexSet(i, n)
    img = zero(cod)
    for (m, v) in zip(inds, cfac.seq)
      sign, mult_ind = _wedge(omi, m)
      is_zero(sign) && continue
      res_ind = linear_index(mult_ind)
      img += sign*v*cod[res_ind]
    end
    push!(img_gens, img)
  end
  return hom(dom, cod, img_gens)
end

function can_compute(fac::HomogKoszulComplexMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  isone(p) || return false
  return 0 <= first(i) < length(chain_factory(self).seq)
end

### The concrete struct
@attributes mutable struct HomogKoszulComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}

  function HomogKoszulComplex(R::Ring, seq::Vector{<:RingElem}; check::Bool=true)
    @assert is_graded(R)
    @assert all(parent(a) === R for a in seq)
    chain_fac = HomogKoszulComplexChainFactory(R, seq; check)
    map_fac = HomogKoszulComplexMapFactory(R)

    # Assuming d is the dimension of the new complex
    internal_complex = HyperComplex(1, chain_fac, map_fac, [:cochain],
                                    upper_bounds=Union{Int, Nothing}[length(seq)],
                                    lower_bounds=Union{Int, Nothing}[0]
                                   )
    # Assuming that ChainType and MorphismType are provided by the input
    return new{FreeMod{elem_type(R)}, FreeModuleHom}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::HomogKoszulComplex) = c.internal_complex

