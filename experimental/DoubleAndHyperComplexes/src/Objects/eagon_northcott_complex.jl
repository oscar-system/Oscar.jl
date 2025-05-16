### Production of the chains
struct ENCChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  A::MatrixElem
  exps::Vector{Vector{Vector{Int}}} # a hashtable for indexing
  F::ChainType # the original free module for the Koszul complex
  sym_powers::Dict{Int, ChainType}
  ext_powers::Dict{Int, ChainType}

  function ENCChainFactory(
      A::MatrixElem{T},
      F::FreeMod
    ) where {T<:RingElem}
    n = ncols(A)
    k = nrows(A)
    exps = Vector{Vector{Vector{Int}}}()
    for i in 0:(n-k+1)
      push!(exps, [a.c for a in WeakCompositions(i, k)])
    end
    ChainType = FreeMod{T}
    return new{ChainType}(A, exps, F, Dict{Int, ChainType}(), Dict{Int, ChainType}())
  end
end

function (fac::ENCChainFactory)(self::AbsHyperComplex, I::Tuple)
  i = first(I)
  R = base_ring(fac.A)
  if is_zero(i)
    return is_graded(R) ? graded_free_module(R, [zero(grading_group(R))]) : FreeMod(R, 1)
  end
  Sd = is_graded(R) ? graded_free_module(R, [zero(grading_group(R)) for _ in 1:length(fac.exps[i])]) : FreeMod(R, length(fac.exps[i]))
  n = ncols(fac.A)
  k = nrows(fac.A)
  wedge_power, _ = exterior_power(fac.F, k + i - 1)
  fac.ext_powers[i] = wedge_power
  fac.sym_powers[i] = Sd
  return tensor_product(Sd, wedge_power)
end

function can_compute(fac::ENCChainFactory, self::AbsHyperComplex, I::Tuple)
  is_one(length(I)) || return false
  i = first(I)
  i >= 0 || return false
  return i <= ncols(fac.A) - nrows(fac.A) + 1
end

### Production of the morphisms
struct ENCMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  function ENCMapFactory()
    return new{FreeModuleHom}()
  end
end

function (fac::ENCMapFactory)(self::AbsHyperComplex, p::Int, I::Tuple)
  cfac = chain_factory(self)
  A = cfac.A
  F = cfac.F
  n = ncols(A)
  k = nrows(A)
  i = first(I)
  dom = self[i]
  cod = self[i-1]

  if isone(i)
    # TODO: Do we need to take signs into account?
    dets = [det(A[:, data(ind)]) for ind in combinations(n, k)]
    return hom(dom, cod, [a*cod[1] for a in dets])
  end

  ext_dom = exterior_power(F, k + i - 1) # both are cached in F
  ext_cod = exterior_power(F, k + i - 2)

  img_gens = elem_type(cod)[]
  dom_dec = tensor_generator_decompose_function(dom)
  dom_sym_power = cfac.sym_powers[i]
  dom_ext_power = cfac.ext_powers[i]
  cod_sym_power = cfac.sym_powers[i-1]
  cod_ext_power = cfac.ext_powers[i-1]
  for g in gens(dom)
    img_gen = zero(cod)
    a, b = dom_dec(g)
    @assert parent(b) === dom_ext_power
    @assert parent(a) === dom_sym_power
    alpha = cfac.exps[i][first(coordinates(a))[1]]
    for (k, e) in enumerate(alpha)
      is_zero(e) && continue
      beta = copy(alpha)
      beta[k] -= 1
      aa = cod_sym_power[findfirst(==(beta), cfac.exps[i-1])]
      bb = _contract(b, A[k, :]; parent=cod_ext_power)
      img_gen = img_gen + tensor_pure_function(cod)((aa, bb))
    end
    push!(img_gens, img_gen)
  end
  return hom(dom, cod, img_gens)
end

function can_compute(fac::ENCMapFactory, self::AbsHyperComplex, p::Int, I::Tuple)
  is_one(p) || return false
  is_one(length(I)) || return false
  i = first(I)
  0 < i || return false
  cfac = chain_factory(self)
  return i <= ncols(cfac.A) - nrows(cfac.A) + 1
end

### The concrete struct
@attributes mutable struct EagonNorthcottComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType}
  internal_complex::HyperComplex{ChainType, MorphismType}

  function EagonNorthcottComplex(
      A::MatrixElem{T};
      F::FreeMod{T}=begin
        R = base_ring(A)
        F = is_graded(R) ? graded_free_module(R, [zero(grading_group(R)) for _ in 1:ncols(A)]) : FreeMod(R, ncols(A))
      end
    ) where {T}
    chain_fac = ENCChainFactory(A, F)
    map_fac = ENCMapFactory()

    internal_complex = HyperComplex(1, chain_fac, map_fac, [:chain];
                                    upper_bounds = Union{Int, Nothing}[ncols(A) - nrows(A) + 1],
                                    lower_bounds = Union{Int, Nothing}[0]
                                   )
    return new{FreeMod{T}, FreeModuleHom}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::EagonNorthcottComplex) = c.internal_complex

original_module(c::EagonNorthcottComplex) = chain_factory(c).F
matrix(c::EagonNorthcottComplex) = chain_factory(c).A

# user facing constructors
@doc raw"""
    eagon_northcott_complex(A::MatrixElem{T}; F::FreeMod{T}) where {T}

Given an ``m \times n``-matrix ``A`` over a commutative ring ``R``, construct
the Eagon-Northcott complex associated to it [EN62](@cite).

A free module ``F`` over the same ring as ``A`` can be passed so that the rows
of ``A`` are interpreted as elements in the dual of ``F`` in the construction
of the complex.
"""
function eagon_northcott_complex(
    A::MatrixElem{T};
    F::FreeMod{T}=begin
      R = base_ring(A)
      F = is_graded(R) ? graded_free_module(R, [zero(grading_group(R)) for _ in 1:ncols(A)]) : FreeMod(R, ncols(A))
    end
  ) where {T}
  return EagonNorthcottComplex(A; F)
end

# return the exterior power used for the i-th entry
function _exterior_power(c::EagonNorthcottComplex, i::Int)
  c[i] # fill the cache
  return chain_factory(c).ext_powers[i]
end

function _symmetric_power(c::EagonNorthcottComplex, i::Int)
  c[i] # fill the cache
  return chain_factory(c).sym_powers[i]
end
