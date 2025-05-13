### Production of the chains
struct InducedENCChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  enc::EagonNorthcottComplex
  I::ChainType
  ext_powers::Dict{Int, ChainType}

  function InducedENCChainFactory(enc::EagonNorthcottComplex, I::SubquoModule)
    @assert ambient_free_module(I) === original_module(enc)
    ChainType = typeof(I)
    return new{ChainType}(enc, I, Dict{Int, ChainType}())
  end
end

function (fac::InducedENCChainFactory)(self::AbsHyperComplex, ind::Tuple)
  i = first(ind)
  enc = fac.enc
  is_zero(i) && return sub(enc[0], [enc[0][1]])[1]
  F0 = original_module(enc)
  n = ngens(F0)
  k = nrows(matrix(enc))
  I0 = fac.I
  k + i - 1 > ngens(I0) && return sub(enc[i], elem_type(enc[i])[])[1]
  amb_ext_power = _exterior_power(enc, i)
  new_gens = elem_type(amb_ext_power)[]
  for as in combinations(gens(I0), k+i-1)
    push!(new_gens, wedge([repres(a) for a in as]; parent=amb_ext_power))
  end
  Ip, _ = sub(amb_ext_power, new_gens)
  fac.ext_powers[i] = Ip
  result = tensor_product(_symmetric_power(enc, i), Ip)
  return result
end

function can_compute(fac::InducedENCChainFactory, self::AbsHyperComplex, i::Tuple)
  return can_compute_index(fac.enc, i)
end

### Production of the morphisms
struct InducedENCMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  function InducedENCMapFactory(enc::EagonNorthcottComplex, I::SubquoModule)
    return new{SubQuoHom}()
  end
end

function (fac::InducedENCMapFactory)(self::AbsHyperComplex, p::Int, I::Tuple)
  i = first(I)
  cfac = chain_factory(self)
  enc = cfac.enc
  dom = self[i]
  cod = self[i-1]

  if i == 1 # evaluation of the determinants on the generators
    A = matrix(enc)
    R = base_ring(A)
    I0 = chain_factory(self).I
    B = matrix_space(R, ncols(A), ngens(I0))([repres(g)[i] for i in 1:ncols(A), g in gens(I0)])
    C = A*B
    dets = [det(C[:, data(ind)]) for ind in combinations(ncols(C), nrows(C))]
    return hom(dom, cod, [a*cod[1] for a in dets])
  end

  # TODO: This can probably be sped up by lifting directly on the exterior powers.
  img_gens = elem_type(cod)[]
  ambient_map = map(enc, i)
  for g in gens(dom)
    rep = repres(g)
    image_g = ambient_map(rep)
    push!(img_gens, cod(image_g))
  end
  return hom(dom, cod, img_gens)
end

function can_compute(fac::InducedENCMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  return can_compute_map(chain_factory(self).enc, p, i)
end

### The concrete struct
@attributes mutable struct InducedENC{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType}
  internal_complex::HyperComplex{ChainType, MorphismType}

  function InducedENC(enc::EagonNorthcottComplex, I::SubquoModule)
    chain_fac = InducedENCChainFactory(enc, I)
    map_fac = InducedENCMapFactory(enc, I)

    internal_complex = HyperComplex(1, chain_fac, map_fac, [:chain];
                                    upper_bounds=enc.internal_complex.upper_bounds,
                                    lower_bounds=enc.internal_complex.lower_bounds
                                   )
    return new{typeof(I), SubQuoHom}(internal_complex)
  end
end

### user facing constructors
@doc raw"""
    induced_eagon_northcott_complex(enc::EagonNorthcottComplex, I::SubquoModule)

Given an Eagon-Northcott complex `enc` for a matrix ``A`` on a free ``R``-module ``F``
there is a natural induced complex for any submodule ``I ⊂ F``.

Namely, if ``v`` is one of the rows of ``A``, interpreted as an element of ``Hom(F, R¹)``,
then contraction with ``v`` provides natural maps ``⋀ ᵖ F → ⋀ ᵖ⁻¹ F`` which restrict
to the submodules ``⋀ ᵖ I → ⋀ ᵖ⁻¹ I``. As the Eagon-Northcott complex is build up from
such contractions, we obtain a natural restricted complex which is constructed by this
method.
"""
function induced_eagon_northcott_complex(enc::EagonNorthcottComplex, I::SubquoModule)
  return InducedENC(enc, I)
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::InducedENC) = c.internal_complex
