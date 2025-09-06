########################################################################
# Tensor products of two complexes as a double complex
########################################################################

# As explained in `Types.jl`, we need to pass "factories" for the production of 
# the entries of a double complex and the morphisms between them to the 
# constructor. We first do the factory for the entries using the already 
# existing function for tensor products. 

# The concrete type contains all information necessary for the production.
mutable struct TensorProductFactory{ChainType} <: ChainFactory{ChainType}
  C1::ComplexOfMorphisms{<:ChainType}
  C2::ComplexOfMorphisms{<:ChainType}
end

# Then we override the call syntax as requested by the interface above.
function (fac::TensorProductFactory{<:ModuleFP})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  M = fac.C1[i]
  N = fac.C2[j]
  if iszero(M) || iszero(N)
    R = base_ring(M)
    is_graded(R) && return graded_free_module(R, 0)
    return FreeMod(R, 0)
  end
  return tensor_product(fac.C1[i], fac.C2[j])
end

function can_compute(fac::TensorProductFactory, dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  return i in range(fac.C1) && j in range(fac.C2)
end

# Now we move on to the factories for the maps. Here we have to provide 
# either one for the horizontal and the vertical morphisms.
mutable struct VerticalTensorMapFactory{MapType} <: ChainMorphismFactory{MapType}
  C1::ComplexOfMorphisms
  C2::ComplexOfMorphisms
end

# Again, we override the call syntax as required by the interface. Note that 
# here we need to make use of the first argument, the double complex itself, 
# in order to access the correct domains and codomains. 
#
# On input (dc, i, j) this must produce the morphism 
#
#   id ⊗ ψⱼ : Cᵢ⊗ Dⱼ → Cᵢ⊗ Dⱼ₊₋₁
#
# induced by the (co-)boundary map ψⱼ : Dⱼ → Dⱼ₊₋1 of the second complex 
# with the sign depending on the `vertical_direction` of `dc`.
function (fac::VerticalTensorMapFactory{<:ModuleFPHom})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  dom = dc[i, j]
  cod_ind = (i, j + (fac.C2.typ == :chain ? -1 : +1))
  cod = dc[cod_ind]
  (iszero(dom) || iszero(cod)) && return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)])
  return tensor_product(dom, cod, [id_hom(fac.C1[i]), map(fac.C2, j)])
end

function can_compute(fac::VerticalTensorMapFactory, dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  return i in range(fac.C1) && j in map_range(fac.C2)
end

# Same for the horizontal maps.
# On input (dc, i, j) this must produce the morphism 
#
#   φᵢ ⊗ id : Cᵢ⊗ Dⱼ → Cᵢ₊₋₁⊗ Dⱼ 
#
# induced by the (co-)boundary map φᵢ : Cᵢ → Cᵢ₊₋1 of the first complex 
# with the sign depending on the `horizontal_direction` of `dc`.
mutable struct HorizontalTensorMapFactory{MapType} <: ChainMorphismFactory{MapType}
  C1::ComplexOfMorphisms
  C2::ComplexOfMorphisms
end

function (fac::HorizontalTensorMapFactory{<:ModuleFPHom})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  dom = dc[i, j]
  cod_ind = (i + (fac.C1.typ == :chain ? -1 : +1), j)
  cod = dc[cod_ind]
  (iszero(dom) || iszero(cod)) && return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)])
  return tensor_product(dom, cod, [map(fac.C1, i), id_hom(fac.C2[j])])
end

function can_compute(fac::HorizontalTensorMapFactory, dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  return i in map_range(fac.C1) && j in range(fac.C2)
end

# The user facing constructor for the tensor product of complexes 
# now takes the following rather simple form:
@doc raw"""
    tensor_product(C::ComplexOfMorphisms{ChainType}, D::ComplexOfMorphisms{ChainType}) where {ChainType}

Create the tensor product of two complexes `C` and `D` as a double complex.

In order for the generic implementation to work for your specific `ChainType` the following 
needs to be implemented.

  * `morphism_type(ChainType)` must produce the type of morphisms between objects of type `ChainType`;
  * the call signature for `function (fac::TensorProductFactory{ChainType})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)` needs to be overwritten for your specific instance of `ChainType` to produce the `(i, j)`-th entry of the double complex, i.e. the tensor product of `C[i]` and `D[j]`;
  * the call signature for `function (fac::HorizontalTensorMapFactory{ChainType})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)` needs to be overwritten to produce the map on tensor products `C[i] ⊗ D[j] → C[i±1] ⊗ D[j]` induced by the (co-)boundary map on `C` (the sign depending on the `typ` of `C`);
  * similarly for the `VerticalTensorMapFactory`.


See the file `experimental/DoubleComplex/src/tensor_products.jl` for examples.
"""
function tensor_product(C1::ComplexOfMorphisms{ChainType}, C2::ComplexOfMorphisms{ChainType}) where {ChainType}
  MorphismType = morphism_type(ChainType)
  result = DoubleComplexOfMorphisms(TensorProductFactory{ChainType}(C1, C2),
                                    HorizontalTensorMapFactory{MorphismType}(C1, C2),
                                    VerticalTensorMapFactory{MorphismType}(C1, C2),
                                    horizontal_direction=typ(C1),
                                    vertical_direction=typ(C2)
                                   )
  # TODO: Find out how to access the information that a complex can not extend.
  #(is_complete(C1) && is_complete(C2)) || error("tensor product is only implemented for complete complexes at the moment")
  r1 = range(C1)
  result.right_bound = (typ(C1) == :chain ? first(r1) : last(r1))
  result.left_bound = (typ(C1) == :cochain ? first(r1) : last(r1))
  r2 = range(C2)
  result.upper_bound = (typ(C2) == :chain ? first(r2) : last(r2))
  result.lower_bound = (typ(C2) == :cochain ? first(r2) : last(r2))

  return result
end

### missing functionality
# A fallback in case nothing more specific is known
morphism_type(::Type{T}) where {T<:ModuleFP} = ModuleFPHom

zero(phi::ModuleFPHom) = hom(domain(phi), codomain(phi), [zero(codomain(phi)) for i in 1:ngens(domain(phi))])
zero_morphism(dom::ModuleFP, cod::ModuleFP) = hom(dom, cod, [zero(cod) for i in 1:ngens(dom)])


########################################################################
# Tensor products of more than two chain complexes as a hypercomplex
########################################################################
struct TensorProductChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  list::Vector{<:ComplexOfMorphisms}

  function TensorProductChainFactory(list::Vector{T}) where {ChainType, T<:ComplexOfMorphisms{ChainType}}
    return new{ChainType}(list)
  end
end

function (fac::TensorProductChainFactory{T})(HC::AbsHyperComplex, i::Tuple) where {T<:ModuleFP}
  d = length(i)
  if all(k->i[k] in range(fac.list[k]), 1:d)
    return tensor_product([fac.list[k][i[k]] for k in 1:d])
  end
  C = first(list)
  k = first(range(C))
  R = base_ring(C[k])
  return FreeMod(R, 0)
end

function can_compute(fac::TensorProductChainFactory, HC::AbsHyperComplex, i::Tuple)
  return true
end

struct TensorProductMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  list::Vector{<:ComplexOfMorphisms}

  function TensorProductMapFactory(list::Vector{T}) where {ChainType, T<:ComplexOfMorphisms{ChainType}}
    MorphismType = morphism_type(ChainType)
    return new{MorphismType}(list)
  end
end

function (fac::TensorProductMapFactory{T})(HC::AbsHyperComplex, p::Int, i::Tuple) where {T<:ModuleFPHom}
  d = length(i)
  v = collect(i)
  v[p] = (direction(HC, p) == :chain ? v[p] - 1 : v[p] + 1)
  j = Tuple(v)
  if all(k->i[k] in range(fac.list[k]), 1:d) && all(k->j[k] in range(fac.list[k]), 1:d)
    map_vec = T[id_hom(fac.list[k][i[k]]) for k in 1:p-1]
    push!(map_vec, map(fac.list[p], i[p]))
    map_vec = vcat(map_vec, T[id_hom(fac.list[k][i[k]]) for k in p+1:d])
    return hom_tensor(HC[i], HC[j], map_vec)
  end
  dom = HC[i]
  cod = HC[j]
  return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)])
end

function can_compute(fac::TensorProductMapFactory, HC::AbsHyperComplex, p::Int, i::Tuple)
  return true
end

function tensor_product(list::Vector{T}) where {T<:ComplexOfMorphisms}
  d = length(list)
  upper_bounds = [(typ(C) == :chain ? last(range(C)) : first(range(C))) for C in list]
  lower_bounds = [(typ(C) == :cochain ? last(range(C)) : first(range(C))) for C in list]
  typs = [typ(C) for C in list]

  chain_factory = TensorProductChainFactory(list)
  map_factory = TensorProductMapFactory(list)
  return HyperComplex(d, chain_factory, map_factory, typs, upper_bounds=upper_bounds, 
                      lower_bounds=lower_bounds
                     )
end

########################################################################
# Tensor products of hypercomplexes                                    #
########################################################################
struct HCTensorProductChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  factors::Vector{<:AbsHyperComplex}
  mult_map_cache::Dict{<:Tuple, <:Map}

  function HCTensorProductChainFactory(
      ::Type{CT},
      factors::Vector{<:AbsHyperComplex}
    ) where {CT}
    mult_map_cache = Dict{Tuple, Map}()
    return new{CT}(factors, mult_map_cache)
  end
end

function (fac::HCTensorProductChainFactory{ChainType})(c::AbsHyperComplex, I::Tuple) where {ChainType}
  i = collect(I)
  # decompose I into its individual indices
  j = Vector{Vector{Int}}()
  k = 0
  for f in fac.factors
    push!(j, i[k+1:k+dim(f)])
    k = k + dim(f)
  end
  factors = [fac.factors[k][Tuple(a)] for (k, a) in enumerate(j)]
  #any(iszero, factors) && return zero_object(first(factors))[1]
  tmp_res = _tensor_product(factors...)
  @assert tmp_res isa Tuple{<:ChainType, <:Map} "the output of `tensor_product` does not have the anticipated format; see the source code for details"
  # If you got here because of the error message above:
  # This is supposed to be generic code and it attempts to build the tensor product 
  # of `factors` using the method of `tensor_product` for your specific type. 
  # The convention in Oscar is that this should return a pair `(M, mult)` consisting 
  # of the actual tensor product `M` and a multiplication map `mult` which takes a 
  # `Tuple` of elements in the `factors` to its tensor product in `M`. 
  # Unfortunately, it can not be assumed that the original method for `tensor_product` 
  # produces this output and it is an ongoing effort to streamline this throughout 
  # Oscar. Since such changes tend to happen with severe delay, if ever, we provide 
  # a little workaround here. If this code does not run for your type of chains, 
  # you may try two things:
  #
  #   1) Overwrite the method for `_tensor_product` below for your type and wrap the 
  #      original method for `tensor_product` so that the anticipated output is produced
  #
  #   2) If that does not work, you may also overwrite the whole method for production 
  #      of the tensor products here.
  result, mult_map = tmp_res
  fac.mult_map_cache[I] = mult_map
  return result
end

# By default, we assume that `tensor_product` produces the correct output
function _tensor_product(v::Any...)
  return tensor_product(v...)
end

# If not, we can wrap it here so that the correct output is produced
function _tensor_product(v::ModuleFP...)
  return tensor_product(v...; task=:with_map)
end

function can_compute(fac::HCTensorProductChainFactory, c::AbsHyperComplex, I::Tuple)
  i = collect(I)

  # decompose I into its individual indices
  j = Vector{Vector{Int}}()
  k = 0
  for f in fac.factors
    push!(j, i[k+1:k+dim(f)])
    k = k + dim(f)
  end
  return all(k->can_compute_index(fac.factors[k], Tuple(j[k])), 1:length(j))
end

### Production of the morphisms 
struct HCTensorProductMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  factors::Vector{<:AbsHyperComplex}

  function HCTensorProductMapFactory(
      ::Type{MT},
      factors::Vector{<:AbsHyperComplex}
    ) where {MT}
    return new{MT}(factors)
  end
end

function zero_map(M::ModuleFP, N::ModuleFP)
  base_ring(M) === base_ring(N) || error("base rings must coincide")
  return hom(M, N, elem_type(N)[zero(N) for i in 1:ngens(M)], check=false)
end

function (fac::HCTensorProductMapFactory{MorphismType})(c::AbsHyperComplex, p::Int, I::Tuple) where {MorphismType}
  next = collect(I)
  next[p] = next[p] + (direction(c, p) == :chain ? -1 : 1)
  J = Tuple(next)

  if iszero(c[I]) || iszero(c[J])
    return zero_map(c[I], c[J])
  end
  
  i = collect(I)
  # decompose I into its individual indices
  j = Vector{Vector{Int}}()
  k = 0
  for f in fac.factors
    push!(j, i[k+1:k+dim(f)])
    k = k + dim(f)
  end
  maps = MorphismType[]
  d = length.(j)
  k = findfirst(k->sum(d[1:k]; init=0) >= p, 1:length(j))
  k === nothing && error("direction out of bounds")
  q = p - sum(d[1:k-1]; init = 0)
  maps = MorphismType[id_hom(fac.factors[l][Tuple(j[l])]) for l in 1:k-1]
  push!(maps, map(fac.factors[k], q, Tuple(j[k]))) 
  maps = vcat(maps, MorphismType[id_hom(fac.factors[l][Tuple(j[l])]) for l in k+1:length(j)])
  return tensor_product(c[I], c[J], maps)
end

### The concrete struct
@attributes mutable struct HCTensorProductComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}
  factors::Vector{<:AbsHyperComplex}

  function HCTensorProductComplex(
      factors::Vector{<:AbsHyperComplex}
    )
    CT = reduce(typejoin, chain_type.(factors))
    MT = reduce(typejoin, morphism_type.(factors))
    chain_fac = HCTensorProductChainFactory(CT, factors)
    map_fac = HCTensorProductMapFactory(MT, factors)

    d = sum(dim(c) for c in factors; init=0)
    dir = vcat([Symbol[direction(c, p) for p in 1:dim(c)] for c in factors]...)
    upper_bounds = vcat([Union{Int, Nothing}[(has_upper_bound(c, p) ? upper_bound(c, p) : nothing) for p in 1:dim(c)] for c in factors]...)
    lower_bounds = vcat([Union{Int, Nothing}[(has_lower_bound(c, p) ? lower_bound(c, p) : nothing) for p in 1:dim(c)] for c in factors]...)

    internal_complex = HyperComplex(d, chain_fac, map_fac, 
                                    dir, upper_bounds=upper_bounds,
                                    lower_bounds=lower_bounds
                                   )
    # Assuming that ChainType and MorphismType are provided by the input
    return new{CT, MT}(internal_complex, factors)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::HCTensorProductComplex) = c.internal_complex

### Additional functionality
factors(c::HCTensorProductComplex) = c.factors

### User facing constructor
function tensor_product(factors::Vector{<:AbsHyperComplex})
  return HCTensorProductComplex(factors)
end

function tensor_product(c::AbsHyperComplex, cs::AbsHyperComplex...)
  return tensor_product([c, cs...])
end

function tensor_product(M::ModuleFP{T}, Ms::ModuleFP{T}...) where {U<:MPolyComplementOfPrimeIdeal, T<:MPolyLocRingElem{<:Any, <:Any, <:Any, <:Any, U}}
  R = base_ring(M)
  @assert all(N->base_ring(N)===R, Ms) "modules must be defined over the same ring"
  return tensor_product([free_resolution(SimpleFreeResolution, N)[1] for N in (M, Ms...)])
end
