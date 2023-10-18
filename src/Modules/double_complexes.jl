### Abstract type and interface for double complexes
abstract type AbsDoubleComplexOfMorphisms{ChainType, MapType} end

### asking for the entries of the complex
getindex(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = underlying_double_complex(dc)[i, j]

### asking for horizontal and vertical maps
horizontal_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = horizontal_map(underlying_double_complex(dc), i, j)
vertical_map(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) = vertical_map(underlying_double_complex(dc), i, j)
vertical_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = vertical_map(dc, t...)
horizontal_map(dc::AbsDoubleComplexOfMorphisms, t::Tuple) = horizontal_map(dc, t...)


### asking for the typ of the row and the column complexes
horizontal_typ(dc::AbsDoubleComplexOfMorphisms) = horizontal_typ(underlying_double_complex(dc))
vertical_typ(dc::AbsDoubleComplexOfMorphisms) = vertical_typ(underlying_double_complex(dc))

### asking for known bounds of the row and column complexes
# These are also used to filter legitimate requests to entries and maps, 
# i.e. one can not ask for an entry or a map of a double complex outside 
# of these bounds. If no bound is given in a specific direction, then every 
# request is considered legitimate there. 
#
# Note that at the moment only rectangular shapes are available to describe bounds.
is_horizontally_bounded(dc::AbsDoubleComplexOfMorphisms) = has_horizontal_upper_bound(dc) && has_horizontal_lower_bound(dc)
is_vertically_bounded(dc::AbsDoubleComplexOfMorphisms) = has_vertical_upper_bound(dc) && has_vertical_lower_bound(dc)
is_bounded(dc::AbsDoubleComplexOfMorphisms) = is_horizontally_bounded(dc) && is_vertically_bounded(dc)

has_horizontal_upper_bound(D::AbsDoubleComplexOfMorphisms) = has_horizontal_upper_bound(underlying_double_complex(D))
has_horizontal_lower_bound(D::AbsDoubleComplexOfMorphisms) = has_horizontal_lower_bound(underlying_double_complex(D))
has_vertical_upper_bound(D::AbsDoubleComplexOfMorphisms)   = has_vertical_upper_bound(underlying_double_complex(D))
has_vertical_lower_bound(D::AbsDoubleComplexOfMorphisms)   = has_vertical_lower_bound(underlying_double_complex(D))

is_complete(dc::AbsDoubleComplexOfMorphisms) = is_complete(underlying_double_complex(dc))

# The concrete architecture of double complexes is lazy by default. 
# Hence the constructor needs to be provided with the means to produce 
# the entries of a double complex on request. This is achieved by passing
# certain "factories" to the constructor which carry out this production 
# of the entries and the maps on request. In what follows we specify 
# abstract types for these factories and their interface.

# An abstract type to produce the chains in the (i,j)-th entry of a double
# complex. The interface is formulated for this abstract type, but the 
# user needs to implement a concrete type for any concrete implementation 
# of a double complex.
abstract type ChainFactory{ChainType} end

# Produce the t = (i, j)-th entry which will then be cached.
function (fac::ChainFactory)(dc::AbsDoubleComplexOfMorphisms, t::Tuple)
  return fac(dc, t...)
end

# A dummy placeholder which must be overwritten.
# The first argument will always be the actual double complex itself, 
# so that the body of the function has access to all data already generated 
# and the other functionality available to this double complex. 
function (fac::ChainFactory{ChainType})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)::ChainType where {ChainType}
  error("production of the ($i, $j)-th chain not implemented")
end

# An abstract type to produce the chain maps going out of the 
# (i, j)-th entry either in the vertical or the horizontal direction.
# The user needs to implement a concrete instance which then knows 
# in particular whether it's supposed to produce the vertical 
# or the horizontal maps (which are then to be cached).
abstract type ChainMorphismFactory{MorphismType} end

function (fac::ChainMorphismFactory)(dc::AbsDoubleComplexOfMorphisms, t1::Tuple)
  return fac(dc, t1...)
end

# A dummy placeholder which must be overwritten; see below.
function (fac::ChainMorphismFactory{MorphismType})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)::MorphismType where {MorphismType}
  error("could not construct morphism from ($i, $j)")
end

# A minimal concrete type realizing a double complex.
#
# The design is lazy by default. All entries are produced on 
# request and then cached in dictionaries. For the production 
# the user has to provide "factories" in the above sense. 
mutable struct DoubleComplexOfMorphisms{ChainType, MorphismType<:Map} <: AbsDoubleComplexOfMorphisms{ChainType, MorphismType}
  chains::Dict{Tuple{Int, Int}, <:ChainType}
  horizontal_maps::Dict{Tuple{Int, Int}, <:MorphismType}
  vertical_maps::Dict{Tuple{Int, Int}, <:MorphismType}

  # Possible functions to produce new chains and maps in case of an incomplete complex
  chain_factory::ChainFactory{<:ChainType}
  horizontal_map_factory::ChainMorphismFactory{<:MorphismType}
  vertical_map_factory::ChainMorphismFactory{<:MorphismType}

  # Information about the nature of the complex
  is_complete::Bool
  horizontal_typ::Symbol
  vertical_typ::Symbol

  # Information about boundedness and completeness of the complex
  horizontal_upper_bound::Int
  horizontal_lower_bound::Int
  vertical_upper_bound::Int
  vertical_lower_bound::Int

  function DoubleComplexOfMorphisms(
      chain_factory::ChainFactory{ChainType}, 
      horizontal_map_factory::ChainMorphismFactory{MorphismType}, 
      vertical_map_factory::ChainMorphismFactory{MorphismType};
      horizontal_typ::Symbol=:chain,
      vertical_typ::Symbol=:chain,
      horizontal_upper_bound::Union{Int, Nothing} = nothing,
      horizontal_lower_bound::Union{Int, Nothing} = nothing,
      vertical_upper_bound::Union{Int, Nothing} = nothing,
      vertical_lower_bound::Union{Int, Nothing} = nothing,
      is_complete::Bool=false
    ) where {ChainType, MorphismType}
    result = new{ChainType, MorphismType}(Dict{Tuple{Int, Int}, ChainType}(),
                                          Dict{Tuple{Int, Int}, MorphismType}(),
                                          Dict{Tuple{Int, Int}, MorphismType}(),
                                          chain_factory, horizontal_map_factory, 
                                          vertical_map_factory, is_complete, 
                                          horizontal_typ, vertical_typ)
    horizontal_upper_bound !== nothing && (result.horizontal_upper_bound = horizontal_upper_bound)
    horizontal_lower_bound !== nothing && (result.horizontal_lower_bound = horizontal_lower_bound)
    vertical_upper_bound !== nothing && (result.vertical_upper_bound = vertical_upper_bound)
    vertical_lower_bound !== nothing && (result.vertical_lower_bound = vertical_lower_bound)
    return result
  end
end

### Basic getters
function is_complete(D::DoubleComplexOfMorphisms)
  return D.is_complete
end

function horizontal_typ(D::DoubleComplexOfMorphisms)
  return D.horizontal_typ
end

function vertical_typ(D::DoubleComplexOfMorphisms)
  return D.vertical_typ
end

function chain_cache(D::DoubleComplexOfMorphisms) 
  return D.chains
end

function horizontal_map_cache(D::DoubleComplexOfMorphisms)
  return D.horizontal_maps
end

function vertical_map_cache(D::DoubleComplexOfMorphisms)
  return D.vertical_maps
end

### Boundedness
has_horizontal_upper_bound(D::DoubleComplexOfMorphisms) = isdefined(D, :horizontal_upper_bound)
has_horizontal_lower_bound(D::DoubleComplexOfMorphisms) = isdefined(D, :horizontal_lower_bound)
has_vertical_upper_bound(D::DoubleComplexOfMorphisms)   = isdefined(D, :vertical_upper_bound)
has_vertical_lower_bound(D::DoubleComplexOfMorphisms)   = isdefined(D, :vertical_lower_bound)

function horizontal_range(D::DoubleComplexOfMorphisms) 
  is_horizontally_bounded(D) || error("complex is not known to be horizontally bounded")
  if horizontal_typ(D) == :chain
    return horizontal_upper_bound(D):-1:horizontal_lower_bound(D)
  elseif horizontal_typ(D) == :cochain
    return horizontal_lower_bound(D):horizontal_upper_bound(D)
  end
  error("typ not recognized")
end

function vertical_range(D::DoubleComplexOfMorphisms) 
  is_vertically_bounded(D) || error("complex is not known to be vertically bounded")
  if vertical_typ(D) == :chain
    return vertical_upper_bound(D):-1:vertical_lower_bound(D)
  elseif vertical_typ(D) == :cochain
    return vertical_lower_bound(D):vertical_upper_bound(D)
  end
  error("typ not recognized")
end

function horizontal_upper_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :horizontal_upper_bound) || error("complex has no known horizontal upper bound")
  return D.horizontal_upper_bound
end

function horizontal_lower_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :horizontal_lower_bound) || error("complex has no known horizontal lower bound")
  return D.horizontal_lower_bound
end

function vertical_upper_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :vertical_upper_bound) || error("complex has no known vertical upper bound")
  return D.vertical_upper_bound
end

function vertical_lower_bound(D::DoubleComplexOfMorphisms) 
  isdefined(D, :vertical_lower_bound) || error("complex has no known vertical lower bound")
  return D.vertical_lower_bound
end

### User facing functionality
function getindex(D::DoubleComplexOfMorphisms, t::Tuple)
  return getindex(D, t...)
end

function getindex(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  (has_horizontal_upper_bound(D) && i<=horizontal_upper_bound(D)) || error("index out of range")
  (has_horizontal_lower_bound(D) && i>=horizontal_lower_bound(D)) || error("index out of range")
  (has_vertical_upper_bound(D) && j<=vertical_upper_bound(D)) || error("index out of range")
  (has_vertical_lower_bound(D) && j>=vertical_lower_bound(D)) || error("index out of range")
  haskey(D.chains, (i, j)) && return D.chains[(i, j)]
  new_chain = D.chain_factory(D, i, j)
  D.chains[(i, j)] = new_chain
  return new_chain
end

function vertical_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  (has_horizontal_upper_bound(D) && i<=horizontal_upper_bound(D)) || error("index out of range")
  (has_horizontal_lower_bound(D) && i>=horizontal_lower_bound(D)) || error("index out of range")
  if horizontal_typ(D) == :chain
    (has_vertical_upper_bound(D) && j<=vertical_upper_bound(D)) || error("index out of range")
    (has_vertical_lower_bound(D) && j>vertical_lower_bound(D)) || error("index out of range")
  else
    (has_vertical_upper_bound(D) && j<vertical_upper_bound(D)) || error("index out of range")
    (has_vertical_lower_bound(D) && j>=vertical_lower_bound(D)) || error("index out of range")
  end
  haskey(vertical_map_cache(D), (i, j)) && return vertical_map_cache(D)[(i, j)]
  new_map = D.vertical_map_factory(D, i, j)
  vertical_map_cache(D)[(i, j)] = new_map
  return new_map
end

function horizontal_map(D::DoubleComplexOfMorphisms, i::Int, j::Int)
  (has_vertical_upper_bound(D) && j<=vertical_upper_bound(D)) || error("index out of range")
  (has_vertical_lower_bound(D) && j>=vertical_lower_bound(D)) || error("index out of range")
  if horizontal_typ(D) == :chain
    (has_horizontal_upper_bound(D) && i<=horizontal_upper_bound(D)) || error("index out of range")
    (has_horizontal_lower_bound(D) && i>horizontal_lower_bound(D)) || error("index out of range")
  else
    (has_horizontal_upper_bound(D) && i<horizontal_upper_bound(D)) || error("index out of range")
    (has_horizontal_lower_bound(D) && i>=horizontal_lower_bound(D)) || error("index out of range")
  end
  haskey(horizontal_map_cache(D), (i, j)) && return horizontal_map_cache(D)[(i, j)]
  new_map = D.horizontal_map_factory(D, i, j)
  horizontal_map_cache(D)[(i, j)] = new_map
  return new_map
end


### Concrete implementation for tensor products of complexes

# As explained above, we need to pass "factories" for the production of 
# the entries of a double complex and the morphisms between them to the 
# constructor. We first do the factory for the entries using the already 
# existing function for tensor products. 

# The concrete type contains all information necessary for the production.
mutable struct TensorProductFactory{ChainType} <: ChainFactory{ChainType}
  C1::ComplexOfMorphisms{<:ChainType}
  C2::ComplexOfMorphisms{<:ChainType}
end

# Then we override the call syntax as requested by the interface above.
function (fac::TensorProductFactory{<:ModuleFP})(dc::DoubleComplexOfMorphisms, i::Int, j::Int)
  M = fac.C1[i]
  N = fac.C2[j]
  if iszero(M) || iszero(N)
    R = base_ring(M)
    is_graded(R) && return graded_free_module(R, 0)
    return FreeMod(R, 0)
  end
  return tensor_product(fac.C1[i], fac.C2[j])
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
function (fac::VerticalTensorMapFactory{<:ModuleFPHom})(dc::DoubleComplexOfMorphisms, i::Int, j::Int)
  dom = dc[i, j]
  cod_ind = (i, j + (fac.C2.typ == :chain ? -1 : +1))
  cod = dc[cod_ind]
  (iszero(dom) || iszero(cod)) && return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)])
  return tensor_product(dom, cod, [identity_map(fac.C1[i]), map(fac.C2, j)])
end

# Same for the horizontal maps.
mutable struct HorizontalTensorMapFactory{MapType} <: ChainMorphismFactory{MapType}
  C1::ComplexOfMorphisms
  C2::ComplexOfMorphisms
end

function (fac::HorizontalTensorMapFactory{<:ModuleFPHom})(dc::DoubleComplexOfMorphisms, i::Int, j::Int)
  dom = dc[i, j]
  cod_ind = (i + (fac.C1.typ == :chain ? -1 : +1), j)
  cod = dc[cod_ind]
  (iszero(dom) || iszero(cod)) && return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)])
  return tensor_product(dom, cod, [map(fac.C1, i), identity_map(fac.C2[j])])
end

# The user facing constructor for the tensor product of complexes 
# now takes the following rather simple form:
function tensor_product(C1::ComplexOfMorphisms{<:ModuleFP}, C2::ComplexOfMorphisms{<:ModuleFP})
  result = DoubleComplexOfMorphisms(TensorProductFactory{ModuleFP}(C1, C2),
                                    HorizontalTensorMapFactory{ModuleFPHom}(C1, C2),
                                    VerticalTensorMapFactory{ModuleFPHom}(C1, C2),
                                    horizontal_typ=typ(C1),
                                    vertical_typ=typ(C2)
                                   )
  r1 = range(C1)
  result.horizontal_upper_bound = (typ(C1) == :chain ? first(r1) : last(r1))
  result.horizontal_lower_bound = (typ(C1) == :cochain ? first(r1) : last(r1))
  r2 = range(C2)
  result.vertical_upper_bound = (typ(C2) == :chain ? first(r2) : last(r2))
  result.vertical_lower_bound = (typ(C2) == :cochain ? first(r2) : last(r2))

  if is_complete(C1) && is_complete(C2)
    result.is_complete = true
  end
  return result
end

### Generic functionality
function total_complex(dc::DoubleComplexOfMorphisms{ChainType, MorphismType}) where {ChainType, MorphismType}
  is_bounded(dc) || error("computation of total complexes is only implemented for bounded double complexes")
  vertical_typ(dc) == horizontal_typ(dc) || error("horizontal and vertical typs must be the same")
  if vertical_typ(dc) == :chain
    return _total_chain_complex(dc)
  else
    return _total_cochain_complex(dc)
  end
end

function _total_chain_complex(dc::DoubleComplexOfMorphisms{ChainType, MorphismType}) where {ChainType, MorphismType}
  r1 = horizontal_range(dc)
  r2 = vertical_range(dc)
  r_tot = (first(r1) + first(r2)):-1:(last(r1) + last(r2))
  new_chains = ChainType[]
  new_maps = MorphismType[]

  # First round for initialization
  index_pairs = [I for I in Iterators.product(r1, r2) if sum(I) == first(r_tot)]
  summands = [dc[I] for I in index_pairs]
  new_chain, inc, pr = direct_sum(summands...; task=:both)
  last_inc = inc
  last_pr = pr
  last_chain = new_chain
  last_index_pairs = index_pairs
  push!(new_chains, new_chain)
  for k in r_tot
    k == first(r_tot) && continue
    index_pairs = [I for I in Iterators.product(r1, r2) if sum(I) == k]
    summands = [dc[I] for I in index_pairs]
    new_chain, inc, pr = direct_sum(summands...; task=:both)
    @assert all(f->domain(f) === new_chain, pr)
    @assert all(f->codomain(f) === new_chain, inc)
    push!(new_chains, new_chain)
    # assemble the map
    boundary_map = hom(last_chain, new_chain, elem_type(new_chain)[zero(new_chain) for i in 1:ngens(last_chain)])
    for (l, I) in enumerate(last_index_pairs)
      p = last_pr[l]
      @assert domain(p) === last_chain

      if I[2] > vertical_lower_bound(dc)
        vert = compose(p, vertical_map(dc, I))
        m = findfirst(x->x==(I[1], I[2]-1), index_pairs)
        @assert m !== nothing
        @assert codomain(vert) === dc[I[1], I[2]-1] === domain(inc[m])
        boundary_map = boundary_map + compose(vert, inc[m])
      end

      if I[1] > horizontal_lower_bound(dc)
        horz = compose(p, horizontal_map(dc, I))
        n = findfirst(x->x==(I[1]-1, I[2]), index_pairs)
        @assert n !== nothing
        @assert codomain(horz) === dc[I[1]-1, I[2]] === domain(inc[n])
        boundary_map = boundary_map + (-1)^I[2]*compose(horz, inc[n])
      end
    end
    push!(new_maps, boundary_map)
    last_index_pairs = index_pairs
    last_chain = new_chain
    last_inc = inc
    last_pr = pr
  end
  return ComplexOfMorphisms(ChainType, new_maps, seed=last(r_tot))
end

function Base.:*(k::Int, f::ModuleFPHom)
  return base_ring(codomain(f))(k)*f
end

function _total_cochain_complex(dc::DoubleComplexOfMorphisms{ChainType, MorphismType}) where {ChainType, MorphismType}
  error("total complex of double cochain complexes currently not implemented")
end

### Missing functionality for complexes
typ(C::ComplexOfMorphisms) = C.typ

### Missing functionality for tensor products
#
# Actually, this was not really missing, but available under a different name. 
# I leave the code snippets here for the moment.
function is_tensor_product(M::ModuleFP)
  !has_attribute(M, :tensor_product) && return false, [M]
  return true, get_attribute(M, :tensor_product)::Tuple
end

function tensor_pure_function(M::ModuleFP)
  success, facs = is_tensor_product(M)
  success || error("not a tensor product")
  return get_attribute(M, :tensor_pure_function)
end

function tensor_generator_decompose_function(M::ModuleFP)
  success, facs = is_tensor_product(M)
  success || error("not a tensor product")
  return get_attribute(M, :tensor_generator_decompose_function)
end
  
function tensor_product(mod::Vector{<:ModuleFP})
  return tensor_product(mod...)
end

function tensor_product(f::ModuleFPHom...)
  return tensor_product(collect(f))
end

function tensor_product(maps::Vector{<:ModuleFPHom};
    domain::ModuleFP=tensor_product(domain.(maps)),
    codomain::ModuleFP=tensor_product(codomain.(maps))
  )
  n = length(maps)
  iszero(n) && error("list of maps must not be empty")
  R = base_ring(domain)
  S = base_ring(codomain)
  R === S || error("tensor product of maps with base change not implemented")
  success, dom_facs = is_tensor_product(domain)
  n == length(dom_facs) || error("number of factors is incompatible")
  !success && error("domain must be a tensor product")
  all(k->dom_facs[k] === Oscar.domain(maps[k]), 1:n) || error("domains not compatible")
  success, cod_facs = is_tensor_product(codomain)
  n == length(cod_facs) || error("number of factors is incompatible")
  !success && error("codomain must be a tensor product")
  all(k->cod_facs[k] === Oscar.codomain(maps[k]), 1:n)|| error("codomains not compatible")

  pure_func_dom = tensor_pure_function(domain)
  inv_pure_func_dom = tensor_generator_decompose_function(domain)
  
  pure_func_cod = tensor_pure_function(codomain)
  inv_pure_func_cod = tensor_generator_decompose_function(codomain)

  # Compute the images of the generators
  img_gens = elem_type(codomain)[]
  ngens_cod_fac = ngens.(cod_facs)
  cod_ranges = Iterators.product([1:r for r in ngens_cod_fac]...)
  # Iterate over the generators of the domain
  for e in gens(domain)
    # decompose it into a product of elementary tensors
    x = inv_pure_func_dom(e)
    # map each factor
    y = [maps[k](x[k]) for k in 1:n]
    # initialize a variable for the image
    z = zero(codomain)
    # iterate over all the tuples of indices
    for I in cod_ranges
      # get the corresponding tuple of factors for this collection of generators
      v = Tuple([cod_facs[k][I[k]] for k in 1:n])
      # map it to its pure tensor
      w = pure_func_cod(v)
      # gather the product of the coefficients
      coeff = prod(y[k][I[k]] for k in 1:n; init=one(S))
      # add it to the intermediate result
      z = z + coeff * w
    end
    push!(img_gens, z)
  end
  return hom(domain, codomain, img_gens)
end

tensor_product(dom::ModuleFP, cod::ModuleFP, maps::Vector{<:ModuleFPHom}) = tensor_product(maps, domain=dom, codomain=cod)



