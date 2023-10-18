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
is_complete(C::ComplexOfMorphisms) = C.complete

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

