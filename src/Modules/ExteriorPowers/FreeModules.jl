########################################################################
# Exterior powers of free modules
#
# For F = Rⁿ we provide methods to create ⋀ ᵖF for arbitrary 0 ≤ p ≤ n.
# These modules are cached in F and know that they are an exterior 
# power of F. This allows us to implement the wedge product of their 
# elements. 
########################################################################

# User facing constructor for ⋀ ᵖ F.
function exterior_power(F::FreeMod, p::Int; cached::Bool=true)
  @req 0 <= p <= rank(F) "exponent out of bounds"

  if cached
    powers = _exterior_powers(F)
    haskey(powers, p) && return powers[p]
  end

  R = base_ring(F)
  n = rank(F)
  result_ = FreeMod(R, binomial(n, p))

  # In case F was graded, we have to take an extra detour. 
  result = if is_graded(F)
    G = grading_group(F)
    weights = elem_type(G)[]
    for ind in OrderedMultiIndexSet(p, n)
      push!(weights, sum(_degree_fast(F[i]) for i in indices(ind); init=zero(G)))
    end
    grade(result_, weights)
  else
    result_
  end

  # Create the multiplication map
  function my_mult(u::Tuple{Vararg{FreeModElem}})
    isempty(u) && return result[1] # only the case p=0
    @req all(x -> parent(x) === F, u) "elements must live in the same module"
    @req length(u) == p "need a $p-tuple of elements"
    return wedge(collect(u); parent=result)
  end
  function my_mult(u::FreeModElem...)
    return my_mult(u)
  end

  function my_decomp(u::FreeModElem)
    @req parent(u) === result "element does not belong to the correct module"
    k = findfirst(==(u), gens(result))
    @req !isnothing(k) "element must be a generator of the module"
    ind = ordered_multi_index(k, p, n)
    e = gens(F)
    return Tuple(e[i] for i in indices(ind))
  end

  mult_map = MapFromFunc(Hecke.TupleParent(Tuple([zero(F) for f in 1:p])), result, my_mult, my_decomp)
  inv_mult_map = MapFromFunc(result, domain(mult_map), my_decomp, my_mult)
  @assert domain(mult_map) === parent(Tuple(zero(F) for i in 1:p)) "something went wrong with the parents"

  # Store the map in the attributes
  set_attribute!(result, :multiplication_map, mult_map)
  set_attribute!(result, :wedge_pure_function, mult_map)
  set_attribute!(result, :wedge_generator_decompose_function, inv_mult_map)
  set_attribute!(result, :is_exterior_power, (F, p))

  cached && (_exterior_powers(F)[p] = (result, mult_map))

  # Set the variable names for printing
  orig_symb = String.(symbols(F))
  new_symb = Symbol[]
  if iszero(p)
    new_symb = [Symbol("1")]
  else
    for ind in OrderedMultiIndexSet(p, n)
      symb_str = orig_symb[ind[1]]
      for i in 2:p
        symb_str = symb_str * (is_unicode_allowed() ? "∧" : "^") * orig_symb[ind[i]]
      end
      push!(new_symb, Symbol(symb_str))
    end
  end
  result.S = new_symb

  set_attribute!(result, :show => show_exterior_product)

  return result, mult_map
end

function symbols(F::FreeMod)
  return F.S
end


########################################################################
# Koszul homology
########################################################################

function koszul_complex(v::FreeModElem; cached::Bool=true)
  F = parent(v)
  R = base_ring(F)
  n = rank(F)
  ext_powers = [exterior_power(F, p, cached=cached)[1] for p in 0:n]
  boundary_maps = [wedge_multiplication_map(ext_powers[i+1], ext_powers[i+2], v) for i in 0:n-1]
  Z = is_graded(F) ? graded_free_module(R, []) : free_module(R, 0)
  pushfirst!(boundary_maps, hom(Z, domain(first(boundary_maps)), elem_type(domain(first(boundary_maps)))[]))
  push!(boundary_maps, hom(codomain(last(boundary_maps)), Z, [zero(Z) for i in 1:ngens(codomain(last(boundary_maps)))]))
  return chain_complex(boundary_maps, seed=-1, check=false)
end

function koszul_complex(v::FreeModElem, M::ModuleFP; cached::Bool=true)
  K = koszul_complex(v, cached=cached)
  KM = tensor_product(K, M)
  return KM
end

function koszul_homology(v::FreeModElem, i::Int; cached::Bool=true)
  F = parent(v)
  n = rank(F)

  # Catch the edge cases
  if i == n # This captures the homological degree zero due to the convention of the chain_complex constructor
    phi = wedge_multiplication_map(exterior_power(F, 0, cached=cached)[1], F, v)
    return kernel(phi)[1]
  end

  if iszero(i) # Homology at the last entry of the complex.
    phi = wedge_multiplication_map(exterior_power(F, n-1)[1], exterior_power(F, n, cached=cached)[1], v)
    return cokernel(phi)[1]
  end

  ext_powers = [exterior_power(F, n-p, cached=cached)[1] for p in i-1:i+1]
  boundary_maps = [wedge_multiplication_map(ext_powers[p+1], ext_powers[p], v) for p in 2:-1:1]
  K = chain_complex(boundary_maps, check=false)
  return homology(K, 1)
end

function koszul_homology(v::FreeModElem, M::ModuleFP, i::Int; cached::Bool=true)
  F = parent(v)
  n = rank(F)

  # Catch the edge cases
  if i == n # This captures the homological degree zero due to the convention of the chain_complex constructor
    phi = wedge_multiplication_map(exterior_power(F, 0, cached=cached)[1], F, v)
    K = chain_complex([phi], check=false)
    KM = tensor_product(K, M)
    return kernel(map(KM, 1))[1]
  end

  if iszero(i) # Homology at the last entry of the complex.
    phi = wedge_multiplication_map(exterior_power(F, n-1)[1], exterior_power(F, n, cached=cached)[1], v)
    K = chain_complex([phi], check=false)
    KM = tensor_product(K, M)
    return cokernel(map(K, 1)) # TODO: cokernel does not seem to return a map by default. Why?
  end

  ext_powers = [exterior_power(F, n-p, cached=cached)[1] for p in i-1:i+1]
  boundary_maps = [wedge_multiplication_map(ext_powers[p+1], ext_powers[p], v) for p in 2:-1:1]
  K = chain_complex(boundary_maps, check=false)
  KM = tensor_product(K, M)
  return homology(KM, 1)
end

function koszul_dual(F::FreeMod; cached::Bool=true)
  success, M, p = _is_exterior_power(F)
  !success && error("module must be an exterior power of some other module")
  return exterior_power(M, rank(M) - p, cached=cached)[1]
end

function koszul_duals(v::Vector{T}; cached::Bool=true) where {T<:FreeModElem}
  isempty(v) && error("list of elements must not be empty")
  allequal(parent, v) || error("parent mismatch")

  F = parent(first(v))
  success, M, p = _is_exterior_power(F)
  n = rank(M)
  success || error("element must be an exterior product")
  k = [findfirst(==(u), gens(F)) for u in v]
  any(isnothing, k) && error("elements must be generators of the module")
  ind = [ordered_multi_index(r, p, n) for r in k]
  comp = [ordered_multi_index([i for i in 1:n if !(i in indices(I))], n) for I in ind]
  lin_ind = linear_index.(comp)
  F_dual = koszul_dual(F, cached=cached)
  results = [F_dual[j] for j in lin_ind]
  for r in 1:length(v)
    sign, _ = _wedge(ind[r], comp[r])
    isone(sign) || (results[r] = -results[r])
  end
  return results 
end

function koszul_dual(v::FreeModElem; cached::Bool=true)
  return first(koszul_duals([v], cached=cached))
end


