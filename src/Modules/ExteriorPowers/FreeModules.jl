########################################################################
# Exterior powers of free modules
#
# For F = Rⁿ we provide methods to create ⋀ ᵖF for arbitrary 0 ≤ p ≤ n.
# These modules are cached in F and know that they are an exterior 
# power of F. This allows us to implement the wedge product of their 
# elements. 
########################################################################

# We need to cache eventually created exterior powers.
@attr Dict{Int, Tuple{T, <:Map}} function _exterior_powers(F::T) where {T<:FreeMod}
  return Dict{Int, Tuple{typeof(F), Map}}()
end

# User facing constructor for ⋀ ᵖ F.
function exterior_power(F::FreeMod, p::Int; cached::Bool=true)
  (p < 0 || p > rank(F)) && error("index out of bounds")

  if cached
    powers = _exterior_powers(F)
    haskey(powers, p) && return powers[p]::Tuple{typeof(F), <:Map}
  end

  R = base_ring(F)
  n = rank(F)
  result = FreeMod(R, binomial(n, p))

  # In case F was graded, we have to take an extra detour. 
  if is_graded(F)
    G = grading_group(F)
    weights = elem_type(G)[]
    for ind in OrderedMultiIndexSet(p, n)
      push!(weights, sum(degree(F[i]) for i in indices(ind); init=zero(G)))
    end
    result = grade(result, weights)
  end

  # Create the multiplication map
  function my_mult(u::FreeModElem...)
    isempty(u) && return result[1] # only the case p=0
    @assert all(x->parent(x)===F, u) "elements must live in the same module"
    @assert length(u) == p "need a $p-tuple of elements"
    return wedge(collect(u), parent=result)
  end
  function my_mult(u::Tuple)
    return my_mult(u...)
  end

  function my_decomp(u::FreeModElem)
    parent(u) === result || error("element does not belong to the correct module")
    k = findfirst(x->x==u, gens(result)) 
    k === nothing && error("element must be a generator of the module")
    ind = ordered_multi_index(k, p, n)
    e = gens(F)
    return Tuple(e[i] for i in indices(ind))
  end

  mult_map = MapFromFunc(Hecke.TupleParent(Tuple([zero(F) for f in 1:p])), result, my_mult, my_decomp)
  @assert domain(mult_map) === parent(Tuple(zero(F) for i in 1:p)) "something went wrong with the parents"

  set_attribute!(result, :multiplication_map, mult_map)
  set_attribute!(result, :is_exterior_power, (F, p))

  cached && (_exterior_powers(F)[p] = (result, mult_map))

  return result, mult_map
end

function multiplication_map(M::FreeMod)
  has_attribute(M, :multiplication_map) || error("module is not an exterior power")
  return get_attribute(M, :multiplication_map)::Map
end

# User facing method to ask whether F = ⋀ ᵖ M for some M.
# This returns a triple `(true, M, p)` in the affirmative case
# and `(false, F, 0)` otherwise.
function is_exterior_power(F::FreeMod)
  !has_attribute(F, :is_exterior_power) && return false, F, 0
  orig_mod, p = get_attribute(F, :is_exterior_power)::Tuple{typeof(F), Int}
  return true, orig_mod, p
end

# Given two exterior powers F = ⋀ ᵖM and G = ⋀ ʳM and an element 
# v ∈ ⋀ ʳ⁻ᵖ M this constructs the module homomorphism associated 
# to 
#
#   v ∧ - : F → G,  u ↦ v ∧ u.
#
# We also allow v ∈ M considered as ⋀ ¹M and the same holds in 
# the cases p = 1 and r = 1.
function wedge_multiplication_map(F::FreeMod, G::FreeMod, v::FreeModElem)
  success, orig_mod, p = is_exterior_power(F)
  if !success 
    Fwedge1, _ = exterior_power(F, 1)
    id = hom(F, Fwedge1, gens(Fwedge1))
    tmp = wedge_multiplication_map(Fwedge1, G, v)
    return compose(id, tmp)
  end

  success, orig_mod_2, q = is_exterior_power(G)
  if !success
    Gwedge1, _ = exterior_power(G, 1)
    id = hom(Gwedge1, G, gens(G))
    tmp = wedge_multiplication_map(F, Gwedge1, v)
    return compose(tmp, id)
  end

  orig_mod === orig_mod_2 || error("modules must be exterior powers of the same module")
  H = parent(v)

  # In case v comes from the original module, convert.
  if H === orig_mod
    M, _ = exterior_power(orig_mod, 1)
    w = M(coordinates(v))
    return wedge_multiplication_map(F, G, w)
  end

  success, orig_mod_2, r = is_exterior_power(H)
  success || error("element is not an exterior product")
  orig_mod_2 === orig_mod || error("element is not an exterior product for the correct module")
  p + r == q || error("powers are incompatible")
  
  # map the generators
  img_gens = [wedge(v, e, parent=G) for e in gens(F)]
  return hom(F, G, img_gens)
end

# The wedge product of two or more elements.
function wedge(u::FreeModElem, v::FreeModElem; 
    parent::FreeMod=begin
      success, F, p = is_exterior_power(Oscar.parent(u))
      if !success 
        F = Oscar.parent(u)
        p = 1
      end
      success, _, q = is_exterior_power(Oscar.parent(v))
      !success && (q = 1)
      exterior_power(F, p + q)[1]
    end
  )
  success1, F1, p = is_exterior_power(Oscar.parent(u))
  if !success1
    F = Oscar.parent(u) 
    Fwedge1, _ = exterior_power(F1, 1)
    return wedge(Fwedge1(coordinates(u)), v, parent=parent)
  end

  success2, F2, q = is_exterior_power(Oscar.parent(v))
  if !success2
    F = Oscar.parent(v) 
    Fwedge1, _ = exterior_power(F1, 1)
    return wedge(u, Fwedge1(coordinates(v)), parent=parent)
  end

  F1 === F2 || error("modules are not exterior powers of the same original module")
  n = rank(F1)

  result = zero(parent)
  a_ind = [k for (k, _) in coordinates(u)]
  for i in a_ind
    ind_i = ordered_multi_index(i, p, n)
    a = u[i]
    b_ind = [k for (k, _) in coordinates(v)]
    for j in b_ind
      ind_j = ordered_multi_index(j, q, n)
      sign, k = _wedge(ind_i, ind_j)
      iszero(sign) && continue
      result = result + sign * a * v[j] * parent[linear_index(k)]
    end
  end
  return result
end

function wedge(u::Vector{T};
    parent::FreeMod=begin
      r = 0
      isempty(u) && error("list must not be empty")
      F = Oscar.parent(first(u)) # initialize variable
      for v in u
        success, F, p = is_exterior_power(Oscar.parent(v))
        if !success 
          F = Oscar.parent(v)
          p = 1
        end
        r = r + p
      end
      exterior_power(F, r)[1]
    end
  ) where {T<:FreeModElem}
  isempty(u) && error("list must not be empty")
  isone(length(u)) && return first(u)
  k = div(length(u), 2)
  result = wedge(wedge(u[1:k]), wedge(u[k+1:end]), parent=parent)
  @assert Oscar.parent(result) === parent
  return result
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
  return chain_complex(boundary_maps, seed=-1)
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

  ext_powers = [exterior_power(F, p, cached=cached)[1] for p in i-1:i+1]
  boundary_maps = [wedge_multiplication_map(ext_powers[p], ext_powers[p+1], v) for p in 1:2]
  K = chain_complex(boundary_maps)
  return homology(K, 1)
end

function koszul_homology(v::FreeModElem, M::ModuleFP, i::Int; cached::Bool=true)
  F = parent(v)
  n = rank(F)

  # Catch the edge cases
  if i == n # This captures the homological degree zero due to the convention of the chain_complex constructor
    phi = wedge_multiplication_map(exterior_power(F, 0, cached=cached)[1], F, v)
    K = chain_complex([phi])
    KM = tensor_product(K, M)
    return kernel(map(KM, 1))[1]
  end

  if iszero(i) # Homology at the last entry of the complex.
    phi = wedge_multiplication_map(exterior_power(F, n-1)[1], exterior_power(F, n, cached=cached)[1], v)
    K = chain_complex([phi])
    KM = tensor_product(K, M)
    return cokernel(map(K, 1)) # TODO: cokernel does not seem to return a map by default. Why?
  end

  ext_powers = [exterior_power(F, p, cached=cached)[1] for p in i-1:i+1]
  boundary_maps = [wedge_multiplication_map(ext_powers[p], ext_powers[p+1], v) for p in 1:2]
  K = chain_complex(boundary_maps)
  return homology(K, 1)
end

function koszul_dual(F::FreeMod; cached::Bool=true)
  success, M, p = is_exterior_power(F)
  !success && error("module must be an exterior power of some other module")
  return exterior_power(M, rank(M) - p, cached=cached)[1]
end

function koszul_duals(v::Vector{T}; cached::Bool=true) where {T<:FreeModElem}
  isempty(v) && error("list of elements must not be empty")
  all(u->parent(u) === parent(first(v)), v[2:end]) || error("parent mismatch")

  F = parent(first(v))
  success, M, p = is_exterior_power(F)
  n = rank(M)
  success || error("element must be an exterior product")
  k = [findfirst(x->x==u, gens(F)) for u in v]
  any(x->x===nothing, k) && error("elements must be generators of the module")
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

function induced_map_on_exterior_power(phi::FreeModuleHom{<:FreeMod, <:FreeMod, Nothing}, p::Int;
    domain::FreeMod=exterior_power(Oscar.domain(phi), p)[1],
    codomain::FreeMod=exterior_power(Oscar.codomain(phi), p)[1]
  )
  F = Oscar.domain(phi)
  m = rank(F)
  G = Oscar.codomain(phi)
  n = rank(G)

  imgs = phi.(gens(F))
  img_gens = [wedge(imgs[indices(ind)], parent=codomain) for ind in OrderedMultiIndexSet(p, m)]
  return hom(domain, codomain, img_gens)
end

