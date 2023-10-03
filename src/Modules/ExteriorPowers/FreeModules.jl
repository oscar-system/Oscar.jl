########################################################################
# Exterior powers of free modules
#
# For F = Rⁿ we provide methods to create ⋀ ᵖF for arbitrary 0 ≤ p ≤ n.
# These modules are cached in F and know that they are an exterior 
# power of F. This allows us to implement the wedge product of their 
# elements. 
########################################################################

# We need to cache eventually created exterior powers.
@attr Dict{Int, <:FreeMod} function _exterior_powers(F::FreeMod) 
  return Dict{Int, typeof(F)}()
end

# User facing constructor for ⋀ ᵖ F.
function exterior_power(F::FreeMod, p::Int)
  (p < 0 || p > rank(F)) && error("index out of bounds")
  powers = _exterior_powers(F)
  haskey(powers, p) && return powers[p]::typeof(F)

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

  set_attribute!(result, :is_exterior_power, (F, p))
  powers[p] = result
  return result
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
    Fwedge1 = exterior_power(F, 1)
    id = hom(F, Fwedge1, gens(Fwedge1))
    tmp = wedge_multiplication_map(Fwedge1, G, v)
    return compose(id, tmp)
  end

  success, orig_mod_2, q = is_exterior_power(G)
  if !success
    Gwedge1 = exterior_power(G, 1)
    id = hom(Gwedge1, G, gens(G))
    tmp = wedge_multiplication_map(F, Gwedge1, v)
    return compose(tmp, id)
  end

  orig_mod === orig_mod_2 || error("modules must be exterior powers of the same module")
  H = parent(v)

  # In case v comes from the original module, convert.
  if H === orig_mod
    M = exterior_power(orig_mod, 1)
    w = M(coordinates(v))
    return wedge_multiplication_map(F, G, w)
  end

  success, orig_mod_2, r = is_exterior_power(H)
  success || error("element is not an exterior product")
  orig_mod_2 === orig_mod || error("element is not an exterior product for the correct module")
  p + r == q || error("powers are incompatible")
  
  # map the generators
  img_gens = [wedge(v, e) for e in gens(F)]
  return hom(F, G, img_gens)
end

# The wedge product of two or more elements.
function wedge(u::FreeModElem, v::FreeModElem)
  success1, F1, p = is_exterior_power(parent(u))
  if !success1
    F = parent(u) 
    Fwedge1 = exterior_power(F1, 1)
    return wedge(Fwedge1(coordinates(u)), v)
  end

  success2, F2, q = is_exterior_power(parent(v))
  if !success2
    F = parent(v) 
    Fwedge1 = exterior_power(F1, 1)
    return wedge(u, Fwedge1(coordinates(v)))
  end

  F1 === F2 || error("modules are not exterior powers of the same original module")
  n = rank(F1)

  result = zero(exterior_power(F1, p + q))
  for (i, ind_i) in enumerate(OrderedMultiIndexSet(p, n))
    a = u[i]
    iszero(a) && continue
    for (j, ind_j) in enumerate(OrderedMultiIndexSet(q, n))
      b = v[j]
      iszero(b) && continue
      sign, k = _wedge(ind_i, ind_j)
      iszero(sign) && continue
      result = result + sign * a * b * parent(result)[linear_index(k)]
    end
  end
  return result
end

function wedge(u::Vector{T}) where {T<:FreeModElem}
  isempty(u) && error("list must not be empty")
  isone(length(u)) && return first(u)
  k = div(length(u), 2)
  return wedge(wedge(u[1:k]), wedge(u[k+1:end]))
end


########################################################################
# Koszul homology
########################################################################

function koszul_complex(v::FreeModElem)
  F = parent(v)
  n = rank(F)
  ext_powers = [exterior_power(F, p) for p in 0:n]
  boundary_maps = [wedge_multiplication_map(ext_powers[i+1], ext_powers[i+2], v) for i in 0:n-1]
  return chain_complex(boundary_maps)
end

function koszul_complex(v::FreeModElem, M::ModuleFP)
  K = koszul_complex(v)
  KM = tensor_product(K, M)
  return KM
end

function koszul_homology(v::FreeModElem, i::Int)
  F = parent(v)
  n = rank(F)

  # Catch the edge cases
  if i == n # This captures the homological degree zero due to the convention of the chain_complex constructor
    phi = wedge_multiplication_map(exterior_power(F, 0), F, v)
    return kernel(phi)[1]
  end

  if iszero(i) # Homology at the last entry of the complex.
    phi = wedge_multiplication_map(exterior_power(F, n-1), exterior_power(F, n), v)
    return cokernel(phi)[1]
  end

  ext_powers = [exterior_power(F, p) for p in i-1:i+1]
  boundary_maps = [wedge_multiplication_map(ext_powers[p], ext_powers[p+1], v) for p in 1:2]
  K = chain_complex(boundary_maps)
  return homology(K, 1)
end

function koszul_homology(v::FreeModElem, M::ModuleFP, i::Int)
  F = parent(v)
  n = rank(F)

  # Catch the edge cases
  if i == n # This captures the homological degree zero due to the convention of the chain_complex constructor
    phi = wedge_multiplication_map(exterior_power(F, 0), F, v)
    K = chain_complex([phi])
    KM = tensor_product(K, M)
    return kernel(map(KM, 1))[1]
  end

  if iszero(i) # Homology at the last entry of the complex.
    phi = wedge_multiplication_map(exterior_power(F, n-1), exterior_power(F, n), v)
    K = chain_complex([phi])
    KM = tensor_product(K, M)
    return cokernel(map(K, 1)) # TODO: cokernel does not seem to return a map by default. Why?
  end

  ext_powers = [exterior_power(F, p) for p in i-1:i+1]
  boundary_maps = [wedge_multiplication_map(ext_powers[p], ext_powers[p+1], v) for p in 1:2]
  K = chain_complex(boundary_maps)
  return homology(K, 1)
end

function koszul_dual(F::FreeMod)
  success, M, p = is_exterior_power(F)
  !success && error("module must be an exterior power of some other module")
  return exterior_power(M, rank(M) - p)
end

function koszul_duals(v::Vector{T}) where {T<:FreeModElem}
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
  F_dual = koszul_dual(F)
  results = [F_dual[j] for j in lin_ind]
  for r in 1:length(v)
    sign, _ = _wedge(ind[r], comp[r])
    isone(sign) || (results[r] = -results[r])
  end
  return results 
end

function koszul_dual(v::FreeModElem)
  return first(koszul_duals([v]))
end

function induced_map_on_exterior_power(phi::FreeModuleHom{<:FreeMod, <:FreeMod, Nothing}, p::Int)
  F = domain(phi)
  m = rank(F)
  G = codomain(phi)
  n = rank(G)

  Fp = exterior_power(F, p)
  Gp = exterior_power(G, p)

  imgs = phi.(gens(F))
  img_gens = [wedge(imgs[indices(ind)]) for ind in OrderedMultiIndexSet(p, m)]
  return hom(Fp, Gp, img_gens)
end

