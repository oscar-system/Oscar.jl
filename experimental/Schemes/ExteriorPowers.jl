########################################################################
# Ordered multiindices of the form 0 < i₁ < i₂ < … < iₚ ≤ n.
#
# These provide a coherent way to enumerate the generators of an 
# exterior power ⋀ ᵖ M of a module M given a chosen set of generators 
# of M.
########################################################################
mutable struct OrderedMultiIndex{IntType<:IntegerUnion}
  i::Vector{IntType}
  n::IntType

  function OrderedMultiIndex(i::Vector{T}, n::T) where {T<:IntegerUnion}
    @assert all(k->i[k]<i[k+1], 1:length(i)-1) "indices must be strictly ordered"
    return new{T}(i, n)
  end
end

function ordered_multi_index(i::Vector{T}, n::T) where {T<:IntegerUnion} 
  return OrderedMultiIndex(i, n)
end

# For an ordered multiindex i = (0 < i₁ < i₂ < … < iₚ ≤ n) this returns the vector (i₁,…,iₚ).
indices(a::OrderedMultiIndex) = copy(a.i)

# For an ordered multiindex i = (0 < i₁ < i₂ < … < iₚ ≤ n) this returns n.
bound(a::OrderedMultiIndex) = a.n

# For an ordered multiindex i = (0 < i₁ < i₂ < … < iₚ ≤ n) this returns p.
length(a::OrderedMultiIndex) = length(a.i)

# For an ordered multiindex i = (0 < i₁ < i₂ < … < iₚ ≤ n) this returns iₖ.
getindex(a::OrderedMultiIndex, k::Int) = a.i[k]

index_type(a::OrderedMultiIndex) = index_type(typeof(a))
index_type(::Type{OrderedMultiIndex{T}}) where {T} = T

# Internal function for "multiplication" of ordered multiindices.
# 
# The for i = (0 < i₁ < i₂ < … < iₚ ≤ n) and j = (0 < j₁ < j₂ < … < jᵣ ≤ n)
# the result is a pair `(sign, k)` with `sign` either 0 in case that 
# iₖ = jₗ for some k and l, or ±1 depending on the number of transpositions 
# needed to put (i₁, …, iₚ, j₁, …, jᵣ) into a strictly increasing order.
function _mult(a::OrderedMultiIndex{T}, b::OrderedMultiIndex{T}) where {T}
  @assert bound(a) == bound(b) "multiindices must have the same bounds"
  p = length(a)
  q = length(b)
  result_indices = vcat(indices(a), indices(b))
  sign = 1

  # bubble sort result_indices and keep track of the sign
  for k in p:-1:1
    l = k
    while l < p + q && result_indices[l] >= result_indices[l+1]
      # in case of a double index return zero
      result_indices[l] == result_indices[l+1] && return 0, result_indices
      c = result_indices[l+1]
      result_indices[l+1] = result_indices[l]
      result_indices[l] = c
      sign = -sign
      l = l+1
    end
  end
  return sign, result_indices
end

# For two ordered multiindices i = (0 < i₁ < i₂ < … < iₚ ≤ n) 
# and j = (0 < j₁ < j₂ < … < jᵣ ≤ n) this returns a pair `(sign, k)` 
# with `sign` either 0 in case that iₖ = jₗ for some k and l, 
# or ±1 depending on the number of transpositions needed to put 
# (i₁, …, iₚ, j₁, …, jᵣ) into a strictly increasing order.
function wedge(a::OrderedMultiIndex{T}, b::OrderedMultiIndex{T}) where {T}
  sign, ind = _mult(a, b)
  iszero(sign) && return sign, OrderedMultiIndex([i for i in 1:length(a) + length(b)], bound(a))
  return sign, OrderedMultiIndex(ind, bound(a))
end

function wedge(a::Vector{T}) where {T <: OrderedMultiIndex}
  isempty(a) && error("list must not be empty")
  isone(length(a)) && return 1, first(a)
  b = first(a)
  rem = a[2:end]
  sign, ind = wedge(rem)
  new_sign, new_ind = wedge(b, ind)
  return sign*new_sign, new_ind
end

function ==(a::OrderedMultiIndex{T}, b::OrderedMultiIndex{T}) where {T}
  bound(a) == bound(b) || return false
  return indices(a) == indices(b)
end

########################################################################
# A data type to facilitate iteration over all ordered multiindices 
# of the form 0 < i₁ < i₂ < … < iₚ ≤ n for fixed 0 < p ≤ n.
#
# Example:
#
#   for i in OrderedMultiIndexSet(3, 5)
#      # do something with i = (0 < i₁ < i₂ < i₃ ≤ 5).
#   end
########################################################################
mutable struct OrderedMultiIndexSet
  n::Int
  p::Int

  function OrderedMultiIndexSet(p::Int, n::Int)
    @assert p <= n
    return new(n, p)
  end
end

bound(I::OrderedMultiIndexSet) = I.n
index_length(I::OrderedMultiIndexSet) = I.p

Base.eltype(I::OrderedMultiIndexSet) = OrderedMultiIndex
Base.length(I::OrderedMultiIndexSet) = binomial(bound(I), index_length(I))

function Base.iterate(I::OrderedMultiIndexSet)
  ind = OrderedMultiIndex([i for i in 1:index_length(I)], bound(I))
  return ind, ind
end

function Base.iterate(I::OrderedMultiIndexSet, state::OrderedMultiIndex)
  bound(I) == bound(state) || error("index not compatible with set")
  ind = indices(state)
  l = length(state)
  while l > 0 && ind[l] == bound(I) - length(state) + l
    l = l - 1
  end
  iszero(l) && return nothing
  ind[l] = ind[l] + 1
  l = l + 1
  while l <= length(state)
    ind[l] = ind[l-1] + 1
    l = l + 1
  end
  result = OrderedMultiIndex(ind, bound(I))
  return result, result
end

function Base.show(io::IO, ind::OrderedMultiIndex)
  i = indices(ind)
  print(io, "0 ")
  for i in indices(ind)
    print(io, "< $i ")
  end
  print(io, "<= $(bound(ind))")
end

# For an ordered multiindex i = (0 < i₁ < i₂ < … < iₚ ≤ n) this 
# returns the number k so that i appears at the k-th spot in the 
# enumeration of all ordered multiindices for this pair 0 < p ≤ n.
function linear_index(ind::OrderedMultiIndex)
  n = bound(ind)
  p = length(ind)
  iszero(p) && return 1
  isone(p) && return ind[1]
  i = indices(ind)
  return binomial(n, p) - binomial(n - first(i) + 1, p) + linear_index(OrderedMultiIndex(i[2:end].-first(i), n-first(i)))
end

# For a pair 0 < p ≤ n return the k-th ordered multiindex in the 
# enumeration of all ordered multiindices (0 < i₁ < i₂ < … < iₚ ≤ n).
function ordered_multi_index(k::Int, p::Int, n::Int)
  (k < 1 || k > binomial(n, p)) && error("index out of range")
  isone(p) && return OrderedMultiIndex([k], n)
  i1 = 1
  bin = binomial(n, p)
  while !(bin - binomial(n - i1, p) > k - 1)
    i1 = i1 + 1
  end
  prev_res = ordered_multi_index(k - bin + binomial(n - i1 + 1, p), p-1, n - i1)
  return OrderedMultiIndex(pushfirst!(indices(prev_res).+i1, i1), n)
end

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
  img_gens = [sum(v[l] * wedge(f, e) for (l, f) in enumerate(gens(H)); init=zero(G)) for e in gens(F)]
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
  for i in OrderedMultiIndexSet(p, n)
    for j in OrderedMultiIndexSet(q, n)
      sign, k = wedge(i, j)
      iszero(sign) && continue
      result = result + sign * u[linear_index(i)] * v[linear_index(j)] * parent(result)[linear_index(k)]
    end
  end
  return result
end

function wedge(u::Vector{T}) where {T<:FreeModElem}
  isone(length(u)) && return first(u)
  return wedge(first(u), wedge(u[2:end]))
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
