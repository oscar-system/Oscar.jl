# We need to cache eventually created exterior powers.
@attr Dict{Int, Tuple{typeof(F), MapFromFunc}} function _exterior_powers(F::ModuleFP)
  return Dict{Int, Tuple{typeof(F), MapFromFunc}}()
end

# User facing method to ask whether F = ⋀ ᵖ M for some M.
# This returns a triple `(true, M, p)` in the affirmative case
# and `(false, F, 0)` otherwise.
function _is_exterior_power(M::ModuleFP)
  if has_attribute(M, :is_exterior_power)
    MM, p = get_attribute(M, :is_exterior_power)::Tuple{typeof(M), Int}
    return (true, MM, p)
  end
  return (false, M, 0)
end

# Printing of exterior powers
function show_exterior_product(io::IO, M::ModuleFP)
  success, F, p = _is_exterior_power(M)
  success || error("module is not an exterior power")
  if is_unicode_allowed()
    print(io, "⋀^$p($F)")
  else
    print(io, "$(ordinal_number_string(p)) exterior power of $F")
  end
end

function show_exterior_product(io::IO, ::MIME"text/html", M::ModuleFP)
  success, F, p = _is_exterior_power(M)
  success || error("module is not an exterior power")
  io = IOContext(io, :compact => true)
  if is_unicode_allowed()
    print(io, "⋀^$p($F)")
  else
    print(io, "$(ordinal_number_string(p)) exterior power of $F")
  end
end

function multiplication_map(M::ModuleFP)
  has_attribute(M, :multiplication_map) || error("module is not an exterior power")
  return get_attribute(M, :multiplication_map)::Map
end

function wedge_pure_function(M::ModuleFP)
  has_attribute(M, :wedge_pure_function) || error("module is not an exterior power")
  return get_attribute(M, :wedge_pure_function)::Map
end

function wedge_generator_decompose_function(M::ModuleFP)
  has_attribute(M, :wedge_generator_decompose_function) || error("module is not an exterior power")
  return get_attribute(M, :wedge_generator_decompose_function)::Map
end

# Given two exterior powers F = ⋀ ᵖM and G = ⋀ ʳM and an element
# v ∈ ⋀ ʳ⁻ᵖ M this constructs the module homomorphism associated
# to
#
#   v ∧ - : F → G,  u ↦ v ∧ u.
#
# We also allow v ∈ M considered as ⋀ ¹M and the same holds in
# the cases p = 1 and r = 1.
function wedge_multiplication_map(F::ModuleFP, G::ModuleFP, v::ModuleFPElem)
  success, orig_mod, p = _is_exterior_power(F)
  if !success
    Fwedge1, _ = exterior_power(F, 1)
    id = hom(F, Fwedge1, gens(Fwedge1); check=false)
    tmp = wedge_multiplication_map(Fwedge1, G, v)
    return compose(id, tmp)
  end

  success, orig_mod_2, q = _is_exterior_power(G)
  if !success
    Gwedge1, _ = exterior_power(G, 1)
    id = hom(Gwedge1, G, gens(G); check=false)
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

  success, orig_mod_2, r = _is_exterior_power(H)
  success || error("element is not an exterior product")
  orig_mod_2 === orig_mod || error("element is not an exterior product for the correct module")
  p + r == q || error("powers are incompatible")

  # map the generators
  img_gens = [wedge(v, e, parent=G) for e in gens(F)]
  return hom(F, G, img_gens; check=false)
end

# The wedge product of two or more elements.
function wedge(u::ModuleFPElem, v::ModuleFPElem;
    parent::ModuleFP=begin
      success, F, p = _is_exterior_power(Oscar.parent(u))
      if !success
        F = Oscar.parent(u)
        p = 1
      end
      success, _, q = _is_exterior_power(Oscar.parent(v))
      !success && (q = 1)
      exterior_power(F, p + q)[1]
    end
  )
  success1, F1, p = _is_exterior_power(Oscar.parent(u))
  if !success1
    F = Oscar.parent(u)
    Fwedge1, _ = exterior_power(F1, 1)
    return wedge(Fwedge1(coordinates(u)), v, parent=parent)
  end

  success2, F2, q = _is_exterior_power(Oscar.parent(v))
  if !success2
    F = Oscar.parent(v)
    Fwedge1, _ = exterior_power(F1, 1)
    return wedge(u, Fwedge1(coordinates(v)), parent=parent)
  end

  F1 === F2 || error("modules are not exterior powers of the same original module")
  n = ngens(F1)

  result = zero(parent)
  for (i, a) in coordinates(u)
    ind_i = combination(n, p, i)
    for (j, b) in coordinates(v)
      ind_j = combination(n, q, j)
      sign, k = _wedge(ind_i, ind_j)
      iszero(sign) && continue
      result = result + sign * a * b * parent[Oscar.linear_index(k, n)]
    end
  end
  return result
end

function wedge(u::Vector{T};
    parent::ModuleFP=begin
      r = 0
      isempty(u) && error("list must not be empty")
      F = Oscar.parent(first(u)) # initialize variable
      for v in u
        success, F, p = _is_exterior_power(Oscar.parent(v))
        if !success
          F = Oscar.parent(v)
          p = 1
        end
        r = r + p
      end
      exterior_power(F, r)[1]
    end
  ) where {T<:ModuleFPElem}
  isempty(u) && error("list must not be empty")
  if isone(length(u))
    Oscar.parent(first(u)) === parent && return first(u)
    return parent(coordinates(first(u)))
  end
  k = div(length(u), 2)
  result = wedge(wedge(u[1:k]), wedge(u[k+1:end]), parent=parent)
  @assert Oscar.parent(result) === parent
  return result
end

function exterior_power(phi::ModuleFPHom{<:ModuleFP, <:ModuleFP, Nothing}, p::Int; 
    cached::Bool=true,
    domain::ModuleFP=exterior_power(domain(phi), p; cached)[1],
    codomain::ModuleFP=exterior_power(codomain(phi), p; cached)[1]
  )
  if ngens(Oscar.domain(phi)) == ngens(Oscar.codomain(phi)) == p
    return hom(domain, codomain, [det(matrix(phi))*codomain[1]]; check=false)
  end

  A = matrix(phi)
  img_gens = elem_type(codomain)[]
  m = ngens(Oscar.domain(phi))
  n = ngens(Oscar.codomain(phi))
  R = base_ring(codomain)
  for (i, I) in enumerate(combinations(m, p))
    dat = Tuple{Int, elem_type(R)}[]
    for (j, J) in enumerate(combinations(n, p))
      c = det(A[data(I), data(J)]) 
      is_zero(c) && continue
      push!(dat, (j, c))
    end
    push!(img_gens, codomain(sparse_row(R, dat)))
  end
  return hom(domain, codomain, img_gens)
end

# with coefficient map
function exterior_power(phi::ModuleFPHom{<:ModuleFP, <:ModuleFP}, p::Int; 
    cached::Bool=true,
    domain::ModuleFP=exterior_power(domain(phi), p; cached)[1],
    codomain::ModuleFP=exterior_power(codomain(phi), p; cached)[1]
  )
  if ngens(Oscar.domain(phi)) == ngens(Oscar.codomain(phi)) == p
    return hom(domain, codomain, [det(matrix(phi))*codomain[1]], base_ring_map(phi))
  end

  A = matrix(phi)
  img_gens = elem_type(codomain)[]
  m = ngens(Oscar.domain(phi))
  n = ngens(Oscar.codomain(phi))
  R = base_ring(codomain)
  for (i, I) in enumerate(combinations(m, p))
    dat = Tuple{Int, elem_type(R)}[]
    for (j, J) in enumerate(combinations(n, p))
      c = det(A[data(I), data(J)]) 
      is_zero(c) && continue
      push!(dat, (j, c))
    end
    push!(img_gens, codomain(sparse_row(R, dat)))
  end
  return hom(domain, codomain, img_gens, base_ring_map(phi))
end

