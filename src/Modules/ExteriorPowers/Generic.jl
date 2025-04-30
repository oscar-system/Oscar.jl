# We need to cache eventually created exterior powers.
@attr Dict{Int, Tuple{typeof(F), MapFromFunc}} function _exterior_powers(F::SparseFPModule)
  return Dict{Int, Tuple{typeof(F), MapFromFunc}}()
end

# User facing method to ask whether F = ⋀ ᵖ M for some M.
# This returns a triple `(true, M, p)` in the affirmative case
# and `(false, F, 0)` otherwise.
function _is_exterior_power(M::SparseFPModule)
  if has_attribute(M, :is_exterior_power)
    MM, p = get_attribute(M, :is_exterior_power)::Tuple{typeof(M), Int}
    return (true, MM, p)
  end
  return (false, M, 0)
end

# Printing of exterior powers
function show_exterior_product(io::IO, M::SparseFPModule)
  success, F, p = _is_exterior_power(M)
  success || error("module is not an exterior power")
  if is_unicode_allowed()
    print(io, "⋀^$p($F)")
  else
    print(io, "$(ordinal_number_string(p)) exterior power of $F")
  end
end

function show_exterior_product(io::IO, ::MIME"text/html", M::SparseFPModule)
  success, F, p = _is_exterior_power(M)
  success || error("module is not an exterior power")
  io = IOContext(io, :compact => true)
  if is_unicode_allowed()
    print(io, "⋀^$p($F)")
  else
    print(io, "$(ordinal_number_string(p)) exterior power of $F")
  end
end

function multiplication_map(M::SparseFPModule)
  has_attribute(M, :multiplication_map) || error("module is not an exterior power")
  return get_attribute(M, :multiplication_map)::Map
end

function wedge_pure_function(M::SparseFPModule)
  has_attribute(M, :wedge_pure_function) || error("module is not an exterior power")
  return get_attribute(M, :wedge_pure_function)::Map
end

function wedge_generator_decompose_function(M::SparseFPModule)
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
function wedge_multiplication_map(F::SparseFPModule, G::SparseFPModule, v::SparseFPModuleElem)
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
function wedge(u::SparseFPModuleElem, v::SparseFPModuleElem; 
    parent::SparseFPModule=begin
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
    ind_i = ordered_multi_index(i, p, n)
    for (j, b) in coordinates(v)
      ind_j = ordered_multi_index(j, q, n)
      sign, k = _wedge(ind_i, ind_j)
      iszero(sign) && continue
      result = result + sign * a * b * parent[linear_index(k)]
    end
  end
  return result
end

function wedge(u::Vector{T};
    parent::SparseFPModule=begin
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
  ) where {T<:SparseFPModuleElem}
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

function induced_map_on_exterior_power(phi::FreeModuleHom{<:FreeMod, <:FreeMod, Nothing}, p::Int;
    domain::FreeMod=exterior_power(Oscar.domain(phi), p)[1],
    codomain::FreeMod=exterior_power(Oscar.codomain(phi), p)[1]
  )
  F = Oscar.domain(phi)
  m = rank(F)
  G = Oscar.codomain(phi)
  n = rank(G)

  is_zero(p) && return hom(domain, codomain, gens(codomain); check=false) # Isomorphism of R^1

  imgs = phi.(gens(F))
  img_gens = [wedge(imgs[indices(ind)], parent=codomain) for ind in OrderedMultiIndexSet(p, m)]
  return hom(domain, codomain, img_gens; check=false)
end

# The induced map on exterior powers
function hom(M::FreeMod, N::FreeMod, phi::FreeModuleHom)
  success, F, p = _is_exterior_power(M)
  @req success "module is not an exterior power"
  success, FF, q = _is_exterior_power(N)
  @req success "module is not an exterior power"
  @req F === domain(phi) "map not compatible"
  @req FF === codomain(phi) "map not compatible"
  @req p == q "exponents must agree"
  return induced_map_on_exterior_power(phi, p; domain=M, codomain=N)
end


