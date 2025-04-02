# By default we cache the Kaehler differentials and their 
# exterior powers in the attributes of the ring. This 
# internal function maintains that structure.
function _kaehler_differentials(R::Ring)
  if !has_attribute(R, :kaehler_differentials)
    set_attribute!(R, :kaehler_differentials, Dict{Int, ModuleFP}())
  end
  return get_attribute(R, :kaehler_differentials)::Dict{Int, <:ModuleFP}
end

function kaehler_differentials(R::Union{MPolyRing, MPolyLocRing}; cached::Bool=true)
  if cached && haskey(_kaehler_differentials(R), 1)
    return _kaehler_differentials(R)[1]
  end
  n = ngens(R)
  result = FreeMod(R, n)
  symb = symbols(R)
  result.S = [Symbol(:d, symb[i]) for i in 1:n]
  set_attribute!(result, :show, show_kaehler_differentials)

  cached && (_kaehler_differentials(R)[1] = result)
  set_attribute!(result, :is_kaehler_differential_module, (R, 1))
  return result
end

function kaehler_differentials(R::MPolyDecRing; cached::Bool=true)
  if cached && haskey(_kaehler_differentials(R), 1)
    return _kaehler_differentials(R)[1]
  end
  n = ngens(R)
  result = graded_free_module(R, [1 for i in 1:n])
  symb = symbols(R)
  result.S = [Symbol(:d, symb[i]) for i in 1:n]
  set_attribute!(result, :show, show_kaehler_differentials)

  cached && (_kaehler_differentials(R)[1] = result)
  set_attribute!(result, :is_kaehler_differential_module, (R, 1))
  return result
end

function kaehler_differentials(R::MPolyQuoRing; cached::Bool=true)
  if cached && haskey(_kaehler_differentials(R), 1)
    return _kaehler_differentials(R)[1]
  end
  n = ngens(R)
  P = base_ring(R)
  OmegaP = kaehler_differentials(P)
  f = gens(modulus(R))
  r = length(f)
  Pr = FreeMod(P, r)
  phi = hom(Pr, OmegaP, [exterior_derivative(a, parent=OmegaP) for a in f])
  phi_res, _, _ = change_base_ring(R, phi)
  M = cokernel(phi_res)
  set_attribute!(M, :show, show_kaehler_differentials)

  cached && (_kaehler_differentials(R)[1] = M)
  set_attribute!(M, :is_kaehler_differential_module, (R, 1))
  return M
end

function kaehler_differentials(R::MPolyQuoRing{<:MPolyDecRingElem}; cached::Bool=true)
  if cached && haskey(_kaehler_differentials(R), 1)
    return _kaehler_differentials(R)[1]
  end
  n = ngens(R)
  P = base_ring(R)
  OmegaP = kaehler_differentials(P)
  f = gens(modulus(R))
  r = length(f)
  Pr = graded_free_module(P, r)
  @assert is_graded(OmegaP)
  phi = hom(Pr, OmegaP, [exterior_derivative(a, parent=OmegaP) for a in f])
  phi_res = _change_base_ring_and_preserve_gradings(R, phi)
  @assert is_graded(codomain(phi_res))
  M = cokernel(phi_res)
  @assert is_graded(M)
  set_attribute!(M, :show, show_kaehler_differentials)

  cached && (_kaehler_differentials(R)[1] = M)
  set_attribute!(M, :is_kaehler_differential_module, (R, 1))
  return M
end

function kaehler_differentials(R::MPolyQuoLocRing; cached::Bool=true)
  if cached && haskey(_kaehler_differentials(R), 1)
    return _kaehler_differentials(R)[1]
  end
  n = ngens(R)
  P = localized_ring(R)
  OmegaP = kaehler_differentials(P)
  f = gens(modulus(R))
  r = length(f)
  Pr = FreeMod(P, r)
  phi = hom(Pr, OmegaP, [exterior_derivative(a, parent=OmegaP) for a in f])
  phi_res, _, _ = change_base_ring(R, phi)
  M = cokernel(phi_res)[1]
  set_attribute!(M, :show, show_kaehler_differentials)

  cached && (_kaehler_differentials(R)[1] = M)
  set_attribute!(M, :is_kaehler_differential_module, (R, 1))
  return M
end

@doc raw"""
    kaehler_differentials(R::Ring; cached::Bool=true)
    kaehler_differentials(R::Ring, p::Int; cached::Bool=true)

For `R` a polynomial ring, an affine algebra, or a localization of these 
over a `base_ring` ``𝕜`` this returns the module of Kaehler 
differentials ``Ω¹(R/𝕜)``, resp. its `p`-th exterior power.
"""
function kaehler_differentials(R::Ring, p::Int; cached::Bool=true)
  isone(p) && return kaehler_differentials(R, cached=cached)
  cached && haskey(_kaehler_differentials(R), p) && return _kaehler_differentials(R)[p]
  result = exterior_power(kaehler_differentials(R, cached=cached), p)[1]
  set_attribute!(result, :show, show_kaehler_differentials)
  set_attribute!(result, :is_kaehler_differential_module, (R, p))
  cached && (_kaehler_differentials(R)[p] = result)
  return result
end

@doc raw"""
    is_kaehler_differential_module(M::ModuleFP)

Internal method to check whether a module `M` was created as 
some ``p``-th exterior power of the Kaehler differentials 
``Ω¹(R/𝕜)`` of some ``𝕜``-algebra ``R``. 

Return `(true, R, p)` in the affirmative case and 
`(false, base_ring(M), 0)` otherwise.
"""
function is_kaehler_differential_module(M::ModuleFP)
  has_attribute(M, :is_kaehler_differential_module) || return false, base_ring(M), 0
  R, p = get_attribute(M, :is_kaehler_differential_module)
  return true, R, p
end

@doc raw"""
    de_rham_complex(R::Ring; cached::Bool=true)

Construct the relative de Rham complex of a ``𝕜``-algebra `R`
as a `ComplexOfMorphisms`.
"""
function de_rham_complex(R::Ring; cached::Bool=true)
  n = ngens(R)
  Omega = [kaehler_differentials(R, p, cached=cached) for p in 0:n]
  res_maps = [MapFromFunc(Omega[i], Omega[i+1], w->exterior_derivative(w, parent=Omega[i+1])) for i in 1:n]
  return ComplexOfMorphisms(typeof(Omega[1]), res_maps, typ=:cochain, seed=0, check=false)
end

# printing of kaehler differentials
function show_kaehler_differentials(io::IO, M::ModuleFP)
  success, F, p = _is_exterior_power(M)
  R = base_ring(F)
  if success 
    if is_unicode_allowed() 
      print(io, "Ω^$p($R)")
    else
      print(io, "\\Omega^$p($R)")
    end
  else
    if is_unicode_allowed() 
      print(io, "Ω^1($R)")
    else
      print(io, "\\Omega^1($R)")
    end
  end
end

function show_kaehler_differentials(io::IO, ::MIME"text/html", M::ModuleFP)
  success, F, p = _is_exterior_power(M)
  R = base_ring(F)
  io = IOContext(io, :compact => true)
  if success 
    if is_unicode_allowed() 
      print(io, "Ω^$p($R)")
    else
      print(io, "\\Omega^$p($R)")
    end
  else
    if is_unicode_allowed() 
      print(io, "Ω^1($R)")
    else
      print(io, "\\Omega^1($R)")
    end
  end
end

# Exterior derivatives
@doc raw"""
    exterior_derivative(f::Union{MPolyRingElem, MPolyLocRingElem, MPolyQuoRingElem, MPolyQuoLocRingElem}; 
                        parent::ModuleFP=kaehler_differentials(parent(f)))

Compute the exterior derivative of an element ``f`` of a ``𝕜``-algebra `R`
as an element of the `kaehler_differentials` of `R`.
"""
function exterior_derivative(f::Union{MPolyRingElem, MPolyLocRingElem, MPolyQuoRingElem, MPolyQuoLocRingElem}; 
    parent::ModuleFP=kaehler_differentials(parent(f))
  )
  R = Oscar.parent(f)
  n = ngens(R)
  dx = gens(parent)
  df = sum(derivative(f, i)*dx[i] for i in 1:n; init=zero(parent))
  return df
end

@doc raw"""
    exterior_derivative(w::ModuleFPElem; parent::ModuleFP=...)

Check whether `parent(w)` is an exterior power ``Ωᵖ(R/𝕜)`` of the module of 
Kaehler differentials of some ``𝕜``-algebra `R` and computes its exterior 
derivative in `parent`. If the latter is not specified, it defaults to 
``Ωᵖ⁺¹(R/𝕜)``, the `kaehler_differentials(R, p+1)`.
"""
function exterior_derivative(w::ModuleFPElem; 
    parent::ModuleFP=begin
      success, R, p = is_kaehler_differential_module(parent(w))
      success || error("not a kaehler differential module")
      kaehler_differentials(R, p+1)
    end
  )
  success, R, p = is_kaehler_differential_module(Oscar.parent(w))
  M = Oscar.parent(w)
  result = zero(parent)
  iszero(p) && return exterior_derivative(w[1], parent=parent)
  for (i, a) in coordinates(w)
    da = exterior_derivative(a)
    result = result + wedge(da, M[i], parent=parent)
  end
  return result
end

function symbols(L::MPolyLocRing)
  return symbols(base_ring(L))
end

function change_base_ring(R::Ring, phi::ModuleFPHom{<:ModuleFP, <:ModuleFP, <:Nothing})
  dom_res, dom_map = change_base_ring(R, domain(phi))
  cod_res, cod_map = change_base_ring(R, codomain(phi))
  result = hom(dom_res, cod_res, cod_map.(phi.(gens(domain(phi)))))
  return result, dom_map, cod_map
end

function derivative(a::MPolyQuoRingElem, i::Int)
  return parent(a)(derivative(lift(a), i))
end
