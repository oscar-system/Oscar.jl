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

########################################################################
# relative Kaehler differentials
########################################################################

# Suppose we are given a family 
#   
#   X₀ ↪ X ↪ 𝔸ⁿ×𝔸¹
#   ↓    ↓     ↓π
#  {0} ↪ 𝔸¹ ≅  𝔸¹
#
# over a line with parameter `t`. Then X will be modeled as Spec(R) 
# for R = 𝕜[x₁,…,xₙ,t]/I. In many cases one is interested in the relative 
# Kaehler diffentials Ω¹_{R/𝕜[t]} over the base (rather than over 𝕜). 
#
# Conceptually, one would probably require to start with the projection 
# map π and construct the relative Kaehler differentials from there. 
# In practice, however, one will basically always work in a polynomial 
# ring for the ambient space, in which some variables are labeled as 
# "parameters", i.e. coming from the base space. For implementations 
# this seems to be the more practical entry point, at least for the 
# first internal functions.
#
# The plan is to start using this in research while not exporting it. 
# When we have gathered some more experience and further opinions on the 
# topic might have matured, we will eventually wrap it up in a user-facing 
# set of methods which can then be exported.

@doc raw"""
    relative_kaehler_differentials(
        R::MPolyAnyRing; base_variables::Vector{T}=elem_type(R)[]
      ) where {T <: RingElem}

Given a polynomial ring `R` (or a quotient/localization/...) and a distinguished set 
of generators (variables) `base_variables` of `R`, compute the module of relative 
Kaehler differentials ``Ω¹_{R/B}`` over the subring ``B ⊂ R`` generated by the 
`base_variables`.
"""
function relative_kaehler_differentials(
    R::MPolyAnyRing; 
    base_variables::Vector{T}=elem_type(R)[] # the parameter variables in the above sense
  ) where {T <: RingElem}
  @assert all(x in gens(R) for x in base_variables) "given polynomials must be generators of the ring"
  @assert _allunique(base_variables)
  return _relative_kaehler_differentials(R, sort!([findfirst(==(y), gens(R)) for y in base_variables]))
end

function _relative_kaehler_differentials(
    R::Union{MPolyRing, MPolyLocRing},
    ind::Vector{Int} # indices of the parameter variables in increasing order
  )
  symb = symbols(R)
  symb = [Symbol(:d, symb[i]) for i in 1:ngens(R) if !(i in ind)]
  result = FreeMod(R, symb)
  #set_attribute!(result, :show, show_kaehler_differentials)
  #set_attribute!(result, :is_kaehler_differential_module, (R, 1))
  return result
end

function _relative_kaehler_differentials(
    R::Union{MPolyQuoRing, MPolyQuoLocRing},
    ind::Vector{Int}
  )
  n = ngens(R)
  m = length(ind)
  ind_comp = [i for i in 1:n if !(i in ind)]
  P = base_ring(R)
  OmegaP = _relative_kaehler_differentials(P, ind)
  f = gens(saturated_ideal(modulus(R)))
  r = length(f)
  Pr = FreeMod(P, r)
  img_gens = [sum(derivative(f, i)*OmegaP[k] for (k, i) in enumerate(ind_comp); init=zero(OmegaP)) for f in f]
  phi = hom(Pr, OmegaP, img_gens)
  phi_res, _, _ = change_base_ring(R, phi)
  #set_attribute!(M, :show, show_kaehler_differentials)
  #set_attribute!(M, :is_kaehler_differential_module, (R, 1))
  return _cokernel(phi_res)
end

# TODO: Needed because of the different behavior of `cokernel` for different rings
_cokernel(f::ModuleFPHom{T}) where {U<:Union{MPolyRingElem, MPolyQuoRingElem}, T<:ModuleFP{U}} = cokernel(f)
_cokernel(f::ModuleFPHom) = cokernel(f)[1]

MPolyAnyRingElem = Union{MPolyRingElem, MPolyQuoRingElem, MPolyLocRingElem, MPolyQuoLocRingElem}

# The contraction given by the Euler vector field on graded modules
# over a standard graded ring.
function _euler_contraction(
    S::MPolyDecRing{<:MPolyAnyRingElem, <:MPolyRing}# S is a polynomial ring over an MPolyAnyRing
  )
  A = base_ring(S)::MPolyAnyRing
  # Build the relative Euler sequence
  W1 = kaehler_differentials(S)
  W0 = kaehler_differentials(S, 0)
  theta = hom(W1, W0, [x*W0[1] for x in gens(S)]; check=false)
end

function _euler_contraction(
    Q::MPolyQuoRing{<:MPolyDecRingElem{<:MPolyAnyRingElem}}, # S is a quotient ring over an MPolyAnyRing
  )
  S = base_ring(Q)::MPolyDecRing
  theta_S = _euler_contraction(S)
  theta_Q = _change_base_ring_and_preserve_gradings(Q, theta_S)
  return theta_Q
end

### Functionality for deformations in a blowup construction.
#
# Say the previous diagram was extended by a blowup 
#
#   Y₀ ↪ Y ↪ ℙʳ×𝔸ⁿ×𝔸¹
#   ↓    ↓     ↓
#   X₀ ↪ X ↪ 𝔸ⁿ×𝔸¹
#   ↓    ↓     ↓π
#  {0} ↪ 𝔸¹ ≅  𝔸¹
#
# with a relative projective space ℙʳ×𝔸ⁿ×𝔸¹ as ambient space and 
# one is interested in the graded module M representing the sheaf of 
# relative differential forms Ω¹_{Y/𝕜[t]}. Then this is *not* the 
# sheaf of relative differential forms over the whole base 𝔸ⁿ×𝔸¹ 
# (for which we have `relative_cotangent_module`)! 
#
# Hence, we need another, specialized method which allows us to 
# specify distinct variables in the ring of the base (i.e. the 
# `coefficient_ring` of the homogeneous coordinate ring of Y), 
# which are treated as parameters. 
#
# How this can be made user-facing in a reasonable way will have 
# to be discussed. But we would like to have this functionality 
# for the time being for internal purposes.
function _relative_kaehler_differentials(
    S::MPolyDecRing{<:MPolyAnyRingElem, <:MPolyRing}, # S is a polynomial ring over an MPolyAnyRing
    ind::Vector{Int} # Indices of variables in the `base_ring` to be omitted
  )
  A = base_ring(S)::MPolyAnyRing
  theta = _euler_contraction(S)
  WS, inc_WS = kernel(theta) # This graded module represents the Kaehler differentials in 
                             # the fiber directions of IP^n_A. 
  WA = _relative_kaehler_differentials(A, ind)
  pb_WA, pb_map = change_base_ring(S, WA)
  _make_zero_grading!(pb_WA)

  result = first(direct_sum(WS, pb_WA))
  return result
end

function _relative_kaehler_differentials(
    Q::MPolyQuoRing{<:MPolyDecRingElem{<:MPolyAnyRingElem}}, # S is a quotient ring over an MPolyAnyRing
    ind::Vector{Int} # Indices of variables in the `base_ring` to be omitted
  )
  S = base_ring(Q)::MPolyDecRing
  A = base_ring(S)::MPolyAnyRing
  theta = _euler_contraction(Q)
  W, inc_W = kernel(theta)
  WA, _ = change_base_ring(Q, _relative_kaehler_differentials(A, ind))
  WA = _make_zero_grading!(WA)
  WS = first(direct_sum(W, WA))
  # WS is a submodule of (⊕ⁿ⁺¹ S[-1]) ⊕ Ω¹_{A/k} with the first summand standing for 
  # the middle term of the Euler sequence. Since Proj(Q) ↪ Proj(S) is a closed 
  # embedding, we can compute the graded module for the Kaehler differentials of 
  # Proj(Q) by means of 
  f = gens(modulus(Q))
  F = ambient_free_module(WS)
  ind_comp = [i for i in 1:ngens(A) if !(i in ind)]
  img_gens = [sum(derivative(f, i)*F[i] for i in 1:ngens(S); init=zero(F)) + 
              sum(_coefficient_derivative(f, k)*F[i+ngens(S)] for (i, k) in enumerate(ind_comp); init=zero(F)) 
              for f in f]
  @assert all(v in WS for v in img_gens) "something went wrong with lifting the jacobian matrix"
  img_gens = vcat(img_gens, relations(WS))
  sub_mod, _ = sub(F, ambient_representatives_generators(WS))
  rel_mod, _ = sub(F, img_gens)
  return quo(sub_mod, rel_mod)
end

### Helper functions

### Take the partial derivative of a polynomial `f` w.r.t to the `i`-th variable 
# in the `coefficient_ring`.
function _coefficient_derivative(f::MPolyRingElem{<:MPolyAnyRingElem}, i::Int)
  S = parent(f)
  ctx = MPolyBuildCtx(S)
  for (e, a) in zip(AbstractAlgebra.exponent_vectors(f), AbstractAlgebra.coefficients(f))
    b = derivative(a, i)
    is_zero(b) && continue
    push_term!(ctx, b, e)
  end
  return finish(ctx)
end

### Given an ungraded module over a graded ring, endow it with the zero grading.
# This is used to process the output of `change_base_ring` when using a graded 
# ring `S` to pull back a module over its `coefficient_ring`. 
function _make_zero_grading!(F::FreeMod)
  is_graded(F) && error("module has a grading already")
  S = base_ring(F)
  is_graded(S) || error("no graded base ring")
  G = grading_group(S)
  o = zero(G)
  F.d = [o for _ in 1:ngens(F)] # TODO: Is there a clean way to do this?
  return F
end

function _make_zero_grading!(M::SubquoModule)
  F = ambient_free_module(M)
  _make_zero_grading!(F)
  @assert is_graded(M)
  return M
end

