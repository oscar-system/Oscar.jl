########################################################################
# Compute direct images of complexes of coherent sheaves
# 
# This is the implementation from arXiv:2602.14657 by F. Groh and 
# M. Zach, building on ideas of J. Weyman. 
# 
# Let `X` be a normal toric variety without torus factors, `S` its Cox 
# ring, and `R` an arbitrary commutative ring. For a (bounded) complex 
# `C*` of graded modules over the ring `S_R = S ⊗ R` one has the 
# associated complex of coherent sheaves `~C*` on `X × Spec R`. For 
# the projection `π : X × Spec R → Spec R` one can compute the direct 
# image `R π_* (~C*) ∈ Dᵇ(R-mod)`. The object `DirectImageComplex` 
# realizes the functor `𝕎 : Cᵇ(S_R-mod) → Cᵇ(R-mod)` which commutes 
# with the sheafification functor `~` and the canonical map to 
# `Dᵇ(R-mod)`. In other words, it makes the following diagram commute:
# 
#        Cᵇ(S_R-mod) ---𝕎-----> Cᵇ(R-mod)
#            |                     |
#            |~                    |
#            V                     V
#        Dᵇ(X×Spec R) --Rπ_*--> Dᵇ(R-mod)
########################################################################

include("generic_direct_images_types.jl")

### Production of the chains
# This samples from a context object which provides raw (truncated) direct images 
# for the twisting sheaves, together with their simplifications up to homotopy
# and the induced maps in the direct limit exhausting the Cech complexes. 
#   We aim at using three different types of context objects in parallel so that 
# the user can decide which one to use depending on their particular application.
struct DirectImageChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  wctx::WeymanCtx
  K::AbsHyperComplex
  ranges::Dict{Int, Vector{Vector{UnitRange{Int64}}}}

  function DirectImageChainFactory(pfctx::Union{ToricCtxWithParams, NewToricCtx, PushForwardCtx}, K::AbsHyperComplex)
    T = elem_type(target_ring(pfctx))
    ranges = Dict{Int, Vector{Vector{UnitRange{Int64}}}}()
    wctx = WeymanCtx(pfctx, Oscar.ReflectedComplex(K))
    return new{FreeMod{T}}(wctx, K, ranges)
  end
end

# Some additional getters to make unification easier. 
relative_dimension(ctx::ToricCtxWithParams) = dim(toric_variety(ctx))
relative_dimension(ctx::PushForwardCtx) = ngens(ctx.S) - rank(grading_group(ctx.S))
target_ring(ctx::ToricCtxWithParams) = ctx.R
target_ring(ctx::PushForwardCtx) = coefficient_ring(ctx.S)
target_ring(ctx::NewToricCtx) = coefficient_ring(cox_ring(toric_variety(ctx)))

# The terms of the direct image complex are (higher) direct images of the 
# twisting sheaves. More precisely, if the complex of graded modules for 
# the input is of the form 
#
#     Kᵖ = ⨁ ᵢS[-αₚᵢ], 
#
# then the direct image complex has terms 
# 
#     Cʳ = ⨁ ₚ₊ₛ₌ᵣ Rˢπ_* Kᵖ = ⨁ ₚ₊ₛ₌ᵣ ⨁ ᵢRˢπ_* 𝒪(-αₚᵢ)
#
# for the projection `π : X × Spec R → R`. 
#   The outer direct sum is for the "macro blocks", the inner one for the 
# "micro blocks". We keep track of the respective ranges, i.e. which range 
# of the indices is covered by which block, to enable fast mapping and 
# projection to the respective summands. 
function (fac::DirectImageChainFactory{ChainType})(self::AbsHyperComplex, I::Tuple) where {T, ChainType <: FreeMod{T}}
  i = first(I)
  wctx = fac.wctx
  ctx = pushforward_ctx(wctx)
  d = relative_dimension(ctx)
  ranges = Vector{Vector{UnitRange{Int64}}}()
  k = 0
  graded_complex = fac.K
  macro_offset = 0
  macro_summands = ChainType[]
  while k <= d
    if !can_compute_index(graded_complex, k + i)
      summand = FreeMod(target_ring(ctx), 0)
      push!(macro_summands, summand)
      push!(ranges, Vector{UnitRange{Int64}}())
      k += 1
      continue
    end
    micro_offset = 0
    micro_ranges = Vector{UnitRange{Int64}}()
    micro_summands = ChainType[]
    for (l, g) in enumerate(gens(graded_complex[i + k]))
      dd = degree(g; check=false)
      coh_mod = cohomology_model(ctx, -dd)
      str = coh_mod[-k]
push!(micro_summands, str)
      push!(micro_ranges, macro_offset+micro_offset+1:macro_offset+micro_offset+ngens(str))
      micro_offset += ngens(str)
    end
    macro_summand = is_empty(micro_summands) ? FreeMod(target_ring(ctx), 0) : direct_sum(micro_summands)[1]
    macro_offset += ngens(macro_summand)
    push!(macro_summands, macro_summand)
    push!(ranges, micro_ranges)
    k += 1
  end
  fac.ranges[i] = ranges
  return is_empty(macro_summands) ? FreeMod(target_ring(ctx), 0) : direct_sum(macro_summands)[1]
end

function can_compute(fac::DirectImageChainFactory, self::AbsHyperComplex, i::Tuple)
  return true
end

### Production of the morphisms 
struct DirectImageMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType} end

#   The real difficulty lies in computing the maps for the direct image. 
# The algorithm is described in our paper arXiv:2602.14657 and it proceeds 
# in a "sparse block format". To this end, we introduce data structures for 
# macro and micro blocks, which can then be assembled to sparse vectors of 
# blocks. This is not user facing and only used within this particular context. 
# Therefore, we do not care to adhere to the usual `parent`-convention in 
# Oscar. However, there are `parent`-like objects which let every block 
# know about its context and its position within the Cech double complex 
# underlying Weyman's construction. 
#
#   Everything happens in the following Cech double complex
#                   φ          φ
#    … → ̌Cᵖ⁺¹(Kⁿ⁺¹) → ̌Cᵖ⁺¹(Kⁿ) → ̌Cᵖ⁺¹(Kⁿ⁻¹) → …
#            ↑∂     φ    ↑∂     φ    ↑
#    … →  ̌Cᵖ(Kⁿ⁺¹)  →  ̌Cᵖ(Kⁿ)  →  ̌Cᵖ(Kⁿ⁻¹)  → …
#            ↑∂     φ    ↑∂     φ    ↑
#    … → ̌Cᵖ⁻¹(Kⁿ⁺¹) → ̌Cᵖ⁻¹(Kⁿ) → ̌Cᵖ⁻¹(Kⁿ⁻¹) → …
#            ↑           ↑           ↑
# The Cech terms are infinitely generated modules over the ring `S_R` and 
# hence modeled as a direct limit of finitely generated modules over this 
# ring. The index set for these direct limit is some form of `ℕᵏ`, somehow 
# tied to the exponent vectors for the denominators in the Cech complex. 
#   Moreover, every strand of a Cech term comes with a structure of a 
# double direct sum, i.e. with an outer and an inner block structure. 
# We wish to keep the index in the direct limit individual to the summands. 
# Hence, every macro block knows its position in the above double complex
# and every micro block knows its macro block and its index in the direct 
# limit. 
function (::DirectImageMapFactory)(self::AbsHyperComplex, p::Int, I::Tuple)
  i = first(I)
  _sign = is_even(i) ? -1 : 1
  @vprint :DirectImages 2 "computing outgoing map from $i\n"
  # fill the cache
  dom = self[i]
  cod = self[i-1]
  (is_zero(dom) || is_zero(cod)) && return hom(dom, cod, elem_type(cod)[zero(cod) for _ in 1:ngens(dom)])
  fac = chain_factory(self)
  wctx = fac.wctx
  img_gens_blocks = Vector{Vector{MacroVec}}() # A dense vector of macro blocks for 
                                               # the resulting map. These blocks might
                                               # be empty, so that the data structure stays 
                                               # effectively sparse. The macro block with 
                                               # the `R⁰π_*`-terms is the first one and it 
                                               # ends with the `Rᵈπ_*`-terms for `d = dim X`.
  ni = -i
  ranges = fac.ranges[i]
  img_gens = elem_type(cod)[]
  running_ind = 0
  gc = graded_complex(wctx)
  rel_d = relative_dimension(pushforward_ctx(wctx))
  # Iterating through the outer direct sum
  for (k, macro_range) in enumerate(ranges)
    p0 = ni - k + 1
    q0 = k - 1
    if !can_compute_index(gc, p0)
      # Domain is zero for this macro block. 
      # This means no generators and hence no images 
      # of generators to be added. We can just skip it.
      continue
    end
    pr_k = canonical_projection(dom, k)
    inj_k = canonical_injection(dom, k)
    block_k = codomain(pr_k)
    mac_mod = get_macro_block!(wctx, p0, q0, :cohomology)
    # Iterate through the inner direct sum
    for (j, micro_range) in enumerate(macro_range)
      pr_j = canonical_projection(block_k, j)
      inj_j = canonical_injection(block_k, j)
      # Iterate over the generators of the micro block
      for (l, g) in enumerate(gens(domain(inj_j)))
        running_ind += 1
        v = MacroVec(mac_mod, j, g) # the generator to be mapped
        w = apply_weyman_differential(v) # all the heavy lifting is done here
        result = zero(cod)
        if is_empty(w) # `w` is a `Vector` of `MacroVec`s
          push!(img_gens, result)
          continue
        end
        p_cod, q_cod = index(first(w)) # every `MacroVec` knows its cohomological 
                                       # position in the Cech double complex
        for (kk, mac_vec) in enumerate(w)
          inj_kk = canonical_injection(cod, kk + q_cod)
          cod_kk = domain(inj_kk)
          for (jj, mic_vec) in micro_vectors(mac_vec)
            inj_jj = canonical_injection(cod_kk, jj)
            result += inj_kk(inj_jj(mic_vec.c))
          end
        end
        push!(img_gens, result)
      end
    end
  end
  return hom(dom, cod, img_gens; check=false)
end

function can_compute(fac::DirectImageMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  return true
end

### The concrete struct
@attributes mutable struct DirectImageComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}

  function DirectImageComplex(
      pfctx::Union{PushForwardCtx, NewToricCtx, ToricCtxWithParams}, 
      K::AbsHyperComplex
    )
    T = elem_type(target_ring(pfctx))
    chain_fac = DirectImageChainFactory(pfctx, K)
    map_fac = DirectImageMapFactory{FreeModuleHom{FreeMod{T}, FreeMod{T}, Nothing}}()

    internal_complex = HyperComplex(1, chain_fac, map_fac, [:chain])
    return new{FreeMod{T}, FreeModuleHom{FreeMod{T}, FreeMod{T}, Nothing}}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::DirectImageComplex) = c.internal_complex

### functionality for the blocks in the weyman complex
is_zero(v::MicroVec) = is_zero(value(v))

function add!(v::MicroVec, w::MicroVec; check::Bool=false)
  return addmul!(v, w, 1; check)
end

function *(a, v::MicroVec)
  return MicroVec(macro_module(v), index(v), direct_limit_index(v), a*value(v))
end

function *(a, v::MacroVec)
  p, q = index(v)
  result = MacroVec(macro_module(v))
  for (_, w) in micro_vectors(v)
    addmul!(result, w, a)
  end
  return result
end

function addmul!(v::MicroVec, w::MicroVec, a; check::Bool=false)
  # When adding two `MicroVec`s one needs to bring everything to a 
  # common denominator; in other words: to a common index in the direct 
  # limit. This is supposed to happen on the fly. 
  @check macro_module(v) === macro_module(w)
  @check index(v) == index(v)
  d = direct_limit_index(v) 
  e = direct_limit_index(w)
  if d == e
    v.c += a*w.c
    return v
  end
  ctx = pushforward_ctx(v)
  alpha = degree(v) 
  p, q = index(macro_module(v))
  if all(a >= b for (a, b) in zip(d, e))
    L = typ(v) == :cochain ? ctx[d, e, alpha][-q] : induced_cohomology_map(ctx, d, e, alpha, -q)
    v.c = L(v.c) + a*w.c
    v.e = e
    return v
  elseif all(a <= b for (a, b) in zip(d, e))
    L = typ(v) == :cochain ? ctx[e, d, alpha][-q] : induced_cohomology_map(ctx, e, d, alpha, -q)
    v.c += a*L(w.c)
    v.e = d
    return v
  end

  de = Int[maximum(a) for a in zip(d, e)]
  L1 = typ(v) == :cochain ? ctx[d, de, alpha][-q] : induced_cohomology_map(ctx, d, de, alpha, -q)
  L2 = typ(v) == :cochain ? ctx[e, de, alpha][-q] : induced_cohomology_map(ctx, e, de, alpha, -q)
  c = L1(v.c) + a*L2(w.c)
  v.e = de
  v.c = c
  return v
end

function add!(mv::MacroVec, v::MicroVec; check::Bool=false)
  return addmul!(mv, v, 1; check)
end

function add!(v::MacroVec, w::MacroVec; check::Bool=false)
  return addmul!(v, w, 1; check)
end

function addmul!(v::MacroVec, w::MacroVec, c; check::Bool=false)
  @assert macro_module(v) === macro_module(w)
  for (j, ww) in micro_vectors(w)
    addmul!(v, ww, c; check)
  end
  return v
end

function addmul!(mv::MacroVec, v::MicroVec, a; check::Bool=false)
  @check macro_module(mv) === macro_module(v)
  is_zero(v) && return mv
  if !isdefined(mv, :micro_vecs)
    mv.micro_vecs = Tuple{Int, MicroVec}[(index(v), a*v)]
    return mv
  end
  j = index(v)
  if is_empty(mv.micro_vecs)
    push!(mv.micro_vecs, (j, a*v))
    return mv
  end
  ind = findfirst(i==j for (i, _) in mv.micro_vecs)
  if isnothing(ind)
    ind = findfirst(i>j for (i, _) in mv.micro_vecs)
    if isnothing(ind)
      push!(mv.micro_vecs, (j, a*v))
      return mv
    end
    insert!(mv.micro_vecs, ind, (j, a*v))
    return mv
  end
  _, w = micro_vectors(mv)[ind]
  w = addmul!(w, v, a; check)
  mv.micro_vecs[ind] = (j, w)
  is_zero(w) && deleteat!(micro_vectors(mv), ind)
  return mv
end

# Depending on the type of the context object, there are different index 
# sets for the direct limit. In each case we want a method to get the 
# zero index. This is provided here. 
_zero_denominator(ctx::PushForwardCtx) = Int[0 for _ in 1:ngens(graded_ring(ctx))]
_zero_denominator(ctx::ToricCtx) = Int[0 for _ in 1:length(cech_complex_generators(ctx))]
_zero_denominator(ctx::NewToricCtx) = Int[0 for _ in 1:ngens(graded_ring(ctx))]
_zero_denominator(ctx::ToricCtxWithParams) = _zero_denominator(ctx.pure_ctx)

# Given a `MacroVec` `v`, its `MicroVec`s may have different denominators
# (i.e. indices in the direct limit). This method computes a common index 
# into which all `MicroVec`s can be casted. 
function _common_denominator(v::MacroVec)
  e = _zero_denominator(pushforward_ctx(v))
  is_zero(v) && return e
  for (_, w) in micro_vectors(v)
    e = Int[maximum(a) for a in zip(e, direct_limit_index(w))]
  end
  return e
end

function _common_denominator(vs::Vector{MacroVec})
  is_empty(vs) && error("list must not be empty")
  v = first(vs)
  e = _zero_denominator(pushforward_ctx(v))
  for v in vs
    for (_, w) in micro_vectors(v)
      e = Int[maximum(a) for a in zip(e, direct_limit_index(w))]
    end
  end
  return e
end

# Brings all `MicroVec`s of a given `MacroVec` to the 
# same index in the direct limit for the Cech double complex.
function bring_to_common_denominator!(
    v::MacroVec, e::Vector{Int}=_common_denominator(v);
    check::Bool=false
  )
  ctx = pushforward_ctx(v)
  p, q = index(v)
  for (j, w) in micro_vectors(v)
    d = direct_limit_index(w)
    e == d && continue
    @check all(e >= d for (e, d) in zip(e, d)) "exponent vector too small"
    @check typ(v) == :cochain "common denominators for single macro blocks are only supported for cochains"
    L = ctx[d, e, degree(w)][-q] 
    w.e = e
    w.c = L(w.c)
  end
  return v
end

function bring_to_common_denominator(
    v::MacroVec, e::Vector{Int}=_common_denominator(v);
    check::Bool=false
  )
  ctx = pushforward_ctx(v)
  p, q = index(v)
  result = MacroVec(macro_module(v); check)
  is_empty(micro_vectors(v)) && return v
  result.micro_vecs = sizehint!(Tuple{Int, MicroVec}[], length(micro_vectors(v)))
  for (j, w) in micro_vectors(v)
    d = direct_limit_index(w)
    if e == d 
      push!(result.micro_vecs, (j, deepcopy(w)))
      continue
    end
    @check all(e >= d for (e, d) in zip(e, d)) "exponent vector too small"
    L = typ(v) == :cochain ? ctx[d, e, degree(w)][-q] : induced_cohomology_map(ctx, d, e, degree(w), -q)
    push!(result.micro_vecs, (j, MicroVec(macro_module(v), j, e, L(w.c); check)))
  end
  return result
end

# In the Cech double complex there are two maps, the horizontal ones 
# induced from those in the input complex `K`, and the vertical ones 
# in the Cech complexes. In addition, there are homotopy maps in the 
# vertical direction which are opposite to the Cech differentials. 
# We realize all these maps on the macro and micro blocks with their 
# individual positions and indices in the direct limit. 
function apply_phi(v::MacroVec; check::Bool=false)
  @assert typ(v) == :cochain "can not apply phi on cohomology"
  ctx = pushforward_ctx(v)
  wctx = weyman_ctx(v)
  p, q = index(v)
  gc = graded_complex(v)
  if !can_compute_map(gc, p)
    return nothing
  end
  f = map(gc, p)
  A = sparse_matrix(f)
  res_mod = get_macro_block!(wctx, p+1, q, :cochain)
  result = MacroVec(res_mod; check)
  for (i, w) in micro_vectors(v)
    is_zero(A[i]) && continue
    e = direct_limit_index(w)
    alpha = degree(w)
    for (j, a) in A[i]
      mm = multiplication_map(ctx, a, e, alpha, -q)
      ww = MicroVec(res_mod, j, e, mm(w.c); check)
      add!(result, ww; check)
    end
  end
  return result
end

function apply_phi(vs::Vector{MacroVec}; check::Bool=false)
  result = MacroVec[]
  for v in vs
    w = apply_phi(v; check)
    isnothing(w) && continue
    push!(result, w)
  end
  return result
end

function apply_cech_map(v::MacroVec; check::Bool=false)
  @assert typ(v) == :cochain "can not apply phi on cohomology"
  ctx = pushforward_ctx(v)
  p, q = index(v)
  wctx = weyman_ctx(v)
  res_mod = get_macro_block!(wctx, p, q+1, :cochain)
  result = MacroVec(res_mod; check)
  for (i, w) in micro_vectors(v)
    e = direct_limit_index(w)
    str = ctx[e, degree(w)]
    !can_compute_map(str, -q) && return nothing # If this happens for one, it happens for all
    partial = map(str, -q)
    add!(result, MicroVec(res_mod, i, e, partial(w.c); check); check)
  end
  return result
end

function apply_cech_map(vs::Vector{MacroVec}; check::Bool=false)
  result = MacroVec[]
  for v in vs
    w = apply_cech_map(v; check)
    isnothing(w) && break
    push!(result, w)
  end
  return result
end

function apply_homotopy_map(v::MacroVec; check::Bool=false)
  @assert typ(v) == :cochain "can not apply phi on cohomology"
  ctx = pushforward_ctx(v)
  wctx = weyman_ctx(v)
  p, q = index(v)
  is_zero(q) && return nothing
  res_mod = get_macro_block!(wctx, p, q-1, :cochain)
  result = MacroVec(res_mod; check)
  for (i, w) in micro_vectors(v)
    e = direct_limit_index(w)
    alpha = degree(w)
    str = ctx[e, alpha]
    # check whether the homotopy map exists
    !can_compute_index(str, -q+1) && return nothing
    h = simplified_strand_homotopy(ctx, e, alpha, -q)
    add!(result, MicroVec(res_mod, i, e, h(w.c); check); check)
  end
  return result
end

function apply_homotopy_map(vs::Vector{MacroVec}; check::Bool=false)
  result = MacroVec[]
  for v in vs
    w = apply_homotopy_map(v; check)
    isnothing(w) && continue
    push!(result, w)
  end
  return result
end


function apply_cohomology_inclusion(v::MacroVec; check::Bool=false)
  @assert typ(v) == :cohomology "can not apply cohomology inclusion to cochain"
  p, q = index(v)
  ctx = pushforward_ctx(v)
  wctx = weyman_ctx(v)
  res_mod = get_macro_block!(wctx, p, q, :cochain)
  result = MacroVec(res_mod; check)
  for (i, w) in micro_vectors(v)
    e = direct_limit_index(w)
    inc = simplified_strand_inclusion(ctx, e, degree(w), -q)
    inc_w = inc(w.c)
    mic_vec = MicroVec(res_mod, i, e, inc_w; check)
    add!(result, mic_vec; check)
  end
  return result
end

is_zero(v::MacroVec) = all(is_zero(v) for (_, v) in micro_vectors(v))

function _has_common_denominator(v::MacroVec)
  is_zero(v) && return true
  _, w = first(micro_vectors(v))
  e = direct_limit_index(w)
  return all(e == direct_limit_index(w) for (_, w) in micro_vectors(v))
end

# Given a `MacroVec`, it can be desirable to split it according to the 
# denominators (indices in the direct limit) used for its `MicroVec`s. 
function _split_by_denominators(v::MacroVec)
  result = Dict{Vector{Int}, MacroVec}()
  is_zero(v) && return result
  p, q = index(v)
  ctx = pushforward_ctx(v)
  for (j, w) in micro_vectors(v)
    e = direct_limit_index(w)
    mv = get!(result, e) do
      MacroVec(ctx, p, q, typ(v))
    end
    add!(mv, MicroVec(mv, j, e, w.c))
  end
  return result
end

function apply_I(v::MacroVec; check::Bool=false)
  @check typ(v) == :cohomology "can not apply I to cochain"
  p, q = index(v)
  sgn = (-1)^(p+q)
  res = Vector{MacroVec}[]
  is_zero(v) && return MacroVec[]

  for (alpha, vh) in sort_by_degree(v)
    # we can assume all `MicroVec`s to have the same exponent vector
    vv = apply_cohomology_inclusion(vh)
    result = [vv]
    for k in q:-1:0
      vv = apply_phi(vv)
      isnothing(vv) && break
      vv = apply_homotopy_map(-sgn*vv)
      isnothing(vv) && break
      push!(result, vv)
    end
    push!(res, reverse(result))
  end
  @check begin
    rr = sum(res)
    e = _common_denominator(rr)
    w = apply_P([bring_to_common_denominator(v, e) for v in rr])
    @assert is_zero(addmul(v, w[end], -1))
    all(is_zero(w[1:end-1]))
  end
  return sum(res)
end

function apply_cohomology_projection(v::MacroVec; check::Bool=false)
  @assert typ(v) == :cochain "can not apply cohomology projection to cohomology"
  p, q = index(v)
  ctx = pushforward_ctx(v)
  wctx = weyman_ctx(v)
  res_mod = get_macro_block!(wctx, p, q, :cohomology)
  result = MacroVec(res_mod; check)
  for (i, w) in micro_vectors(v)
    e = direct_limit_index(w)
    inc = simplified_strand_projection(ctx, e, degree(w), -q)
    add!(result, MicroVec(res_mod, i, e, inc(w.c); check); check)
  end
  return result
end

# The maps `I` and `P` are the inclusions and projections to and from 
# cohomology of the twisting sheaves. These occur for every index 
# in the direct limit of the Cech complex.
function apply_P(v::MacroVec; check::Bool=false)
  @assert typ(v) == :cochain "can not apply P to cohomology"
  p, q = index(v)
  sgn = (-1)^(p+q)
  @check _has_common_denominator(v) "only elements with common denominators allowed"

  # we can assume all `MicroVec`s to have the same exponent vector
  vv = v
  inter = MacroVec[vv]
  for k in q-1:-1:0
    vv = apply_homotopy_map(vv)
    vv = -sgn*vv
    vv = apply_phi(vv)
    push!(inter, vv)
  end
  return MacroVec[apply_cohomology_projection(w; check) for w in reverse(inter)]
end

function apply_P(vs::Vector{MacroVec}; check::Bool=false)
  @assert all(typ(v) == :cochain for v in vs) "can not apply P to cohomology"
  is_empty(vs) && return MacroVec[]
  p, q = index(first(vs))
  sgn = (-1)^(p+q)
  @check all(sum(index(v)) == p+q for v in vs) "not an element in a single cohomological degree"
  @check all(index(vs[k])[2] == index(vs[k+1])[2]-1 for k in 1:length(vs)-1) "blocks must be consecutive on the anti-diagonal, but we found $(index.(vs))"
  @check all(_has_common_denominator(v) for v in vs) "only elements with common denominators allowed"
  if is_zero(vs)
    return MacroVec[]
    result = MacroVec[]
    for v in vs
      p, q = index(v)
      bm = get_macro_block!(weyman_ctx(v), p, q, :cohomology)
      push!(result, MacroVec(bm; check))
    end
    p0, q0 = index(last(vs))
    deg = p0 + q0
    gc = graded_complex(weyman_ctx(first(vs)))
    for q in q0-1:-1:0
      p = deg - q
      !can_compute_index(gc, p) && break
      bm = get_macro_block!(weyman_ctx(last(vs)), p, q, :cohomology)
      push!(result, MacroVec(bm; check))
    end
    return reverse(result)
  end

  vvs = filter(!is_zero, vs)
  e = _common_denominator(first(vvs))
  @check all(e == _common_denominator(v) for v in vvs) "only elements with common denominators allowed"

  # we can assume all `MicroVec`s to have the same exponent vector
  vv = vs[end]
  inter = MacroVec[vv]
  for v in reverse(vs)[2:end]
    vv = apply_homotopy_map(vv)
    vv = apply_phi(vv)
    @check index(vv) == index(v)
    vv = addmul(v, vv, -sgn)
    push!(inter, vv)
  end
  return MacroVec[apply_cohomology_projection(w) for w in reverse(inter)]
end

# Functionality for the convenience of debugging. The user will probably never see this.
function Base.show(io::IO, mv::MacroVec)
  println(io, "macro vector of type $(typ(mv)) with index $(index(mv)) and components")
  for (_, v) in micro_vectors(mv)
    print(io, v)
  end
end

function Base.show(io::IO, mv::MicroVec)
  println(io, "$(index(mv)) -> $(direct_limit_index(mv)): $(mv.c)")
end

function deepcopy_internal(v::MicroVec, d::IdDict)
  return MicroVec(v.macro_mod, deepcopy_internal(v.j, d),
                  deepcopy_internal(v.e, d), deepcopy_internal(v.c, d))
end

# Auxillary further implementations
function deepcopy_internal(v::MacroVec, d::IdDict) 
  result = MacroVec(v.macro_mod)
  mic_vecs = deepcopy_internal(micro_vectors(v), d)
  result.micro_vecs = mic_vecs
  return result
end

function addmul(v::MicroVec, w::MicroVec, a; check::Bool=false)
  return addmul!(deepcopy(v), w, a; check)
end

function +(v::MicroVec, w::MicroVec)
  return addmul!(deepcopy(v), w, 1)
end

function addmul(v::MacroVec, w::MicroVec, a; check::Bool=false)
  return addmul!(deepcopy(v), w, a; check)
end

function +(v::MacroVec, w::MicroVec)
  return addmul!(deepcopy(v), w, 1)
end

function addmul(v::MacroVec, w::MacroVec, a; check::Bool=false)
  return addmul!(deepcopy(v), w, a; check)
end

function +(v::MacroVec, w::MacroVec)
  return addmul!(deepcopy(v), w, 1)
end

function add!(a::Vector{MacroVec}, b::Vector{MacroVec}; check::Bool=false)  
  return addmul!(a, b, 1; check)
end

function addmul!(a::Vector{MacroVec}, b::Vector{MacroVec}, c; check::Bool=false)
  is_empty(a) && return c*b                      
  is_empty(b) && return a                                
  @check begin
    v0 = first(a)                                          
    w0 = first(b)                                          
    @assert weyman_ctx(v0) === weyman_ctx(w0)              
    @assert all(weyman_ctx(v) === weyman_ctx(v0) for v in   a)
    @assert all(weyman_ctx(w) === weyman_ctx(w0) for w in   b)
    deg = sum(index(v0))                                   
    @assert deg == sum(index(w0))                          
    @assert all(deg == sum(index(v)) for v in a)           
    @assert all(deg == sum(index(v)) for v in b)
  end
  for mv in b
    ind = findfirst(index(mv) == index(mw) for mw in a)
    if isnothing(ind) 
      push!(a, c*mv)
      continue
    end
    a[ind] = addmul!(a[ind], mv, c)
  end
  return sort!(a; by=v->index(v)[2])
end

function MacroVec(v::MicroVec; check::Bool=false)
  result = MacroVec(macro_module(v); check)
  result.micro_vecs = [(index(v), v)]
  return result
end

function MicroVec(mb::MacroMod, j::Int, c::FreeModElem; check::Bool=false)
  gc = graded_complex(mb)
  p, q = index(mb)
  alpha = -degrees_of_generators(gc[p])[j]
  e = _minimal_exponent_vector(pushforward_ctx(mb), alpha)
  return MicroVec(mb, j, e, c; check)
end

function MicroVec(mb::MacroMod, j::Int; check::Bool=false)
  gc = graded_complex(mb)
  p, q = index(mb)
  alpha = -degrees_of_generators(gc[p])[j]
  e = _minimal_exponent_vector(pushforward_ctx(mb), alpha)
  return MicroVec(mb, j, e; check)
end

function MicroVec(mb::MacroMod, j::Int, k::Int; check::Bool=false)
  gc = graded_complex(mb)
  p, q = index(mb)
  alpha = -degrees_of_generators(gc[p])[j]
  e = _minimal_exponent_vector(pushforward_ctx(mb), alpha)
  return MicroVec(mb, j, e, k; check)
end

function MacroVec(mb::MacroMod, j::Int, c::FreeModElem; check::Bool=false)
  return MacroVec(MicroVec(mb, j, c); check)
end

function MacroVec(mb::MacroMod, j::Int, e::Vector{Int}, c::FreeModElem; check::Bool=false)
  return MacroVec(MicroVec(mb, j, e, c; check); check)
end

function is_zero(mb::MacroMod)
  p, q = index(mb)
  ctx = pushforward_ctx(mb)
  gc = graded_complex(mb)
  gm = gc[p]
  for (i, v) in enumerate(gens(gm))
    alpha = -degree(v; check=false)
    coh = cohomology_model(ctx, alpha)[-q]
    !is_zero(coh) && return false
  end
  return true
end

function ==(v::MacroVec, w::MacroVec)
  weyman_ctx(v) === weyman_ctx(w) || error("incompatible macro vectors")
  index(v) == index(w) || return false
  return is_zero(addmul(v, w, -1))
end

# The heavy lifting for the differentials in the Weyman complex.
# 
# The procedure is to "go down the staircase" by mapping any 
# generator first by `φ`, then by the homotopy map `h`, 
# then by `φ` again, and so on. On each step of the staircase 
# we collect the intermediate result (with an appropriate sign) 
# and then combine these in a direct sum for the image. 
# Thanks to the individual steps being wrapped in the above 
# functions, it is now relatively easy to write down this 
# procedure.
function apply_weyman_differential(v::MacroVec; check::Bool=false)
  p, q = index(v)
  sgn = (-1)^(p+q)
  u = apply_cohomology_inclusion(sgn*v)
  w = apply_phi(u; check)
  if isnothing(w) 
    # The image is zero, because the complex stops there. 
    # But we need to build a `Vector` of macro vectors
    # for all relevant degrees in order to be compatible with the further 
    # processing and mapping to the direct sums. 
    result = MacroVec[]
    gc = graded_complex(v)
    wctx = weyman_ctx(v)
    for qq in q+1:relative_dimension(pushforward_ctx(wctx))
      pp = p + q - qq + 1
      !can_compute_index(gc, pp) && continue
      push!(result, MacroVec(get_macro_block!(weyman_ctx(v), pp, qq, :cohomology)))
    end
    return result
  end
  result = MacroVec[apply_cohomology_projection(w; check)]
  for k in q:-1:1
    all_zero = true
    wctx = weyman_ctx(v)
    gc = graded_complex(wctx)
    nn = p+q+1
    for kk in k-1:-1:0
      pp = nn-kk
      !can_compute_index(gc, pp) && break
      mb = get_macro_block!(wctx, pp, kk, :cohomology)
      if !is_zero(mb)
        all_zero = false
        break
      end
    end

    if all_zero
      for kk in k-1:-1:0
        pp = nn - kk
        !can_compute_index(gc, pp) && continue
        push!(result, MacroVec(get_macro_block!(weyman_ctx(v), pp, kk, :cohomology)))
      end
      return project_to_weyman_complex!(reverse(result); check)
    end
    w = apply_homotopy_map(w; check)
    isnothing(w) && break
    w = -sgn*w
    w = apply_phi(w; check)
    isnothing(w) && break
    push!(result, apply_cohomology_projection(w; check))
  end
  return project_to_weyman_complex!(reverse(result); check)
end

# The components of a `MacroVec` might live in different 
# strands, which require different treatment each. This 
# decomposes a given `MacroVec` into its homogeneous 
# components. 
function sort_by_degree(mv::MacroVec; check::Bool=false)
  ctx = pushforward_ctx(mv)
  G = grading_group(graded_ring(ctx))
  result = Dict{elem_type(G), MacroVec}()
  is_zero(mv) && return result
  for (j, v) in micro_vectors(mv)
    alpha = degree(v)
    bucket = get!(result, alpha) do
      MacroVec(macro_module(mv); check)
    end
    # things come in the correct order here
    push!(micro_vectors(bucket), (j, v))
  end
  return result
end

# The Weyman complex has as its terms the cohomology of direct 
# sums of twisting sheaves. The models for this cohomology depend
# on the twist, i.e. we need to choose the index in the direct limit
# depending on the twist. Consequently, we need to project an element 
# from the Cech double complex to the correct cohomology model. 
#   When the `NewToricCtx` is used, this is easy, because the cohomology 
# models form a direct limit over the same index set. If not, then we 
# need to use the artificially constructed direct limit over `ℕ` which 
# is described in our paper. 
function project_to_weyman_complex!(vs::Vector{MacroVec}; check::Bool=false)
  #vsc = deepcopy(vs)
  @check all(typ(v) == :cohomology for v in vs) "can not apply P to cohomology"
  is_empty(vs) && return MacroVec[]
  p, q = index(first(vs))
  @check all(sum(index(v)) == p+q for v in vs) "not an element in a single cohomological degree"
  wctx = weyman_ctx(first(vs))
  ctx = pushforward_ctx(wctx)

  if ctx isa ToricCtxWithParams && ctx.pure_ctx isa NewToricCtx # In this case homotopies are compatible with direct limits
    part = MacroVec[]
    for mv in vs
      p, q = index(mv)
      res_mac_vec = MacroVec(get_macro_block!(wctx, p, q, :cohomology); check)
      for (alpha, w) in sort_by_degree(mv)
        e0 = _minimal_exponent_vector(ctx, alpha)
        for (j, v) in micro_vectors(w)
          e = direct_limit_index(v)
          lambda = induced_cohomology_map(ctx, e, e0, alpha, -q)
          lambda_v = lambda(value(v))
          add!(res_mac_vec, MicroVec(macro_module(res_mac_vec), j, e0, lambda_v))
        end
      end
      push!(part, res_mac_vec)
    end
    return part
  end

  # Take the long way with artificially produced direct limits
  @check all(index(vs[k])[2] == index(vs[k+1])[2]-1 for k in 1:length(vs)-1) "blocks must be consecutive on the anti-diagonal"
  @check all(_has_common_denominator(v) for v in vs) "only elements with common denominators allowed"
  vvs = filter(!is_zero, vs)
  is_empty(vvs) && return MacroVec[]
  e = _common_denominator(first(vvs))
  @check all(e == _common_denominator(v) for v in vvs) "only elements with common denominators allowed"
  
  part = MacroVec[]
  for mv in vvs
    p, q = index(mv)
    for (j, v) in micro_vectors(mv)
      img_gens = weyman_projections(wctx, p, q, j, first(e))
      part = add!(part, sum(c.*img_gens[l] for (l, c) in coordinates(value(v))); check)
    end
  end
  return part

  result = MacroVec[]
  p0, q0 = index(first(vs))
  wctx = weyman_ctx(first(vs))
  gc = graded_complex(wctx)
  for k in q0-1:-1:0
    !can_compute_index(gc, p0+q0-k) && break
    mb = get_macro_block!(wctx, p0+q0-k, k, :cohomology)
    pushfirst!(vs, MacroVec(mb; check))
  end

  for (k, v) in enumerate(reverse(vs))
    ctx = pushforward_ctx(v)
    p, q = index(v)
    res_mod = get_macro_block!(weyman_ctx(v), p, q, :cohomology)
    res = MacroVec(res_mod)
    for (alpha, w) in sort_by_degree(v)
      d = _minimal_exponent_vector(ctx, alpha)
      if d == e
        add!(res, w)
        continue
      end
      str_d = simplified_strand(ctx, d, alpha)
      str_e = simplified_strand(ctx, e, alpha)
      LH_inv = induced_cohomology_map(ctx, e, d, alpha, -q)
      wd = MacroVec(macro_module(w))
      wd.micro_vecs = [(j, MicroVec(macro_module(w), j, d, LH_inv(value(ww)))) for (j, ww) in micro_vectors(w)]
      add!(res, wd)
      Icd = apply_I(wd)
      @assert all(d == _common_denominator(v) for v in filter(!is_zero, Icd))
      LIcd = MacroVec[bring_to_common_denominator(w, e) for w in Icd]
      @assert all(e == _common_denominator(v) for v in filter(!is_zero, LIcd))
      PLIcd = apply_P(LIcd)
      @assert all(e == _common_denominator(v) for v in filter(!is_zero, PLIcd))
      head = PLIcd[end]
      @assert w == head
      tail = PLIcd[1:end-1]
      is_zero(tail) && continue
      vs[1:end-k] = addmul!(vs[1:end-k], tail, -1)
    end
    push!(result, res)
  end
  return reverse(result)
end

# When the cohomology models for the twisting sheaves form no natural 
# direct limit, we need to construct an artificial direct limit over 
# `ℕ` where the map `lₖᵏ⁺ⁿ` from index `k` to `k+n` is realized as the 
# composition `lₖ₊ₙ₋₁ᵏ⁺ⁿ∘…∘lₖ₊₁ᵏ⁺² ∘lₖᵏ⁺¹`. This is, of course, expensive 
# and a lot of these constructions should again be cached. The code below 
# is for that. 

### Inclusions and projections to the Weyman complex
function _weyman_inc_dict(wctx::WeymanCtx)
  if !isdefined(wctx, :weyman_inclusions)
    wctx.weyman_inclusions = Dict{Tuple{Int, Int, Int}, Dict}()
  end
  return wctx.weyman_inclusions
end

# Construct the inclusion of the `(p, q, j)`-th micro block into W_{\leq e}^\bullet.
# This returns a `Vector` of `Vector`s of `MacroVec`s which are the images of the 
# generators. 
function weyman_inclusions(wctx::WeymanCtx, p::Int, q::Int, j::Int, e::Int; check::Bool=false)
  inc_dict = get!(_weyman_inc_dict(wctx), (p, q, j)) do
    Dict{Int, Vector}()
  end::Dict{Int, <:Vector}
  return get!(inc_dict, e) do
    _build_weyman_inclusions(wctx, p, q, j, e; check)
  end::Vector{Vector{MacroVec}}
end

function _build_weyman_inclusions(wctx::WeymanCtx, p::Int, q::Int, j::Int, e::Int; check::Bool=false)
  mac_mod = get_macro_block!(wctx, p, q)
  mic_zero = MicroVec(mac_mod, j; check)
  alpha = degree(mic_zero; check=false)
  d = _minimal_exponent_vector(pushforward_ctx(wctx), alpha)
  result = Vector{MacroVec}[]

  # If this is the minimal degree, just return the inclusion
  if first(d) == e
    str = simplified_strand(pushforward_ctx(wctx), d, alpha)
    dom = str[-q]
    res = MacroVec[MacroVec(MicroVec(mac_mod, j, k; check); check) for k in 1:ngens(dom)]
    push!(result, res)
    return result
  end

  # Otherwise compose from what's already there
  low_incs = weyman_inclusions(wctx, p, q, j, e-1)
  for (k, low_img_gens) in enumerate(low_incs)
    # TODO: Replace by a method for `apply_I` on `Vector`s!
    img_gen = sum(apply_I(v) for v in low_img_gens)
    img_gen = MacroVec[bring_to_common_denominator(v, e) for v in img_gen]
    img_gen = apply_P(img_gen)
    push!(result, img_gen)
  end
  return result
end

function _weyman_pr_dict(wctx::WeymanCtx)
  if !isdefined(wctx, :weyman_projections)
    wctx.weyman_projections = Dict{Tuple{Int, Int, Int}, Dict}()
  end
  return wctx.weyman_projections
end

# Construct the inclusion of the `(p, q, j)`-th micro block into W_{\leq e}^\bullet.
# This returns a `Vector` of `Vector`s of `MacroVec`s which are the images of the 
# generators. 
function weyman_projections(wctx::WeymanCtx, p::Int, q::Int, j::Int, e::Int; check::Bool=false)
  pr_dict = get!(_weyman_pr_dict(wctx), (p, q, j)) do
    Dict{Int, Vector}()
  end::Dict{Int, <:Vector}
  return get!(pr_dict, e) do
    _build_weyman_projection(wctx, p, q, j, e; check)
  end::Vector{Vector{MacroVec}}
end

function _build_weyman_projection(wctx::WeymanCtx, p::Int, q::Int, j::Int, e::Int; check::Bool=false)
  mac_mod = get_macro_block!(wctx, p, q, :cohomology)
  mic_zero = MicroVec(mac_mod, j; check)
  alpha = degree(mic_zero)
  d = _minimal_exponent_vector(pushforward_ctx(wctx), alpha)
  # reconstruct the whole exponent vector (fragile!)
  ee = Int[e for _ in 1:length(d)]
  result = Vector{MacroVec}[]
  str_d = simplified_strand(pushforward_ctx(wctx), d, alpha)
  dom_d = str_d[-q]

  # If this is the minimal degree, just return the inclusion
  if first(d) == e
    p, q = index(mac_mod)
    rng = [(p+q-i, i) for i in q-1:-1:0 if can_compute_index(graded_complex(wctx), p+q-i)]
    tail = MacroVec[MacroVec(get_macro_block!(wctx, p, q, :cohomology); check) for (p, q) in reverse(rng)]
    for k in 1:ngens(dom_d)
      r = MacroVec(MicroVec(mac_mod, j, ee, k; check); check)
      push!(result, push!(copy(tail), r))
    end
    return result
  end
  
  if first(d) > e
    error("projection to Weyman complex does not exist for this index")
  end

  # Otherwise, go one step down and compose from what's already there
  str_ee = simplified_strand(pushforward_ctx(wctx), ee, alpha)
  dom_ee = str_ee[-q]
  d = [e-1 for _ in 1:length(d)]
  for k in 1:ngens(dom_ee) #(k, v) in enumerate(gens(dom_ee))
    mic_vec = MicroVec(mac_mod, j, ee, k; check)
    mac_vec = MacroVec(mic_vec; check)
    #interm = MacroVec[]
    interm = [mac_vec]
    img = MacroVec[]
    p, q = index(mac_mod)
    for qq in q:-1:0
      pp = p + q - qq
      if !can_compute_index(graded_complex(wctx), pp)
        break
      end
      top_vec = only(v for v in interm if index(v) == (pp, qq))
      res = MacroVec(macro_module(top_vec); check)
      for (beta, v) in sort_by_degree(top_vec)
        #@assert all(degree(vv) == beta for (_, vv) in micro_vectors(v))
        LH_inv = induced_cohomology_map(pushforward_ctx(wctx), ee, d, beta, -qq)
        #phi = hom(domain(LH_inv), codomain(LH_inv), LH_inv.(gens(domain(LH_inv))))
        #@show phi
        #@assert is_isomorphism(phi)
        vv = MacroVec(macro_module(top_vec); check)
        vv.micro_vecs = [(j, MicroVec(macro_module(top_vec), j, d, LH_inv(value(ww)); check)) for (j, ww) in micro_vectors(v)]
        res = add!(res, vv)
      end
      push!(img, res)
      I_res = apply_I(res)
      LI_res = MacroVec[bring_to_common_denominator(w, ee) for w in I_res]
      PLI_res = apply_P(LI_res)
      interm = addmul!(interm, PLI_res, -1)
      @assert is_zero(only(v for v in interm if index(v) == (pp, qq)))
    end
    push!(result, reverse(img))
  end

  # now `result` holds the projections of the generators in the (p, q, j)-block 
  # in index d = e-1 

  final = Vector{MacroVec}[]
  for mvs in result
    part = MacroVec[]
    for mv in mvs
      p, q = index(mv)
      for (j, v) in micro_vectors(mv)
        img_gens = weyman_projections(wctx, p, q, j, e-1)
        @assert length(img_gens) == ngens(parent(value(v)))
        part = add!(part, sum(c.*img_gens[l] for (l, c) in coordinates(value(v))); check)
      end
    end
    push!(final, part)
  end
  return final
end

########################################################################
# Functionality for NewToricCtx
########################################################################
toric_variety(ctx::NewToricCtx) = ctx.X
graded_ring(ctx::NewToricCtx) = ctx.S


function ring_as_hypercomplex(ctx::NewToricCtx)
  if !isdefined(ctx, :S1)
    S = graded_ring(ctx)
    ctx.S1 = ZeroDimensionalComplex(graded_free_module(S, [zero(grading_group(S))]))
  end
  return ctx.S1
end

function cech_complex_generators(ctx::NewToricCtx)
  if !isdefined(ctx, :cech_gens)
    ctx.cech_gens = gens(irrelevant_ideal(toric_variety(ctx)))
  end
  return ctx.cech_gens::Vector{elem_type(graded_ring(ctx))}
end

# obtain all monomials of a given degree in the Cox ring. 
function all_monomials(ctx::NewToricCtx, alpha::FinGenAbGroupElem)
  return get!(ctx.all_monomial_cache, alpha) do
    S = graded_ring(ctx)
    collect(all_exponents(S, alpha))
  end::Vector{Vector{Int}}
end

# obtain a dictionary for the inverse mapping of the above
function all_monomials_inv(ctx::NewToricCtx, alpha::FinGenAbGroupElem)
  return get!(ctx.all_monomial_inv_dicts, alpha) do
    Dict{Vector{Int}, Int}(e=>i for (i, e) in enumerate(all_monomials(ctx, alpha)))
  end
end

# get the Cox ring with its fine grading
function fine_graded_ring(ctx::NewToricCtx)
  if !isdefined(ctx, :fine_graded_ring)
    S = graded_ring(ctx)
    P = forget_grading(S)
    F = free_abelian_group(ngens(P))
    FS, _ = grade(P, gens(F))
    ctx.fine_graded_ring = FS
  end
  return ctx.fine_graded_ring::MPolyDecRing
end

# get the sample complex Hom(P*, S) for the resolution P* of the `irrelevant_ideal`. 
function sample_complex(ctx::NewToricCtx)
  if !isdefined(ctx, :sample_complex)
    S = fine_graded_ring(ctx)
    X = toric_variety(ctx)
    I = ideal(S, elem_type(S)[S(forget_grading(g)) for g in gens(irrelevant_ideal(X))])
    res, _ = free_resolution(SimpleFreeResolution, I)
    FG = grading_group(S)
    ctx.sample_complex = hom(res, ZeroDimensionalComplex(graded_free_module(S, [zero(FG)])))
  end
  return ctx.sample_complex
end

# selects directly the sector for this degree
function fine_strand(ctx::NewToricCtx, e::FinGenAbGroupElem)
  return fine_strand(ctx, Int[Int(e[i]) for i in 1:ngens(parent(e))])
end

function fine_strand(ctx::NewToricCtx, e::Vector{Int})
  k = sum(((c<0) << (k-1)) for (k, c) in enumerate(e); init=0)
  return fine_strand(ctx, k)
end

# p is considered as a binary vector with its i-th digit ==1 if 
# and only if the imagined exponent vector has its i-th component < 0.
function fine_strand(ctx::NewToricCtx, p::Int)
  return get!(ctx.fine_strands, p) do
    FS = fine_graded_ring(ctx)
    FG = grading_group(FS)
    deg = sum(-((p >> (k-1))%2)*g for (k, g) in enumerate(gens(FG)); init=zero(FG))
    res, _ = strand(sample_complex(ctx), deg; check=false)
    return res
  end
end

function simplified_fine_strand(ctx::NewToricCtx, p::Int)
  return get!(ctx.simplified_fine_strands, p) do
    simplify(fine_strand(ctx, p); with_homotopy_maps=true)
  end
end

function simplified_fine_strand(ctx::NewToricCtx, e::FinGenAbGroupElem)
  return simplified_fine_strand(ctx, Int[e[i] for i in 1:ngens(parent(e))])
end

function simplified_fine_strand(ctx::NewToricCtx, e::Vector{Int})
  k = sum(((c<0) << (k-1)) for (k, c) in enumerate(e); init=0)
  return simplified_fine_strand(ctx, k)
end

function getindex(ctx::NewToricCtx, e::Vector{Int}, alpha::FinGenAbGroupElem)
  return getindex(ctx, first(e), alpha)
end

function getindex(ctx::NewToricCtx, e::Int, alpha::FinGenAbGroupElem)
  G = parent(alpha)
  S = graded_ring(ctx)
  kk = coefficient_ring(S)
  @assert G === grading_group(S)
  return get!(ctx.strands, (e, alpha)) do
    beta = e*sum(degree(x; check=false) for x in gens(S); init=zero(G))
    all_mons = all_monomials(ctx, alpha+beta)
    return DirectSumComplex(kk, [:chain], [fine_strand(ctx, ee.-e) for ee in all_mons])
  end
end

# dirty catch of the edge case with an empty list of summands
function DirectSumComplex(
    R::Ring, dirs::Vector{Symbol}, summands::Vector
  )
  @assert isempty(summands) "non-empty list of summands but no useful type of the entries of that list; try preparing your list in a type-stable way"
  ET = elem_type(R)
  return DirectSumComplex(R, dirs, AbsHyperComplex{OFPModule{ET}, OFPModuleHom}[])
end

function simplified_strand(ctx::NewToricCtx, e::Vector{Int}, alpha::FinGenAbGroupElem)
  return simplified_strand(ctx, first(e), alpha)
end

function simplified_strand(ctx::NewToricCtx, e::Int, alpha::FinGenAbGroupElem)
  return get!(ctx.simplified_strands, (e, alpha)) do
    G = parent(alpha)
    S = graded_ring(ctx)
    kk = coefficient_ring(S)
    @assert G === grading_group(S)
    beta = e*sum(degree(x; check=false) for x in gens(S); init=zero(G))
    all_mons = all_monomials(ctx, alpha+beta)
    summands = [simplified_fine_strand(ctx, ee .- e ) for ee in all_mons]
    @assert length(all_mons) == length(summands)
    return DirectSumComplex(kk, [:chain], summands)
  end
end

function induced_cohomology_map(
    ctx::NewToricCtx, e0::Vector{Int},
    e1::Vector{Int}, alpha::FinGenAbGroupElem,
    i::Int
  )
  return induced_cohomology_map(ctx, first(e0), first(e1), alpha, i)
end

function induced_cohomology_map(
    ctx::NewToricCtx, e0::Int,
    e1::Int, alpha::FinGenAbGroupElem,
    i::Int
  )
  return get!(ctx.induced_cohomology_maps, (e0, e1, alpha, i)) do 
    dom = simplified_strand(ctx, e0, alpha)[i]
    e0 == e1 && return identity_map(dom)
    cod_str = simplified_strand(ctx, e1, alpha)
    cod = cod_str[i]
    S = graded_ring(ctx)
    G = grading_group(S)
    delta = sum(degree(x; check=false) for x in gens(S); init=zero(G))
    if e0 <= e1
      all_mons_dom = all_monomials(ctx, alpha + e0*delta)
      all_mons_cod_inv = all_monomials_inv(ctx, alpha + e1*delta)
      index_mapping = Int[all_mons_cod_inv[ee.+(e1 - e0)] for ee in all_mons_dom]
      img_gens = elem_type(cod)[]
      for (i, k) in enumerate(index_mapping)
        #pr = canonical_projection(dom, i)
        inc = canonical_injection(cod, k)
        #@assert codomain(pr) === domain(inc)
        img_gens = vcat(img_gens, images_of_generators(inc))
      end
      return hom(dom, cod, img_gens)
    elseif e0 > e1
      @assert e1 >= _minimal_exponent_vector(ctx, alpha)
      orig = induced_cohomology_map(ctx, e1, e0, alpha, i)
      @assert ngens(domain(orig)) == ngens(codomain(orig))
      return inv(induced_cohomology_map(ctx, e1, e0, alpha, i))
    else
      error("not implemented")
    end
  end
end

function simplified_strand_homotopy(
    ctx::NewToricCtx, e::Vector{Int}, alpha::FinGenAbGroupElem, p::Int
  )
  return simplified_strand_homotopy(ctx, first(e), alpha, p)
end

function simplified_strand_homotopy(
    ctx::NewToricCtx, e::Int, alpha::FinGenAbGroupElem, p::Int
  )
  return get!(ctx.simplified_strand_homotopies, (e, alpha, p)) do
    dom = ctx[e, alpha][p]
    cod = ctx[e, alpha][p+1]
    S = graded_ring(ctx)
    kk = coefficient_ring(S)
    G = grading_group(S)
    delta = sum(degree(x; check=false) for x in gens(S); init=zero(G))
    all_mons = all_monomials(ctx, alpha+e*delta)
    maps = [homotopy_map(simplified_fine_strand(ctx, ee.-e), p) for ee in all_mons]
    #for (k, phi) in enumerate(maps)
      #@assert domain(phi) === codomain(canonical_projection(dom, k))
      #@assert codomain(phi) === codomain(canonical_projection(cod, k))
    #end
    return _direct_sum(kk, maps; domain=dom, codomain=cod)
  end
end

function simplified_strand_inclusion(
    ctx::NewToricCtx, e::Vector{Int}, alpha::FinGenAbGroupElem, p::Int
  )
  return simplified_strand_inclusion(ctx, first(e), alpha, p)
end

function simplified_strand_inclusion(
    ctx::NewToricCtx, e::Int, alpha::FinGenAbGroupElem, p::Int
  )
  return get!(ctx.simplified_strand_inclusions, (e, alpha, p)) do
    dom = simplified_strand(ctx, e, alpha)[p]
    cod = ctx[e, alpha][p]
    S = graded_ring(ctx)
    kk = coefficient_ring(S)
    G = grading_group(S)
    delta = sum(degree(x; check=false) for x in gens(S); init=zero(G))
    all_mons = all_monomials(ctx, alpha+e*delta)
    maps = FreeModuleHom{FreeMod{elem_type(kk)}, FreeMod{elem_type(kk)}, Nothing}[map_to_original_complex(simplified_fine_strand(ctx, ee.-e))[p] for ee in all_mons]
    #for (k, phi) in enumerate(maps)
      #@assert domain(phi) === codomain(canonical_projection(dom, k))
      #@assert codomain(phi) === codomain(canonical_projection(cod, k))
    #end
    return _direct_sum(kk, maps; domain=dom, codomain=cod)
  end
end

function simplified_strand_projection(
    ctx::NewToricCtx, e::Vector{Int}, alpha::FinGenAbGroupElem, p::Int
  )
  return simplified_strand_projection(ctx, first(e), alpha, p)
end

function simplified_strand_projection(
    ctx::NewToricCtx, e::Int, alpha::FinGenAbGroupElem, p::Int
  )
  return get!(ctx.simplified_strand_projectionss, (e, alpha, p)) do
    cod = simplified_strand(ctx, e, alpha)[p]
    dom = ctx[e, alpha][p]
    S = graded_ring(ctx)
    kk = coefficient_ring(S)
    G = grading_group(S)
    delta = sum(degree(x; check=false) for x in gens(S); init=zero(G))
    all_mons = all_monomials(ctx, alpha+e*delta)
    maps = [map_from_original_complex(simplified_fine_strand(ctx, ee.-e))[p] for ee in all_mons]
    #for (k, phi) in enumerate(maps)
      #@assert domain(phi) === codomain(canonical_projection(dom, k))
      #@assert codomain(phi) === codomain(canonical_projection(cod, k))
    #end
    return _direct_sum(kk, maps; domain=dom, codomain=cod)
  end
end

function cohomology_model(ctx::NewToricCtx, d::FinGenAbGroupElem)
  get!(ctx.cohomology_models, d) do
    simplified_strand(ctx, _minimal_exponent_vector(ctx, d), d)
  end
end

function cohomology_model_inclusion(ctx::NewToricCtx, d::FinGenAbGroupElem, i::Int)
  return simplified_strand_inclusion(ctx, _minimal_exponent_vector(ctx, d), i)
end

function cohomology_model_projection(ctx::NewToricCtx, d::FinGenAbGroupElem, i::Int)
  return simplified_strand_projection(ctx, _minimal_exponent_vector(ctx, d), i)
end

# return the minimal exponent `e`  such that the whole 
# cohomology in degree `m` is contained in the truncated ̌complex for `e`.
function _minimal_exponent_vector(ctx::NewToricCtx, m::FinGenAbGroupElem)
  # TODO: Make this use the already build cache for the support sets! 
  !isnothing(ctx.fixed_exponent_vector) && return ctx.fixed_exponent_vector::Int
  return get!(ctx.exp_vec_cache, m) do
    X = toric_variety(ctx)
    return maximum([optimal_k(ctx, i, m) for i in 0:dim(X)])
  end
end

function getindex(ctx::NewToricCtx, e0::Vector{Int}, e1::Vector{Int}, alpha::FinGenAbGroupElem)
  return getindex(ctx, first(e0), first(e1), alpha)
end

function getindex(ctx::NewToricCtx, e0::Int, e1::Int, alpha::FinGenAbGroupElem)
  return get!(ctx.strand_inclusions, (e0, e1, alpha)) do
    # TODO: This is not lazy!
    @assert e0 <= e1
    S = graded_ring(ctx)
    G = grading_group(S)
    @assert G === parent(alpha)
    delta = sum(degree(x; check=false) for x in gens(S); init=zero(G))
    all_mons_dom = all_monomials(ctx, alpha + e0*delta)
    all_mons_cod_inv = all_monomials_inv(ctx, alpha + e1*delta)
    index_mapping = Int[all_mons_cod_inv[ee.+(e1 - e0)] for ee in all_mons_dom]
    dom = ctx[e0, alpha]
    cod = ctx[e1, alpha]
    p = 0
    map_dict = Dict{Tuple{Int}, FreeModuleHom}()
    while p >= -dim(toric_variety(ctx)) && can_compute_index(dom, p) && can_compute_index(cod, p)
      img_gens = elem_type(cod[p])[]
      for (k, i) in enumerate(index_mapping)
        #pr = canonical_projection(dom[p], k)
        inc = canonical_injection(cod[p], i)
        #@assert codomain(pr) === domain(inc)
        img_gens = vcat(img_gens, images_of_generators(inc))
      end
      map_dict[(p,)] = hom(dom[p], cod[p], img_gens)
      p -= 1
    end
    return MorphismFromDict(dom, cod, map_dict)
  end
end

function sample_multiplication(ctx::NewToricCtx, e0::Vector{Int}, e1::Vector{Int})
  @assert all(i <= j for (i, j) in zip(e0, e1))
  return sample_multiplication(ctx, sum((i < 0) << (k-1) for (k, i) in enumerate(e0); init=0),
                               sum((i < 0) << (k-1) for (k, i) in enumerate(e1); init=0)
                              )
end

function sample_multiplication(ctx::NewToricCtx, e0::Int, e1::Int)
  return get!(ctx.fine_strand_morphisms, (e0, e1)) do
    dom = fine_strand(ctx, e0)
    cod = fine_strand(ctx, e1)
    delta = e0 - e1
    S = fine_graded_ring(ctx)
    #ee0 = Int[-(e0 >> (k-1))%2 for k in 1:ngens(S)]
    #ee1 = Int[-(e1 >> (k-1))%2 for k in 1:ngens(S)]
    #@show ee0, ee1
    delta_exp = Int[-(delta >> (k-1))%2 for k in 1:ngens(S)]
    #@show delta_exp
    mon = prod(x^(-k) for (k, x) in zip(delta_exp, gens(S)); init=one(S))
    #@show mon
    p = 0
    map_dict = Dict{Tuple{Int}, OFPModuleHom}()
    while p >= -dim(toric_variety(ctx)) && can_compute_index(dom, p) && can_compute_index(cod, p)
      inc = inclusion_map(dom)[p]
      pr = projection_map(cod)[p]
      imgs = [mon*inc(v) for v in gens(dom[p])]
      @assert all(degree(v; check=false) == degree(cod) for v in imgs)
      img_gens = elem_type(cod[p])[pr(v) for v in imgs]
      !is_zero(dom[p]) && @assert !is_zero(img_gens)
      map_dict[(p,)] = hom(dom[p], cod[p], img_gens)
      p -= 1
    end
    return MorphismFromDict(dom, cod, map_dict)
  end
end

function multiplication_map(
    ctx::NewToricCtx, 
    p::MPolyDecRingElem,
    e0::Vector{Int}, alpha::FinGenAbGroupElem, 
    j::Int
  )
  return multiplication_map(ctx, p, first(e0), alpha, j)
end

function multiplication_map(
    ctx::NewToricCtx, 
    p::MPolyDecRingElem,
    e0::Int, alpha::FinGenAbGroupElem, 
    j::Int
  )
  cache = get!(ctx.mult_map_cache, (e0, alpha, j)) do
    WeakKeyDict{typeof(p), Map}()
  end
  neg_res = get(cache, -p, nothing)
  !isnothing(neg_res) && return MapFromFunc(domain(neg_res), codomain(neg_res), v->-neg_res(v))
  return get!(cache, p) do
    beta = alpha + degree(p; check=false)
    dom_cplx = ctx[e0, alpha]
    cod_cplx = ctx[e0, beta]
    dom = dom_cplx[j]
    cod = cod_cplx[j]
    S = graded_ring(ctx)
    @assert S === parent(p)
    G = grading_group(S)
    delta = sum(degree(x; check=false) for x in gens(S); init=zero(G))
    all_mons_dom = all_monomials(ctx, alpha + e0*delta)
    all_mons_cod_inv = all_monomials_inv(ctx, beta + e0*delta)
    #img_gens = elem_type(cod)[zero(cod) for _ in 1:ngens(dom)]
    # TODO: Make this efficient by building the morphism directly 
    # from inplace operations on SRows.
    kk = coefficient_ring(S)
    img_gens = sparse_matrix(kk, 0, ngens(cod))
    for _ in 1:ngens(dom)
      push!(img_gens, sparse_row(kk))
    end
    #result = hom(dom, cod, elem_type(cod)[zero(cod) for _ in 1:ngens(dom)])
    for (k, ee) in enumerate(all_mons_dom)
      #pr = canonical_projection(dom, k)
      inc_dom = canonical_injection(dom, k)
      for (c, eee) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
        sample = sample_multiplication(ctx, ee.-e0, (ee + eee).-e0)[j]
        #@assert codomain(pr) === domain(sample)
        inc = canonical_injection(cod, all_mons_cod_inv[ee+eee])
        psi = compose(sample, inc)
        psi_imgs = images_of_generators(psi)
        for (ii, v) in enumerate(images_of_generators(inc_dom))
          jj, _ = only(coordinates(v))
          Hecke.add_scaled_row!(coordinates(psi_imgs[ii]), img_gens[jj], c)
        end
        @assert codomain(sample) === domain(inc)
        #result += compose(pr, c*compose(sample, inc))
      end
    end
    res = hom(dom, cod, elem_type(cod)[cod(v) for v in img_gens])
    #@assert res == result
    res
  end
end

function support_set(ctx::NewToricCtx, p::Int)
  return get!(ctx.support_sets, p) do 
    X = toric_variety(ctx)
    S = graded_ring(ctx)
    n = ngens(S)
    result = Int[]
    for i in 0:2^n-1
      s = fine_strand(ctx, i)
      r = ngens(s[-p])
      if !is_zero(p)
        r -= get!(ctx.fine_strand_map_ranks, (i, -p+1)) do
          rank(matrix(map(s, -p+1)))
        end
      end
      if p != dim(X)
        r -= get!(ctx.fine_strand_map_ranks, (i, -p)) do
          rank(matrix(map(s, -p)))
        end
      end
      !is_zero(r) && push!(result, i)
      # s = simplified_fine_strand(ctx, i)
      # !is_zero(s[-p]) && push!(result, i)
    end
    result
  end
end

function optimal_k(ctx::NewToricCtx, i::Int, alpha::FinGenAbGroupElem)
  return get!(ctx.optimal_ks, (i, alpha)) do
    X = toric_variety(ctx)
    k = 0
    S = cox_ring(X)
    G = grading_group(S)
    n = ngens(S)
    @assert 0 <= i <= dim(X) "index out of bounds"
    @assert parent(alpha) === G "degree does not belong to the grading group"
    Sigma_i = support_set(ctx, i)
    alpha_vec = elem_type(ZZ)[alpha[i] for i in 1:rank(G)]
    phi = map_from_torusinvariant_weil_divisor_group_to_class_group(X)
    A = transpose(matrix(phi))
    for I in Sigma_i
      #expanded = [(I >> j)%2 for j in 0:n-1]
      L_I = zero_matrix(ZZ, n, n)
      for j in 1:n
        L_I[j, j] = is_zero((I >> (j-1))%2) ? -1 : 1
      end
      D = vcat(L_I, A, -A)
      #b = vcat(elem_type(ZZ)[j in I ? -1 : 0 for j in 1:n], alpha_vec, -alpha_vec)
      b = elem_type(ZZ)[is_zero((I >> (j-1))%2) ? 0 : -1 for j in 1:n]
      #@show L_I
      #@show b
      #@show A
      #@show alpha_vec
      #P = polyhedron(QQ, D, b)
      @assert nrows(A) == length(alpha_vec)
      @assert ncols(L_I) == ncols(A)
      #return L_I, b, A, alpha_vec
      P = polyhedron((L_I, b), (A, alpha_vec))
      #@show is_bounded(P)
      #!is_bounded(P) && error("polyhedron not bounded")
      #@show lattice_points(P)
      for j in 1:n
        is_zero((I >> (j-1))%2) && continue
        l = elem_type(ZZ)[i == j ? 1 : 0 for i in 1:n]
        #lp = linear_program(P, l)
        #v1, _ = solve_lp(lp)
        lp = linear_program(P, -l)
        v, _ = solve_lp(lp)
        isnothing(v) && break # empty polyhedron
        is_infinite(v) && error("polyhedron not bounded")
        #@show j, v
        if v > k
          k = Int(floor(v))
        end
      end
    end
    return k
  end::Int
end

relative_dimension(ctx::NewToricCtx) = dim(toric_variety(ctx))
target_ring(ctx::NewToricCtx) = coefficient_ring(cox_ring(toric_variety(ctx)))
