mutable struct CSSPage#{GradedRingType, CoeffRingType}
  parent::Any # The spectral sequence; type declaration below
  page_number::Int
  entries::Dict{Tuple{Int, Int}, ModuleFP}
  maps::Dict{Tuple{Int, Int}, ModuleFPHom}
  lifted_kernel_generators::Dict{Tuple{Int, Int}, Vector{Tuple{Vector{Int}, Vector{Tuple{Int, <:FreeModElem}}}}}

  function CSSPage(css::Any, k::Int)
    return new(css, k)
  end
end

mutable struct CohomologySpectralSequence{GradedRingType, CoeffRingType}
  S::GradedRingType
  A::CoeffRingType
  graded_complex::AbsHyperComplex
  pages::Dict{Int, CSSPage}
  pfctx::PushForwardCtx

  function CohomologySpectralSequence(S::MPolyRing, comp::AbsHyperComplex)
    @assert isone(dim(comp)) "complex must be 1-dimensional"
    A = coefficient_ring(S)
    return new{typeof(S), typeof(A)}(S, A, comp, Dict{Int, CSSPage}())
  end
end

function pushforward_ctx(css::CohomologySpectralSequence)
  if !isdefined(css, :pfctx)
    css.pfctx = PushForwardCtx(graded_ring(css))
  end
  return css.pfctx
end

function pushforward_ctx(cssp::CSSPage)
  return pushforward_ctx(spectral_sequence(cssp))
end

function graded_ring(cssp::CSSPage)
  return graded_ring(spectral_sequence(cssp))
end

function graded_ring(css::CohomologySpectralSequence)
  return css.S
end

# get the k-th page of the spectral sequence
function getindex(css::CohomologySpectralSequence, k::Int)
  k > 0 || error("index out of bounds")
  return get!(css.pages, k) do
    CSSPage(css, k)
  end
end

function getindex(css::CohomologySpectralSequence, k::Int, i::Int, j::Int)
  return (css[k])[i, j]
end

function entries(cssp::CSSPage)
  if !isdefined(cssp, :entries)
    cssp.entries = Dict{Tuple{Int, Int}, ModuleFP}()
  end
  return cssp.entries
end

function getindex(cssp::CSSPage, i::Int, j::Int)
  @assert can_compute_index(graded_complex(spectral_sequence(cssp)), i) "index out of bounds"
  @assert j<=0 "index out of bounds"
  return get!(entries(cssp), (i, j)) do
    produce_entry(cssp, i, j)
  end
end
  
function lifted_kernel_generators(cssp::CSSPage)
  if !isdefined(cssp, :lifted_kernel_generators)
    cssp.lifted_kernel_generators = Dict{Tuple{Int, Int}, 
                                         Vector{Tuple{Vector{Int}, Vector{Tuple{Int, <:FreeModElem}}}}}()
  end
  return cssp.lifted_kernel_generators
end

function lifted_kernel_generators(cssp::CSSPage, i::Int, j::Int)
  return get!(lifted_kernel_generators(cssp), (i, j)) do
    produce_lifted_kernel_generators(cssp, i, j)
  end
end

function produce_entry(cssp::CSSPage, i::Int, j::Int)
  is_one(page_number(cssp)) && return produce_entry_on_initial_page(cssp, i, j)
  css = spectral_sequence(cssp)
  p = page_number(cssp)
  previous_page = css[p-1]
  H = previous_page[i, j]
  Z, inc_Z = can_compute_map(previous_page, i, j) ? kernel(map(previous_page, i, j)) : sub(H, gens(H))
  B, inc_B = can_compute_map(previous_page, i+p-1, j-p+2) ? image(map(previous_page, i+p-1, j-p+2)) : sub(H, elem_type(H)[])
  return quo(Z, B)[1]
end

function produce_entry_on_initial_page(cssp::CSSPage, i::Int, j::Int)
  c = graded_complex(cssp)
  F = c[i]
  degs = degrees_of_generators(F)
  ctx = pushforward_ctx(cssp)
  kk = coefficient_ring(graded_ring(cssp))
  summands = FreeMod{elem_type(kk)}[cohomology_model(ctx, -d)[j] for d in degs]
  is_empty(summands) && return FreeMod(kk, 0)
  return direct_sum(summands)[1]
end

function can_compute_index(cssp::CSSPage, i::Int, j::Int)
  can_compute_index(graded_complex(cssp), i) || return false
  j <= 0 || return false
  S = graded_ring(cssp)
  -j <= ngens(S) - rank(grading_group(S)) || return false
  return true
end

function can_compute_map(cssp::CSSPage, i::Int, j::Int)
  p = page_number(cssp)
  return can_compute_index(cssp, i, j) && can_compute_index(cssp, i-p, j+p-1)
end

function produce_entry_on_second_page(cssp::CSSPage, i::Int, j::Int)
  css = spectral_sequence(cssp)
  p = page_number(cssp)
  #first_page = css[1]
  previous_page = css[p-1]
  H = previous_page[i, j]
  Z, inc_Z = can_compute_map(previous_page, i, j) ? kernel(map(previous_page, i, j)) : sub(H, gens(H))
  B, inc_B = can_compute_map(previous_page, i+p-1, j-p+2) ? image(map(previous_page, i+p-1, j-p+2)) : sub(H, elem_type(H)[])
  return quo(Z, B)[1]
  return SubquoModule(F, ambient_representatives_generators(Z), 
                      ambient_representatives_generators(B))
end


function maps(cssp::CSSPage)
  if !isdefined(cssp, :maps)
    cssp.maps = Dict{Tuple{Int, Int}, ModuleFPHom}()
  end
  return cssp.maps
end

function map(cssp::CSSPage, i::Int, j::Int)
  return get!(maps(cssp), (i, j)) do
    produce_map(cssp, i, j)
  end
end

function produce_lifted_kernel_generators_on_initial_page(cssp::CSSPage, i::Int, j::Int)
  S = graded_ring(cssp)
  ctx = pushforward_ctx(cssp)
  cplx = graded_complex(cssp)
  orig_map = map(cplx, i+1)
  dom_degs = degrees_of_generators(cplx[i+1])
  dom = cssp[i+1, j]
  result = Tuple{Vector{Int}, Vector{Tuple{Int, FreeModElem}}}[]
  is_zero(dom) && return result
  prs = canonical_projections(dom)
  rngs = get_attribute(dom, :ranges)::Vector{UnitRange{Int}}
  for (v_num, v) in enumerate(gens(dom))
    @show v_num
    buckets = Tuple{Int, FreeModElem}[]
    e = Int[]
    ind, _ = only(coordinates(v))
    @show ind
    k = findfirst(ind in r for r in rngs)
    #for (k, pr) in enumerate(prs)
      vv = prs[k](v)
      @assert !is_zero(vv)
      #is_zero(vv) && continue # TODO: Make this quicker!
      d = -dom_degs[k]
      vv2 = cohomology_model_inclusion(ctx, d, j)(vv)
      @assert !is_zero(vv2)
      e = _minimal_exponent_vector(ctx, d)
      @show e
      for (l, p) in coordinates(images_of_generators(orig_map)[k])
        @show l, p
        phi = multiplication_map(ctx, p, e, d, j)
        vvv = phi(vv2)
        @assert !is_zero(vvv)
        ind = findfirst(ll == l for (ll, _) in buckets)
        @show vvv
        if isnothing(ind)
          push!(buckets, (l, vvv))
        else
          error("this should not happen")
          _, www = buckets[ind]
          buckets[ind] = (l, www + vvv)
        end
      #end
      #break
    end
    @assert !is_empty(e)
    push!(result, (e, buckets))
  end
  return result
end

function produce_lifted_kernel_generators(cssp::CSSPage, i::Int, j::Int)
  @show "production of lifted kernel generators on page $(page_number(cssp)) at place $i, $j"
  @assert can_compute_index(graded_complex(cssp), i) "index out of bounds"
  @assert j <= 0 "index out of bounds"

  # Let `p` be the number of the page `cssp`. We have the 
  # incoming boundary map
  #
  #        cssp[i, j]
  #                  ↖ ∂ₚ
  #                    cssp[i0, j0]
  #
  # with i0 = i + p, j0 = j - p + 1. This map posesses a lift ψ
  #                      ̌Cᵢⱼ
  #                   ↙       ↖ ψ
  #        cssp[i, j]            Zₚ = ker ∂ₚ₋₁ ⊂ ̌Cᵢ₀,ⱼ₀
  #                  ↖ ∂ₚ      ↙
  #                    cssp[i0, j0]
  #
  # which takes the generators v of the kernel Zₚ to elements 
  # ψ(v) in the respective entry of the ̌Cech double complex. 
  #
  # This procedure is to produce those elements ψ(v).

  S = graded_ring(cssp)
  n = ngens(S)
  ctx = pushforward_ctx(cssp)

  # In case of a trivial incoming map, return the empty list
  if !(can_compute_index(graded_complex(cssp), i+page_number(cssp)) || 
       j - page_number(cssp) + 1 <= - n + rank(grading_group(S)))
    return [0 for _ in 1:n], Tuple{Int, FreeModElem}[]
  end
  
  # On the first page special initialization is required. 
  if is_one(page_number(cssp))
    return produce_lifted_kernel_generators_on_initial_page(cssp, i, j)
  end
  
  p = page_number(cssp)
  css = spectral_sequence(cssp)

  # We do things by induction, lifting suitable linear combinations 
  # of the `lifted_kernel_generators` on the previous page.
  prev_page = css[p-1]
  # the indices of the domain
  i0 = i + p 
  j0 = j - p + 1

  # the indices of the codomain of the map on the previous page
  ii = i + 1
  jj = j - 1

  @show page_number(prev_page), ii, jj
  v = lifted_kernel_generators(prev_page, ii, jj) # lifts of the kernel generators on the previous page

  # The generators of the kernel Zₚ = ker ∂ₚ₋₁ are extracted from this page.
  u_next = repres.(gens(cssp[i0, j0]))

  orig_cplx = graded_complex(cssp)
  # The degrees of the generators in the original graded complex 
  # determine the exponents of the numerators in the Cech complex 
  # which are required to properly represent the involved cohomology classes. 
  dom_degs = degrees_of_generators(orig_cplx[i0])
  inter_degs = degrees_of_generators(orig_cplx[ii])

  # Initialize the result.
  # For every element ψ(v) we store a pair (e, comps) where `e` is the 
  # exponent vector for the Cech-complex and `comps` is something like 
  # a "sparse block vector", meaning a list of pairs `(k, v_k)` where 
  # `k` indicates the block for the `k`-th generator of the i-th entry 
  # of the original complex and `v_k` is the corresponding element 
  # in the Cech-complex for that block (for exponents e).
  #
  # This is rather complicated, but it prevents us from having to 
  # spell out the full total complex behind this spectral sequence. 
  new_lifted_kernel_gens = Tuple{Vector{Int}, Vector{Tuple{Int, <:FreeModElem}}}[]

  for (u_num, u) in enumerate(u_next)
    @show u_num
    @show u
    # We write `u` in the generators of `Zₚ₋₁` on which the previous 
    # map ∂ₚ₋₁ is defined.
    u_coords = p == 2 ? coordinates(u) : coordinates(u, prev_page[i0, j0])
    @show u_coords

    # To bring everything to the same denominator in the Cech complex 
    # we need to find a common exponent.
    exps = Vector{Int}[v[i][1] for (i, c) in u_coords]
    sup_exp = Int[maximum([e[i] for e in exps]; init=0) for i in 1:n]
    @show exps
    @show sup_exp

    buckets = Tuple{Int, FreeModElem}[]
    # ∂ₚ₋₁(u) = 0, so the linear combination for `u` in the generators of `Zₚ₋₁` 
    # applied to the `lifted_kernel_generators` for ψₚ₋₁ gives us elements 
    # which can be lifted further
    for (l, c) in u_coords
      @show l, c
      e, g = v[l]
      @show e, g
      # bring everything to the same denominators
      if e != sup_exp
        g = Tuple{Int, <:FreeModElem}[(j, ctx[e, sup_exp, -inter_degs[j]][jj](v)) for (j, v) in g]
      end
      for (j, v) in g
        k = findfirst(k==j for (k, _) in buckets)
        if isnothing(k)
          push!(buckets, (j, c*v))
        else
          _, w = buckets[k]
          buckets[k] = (j, w + c*v)
        end
      end
    end

    # Now the buckets hold lifts of the kernel generators on page `p`. 
    # We need to lift them once more through the double complex.
    lifted_buckets = Tuple{Int, FreeModElem}[]
    for (i, v) in buckets
      @show "lifting $i, $v"
      is_zero(v) && continue
      @show sup_exp
      @show _minimal_exponent_vector(ctx, -inter_degs[i])
      @show -inter_degs[i]
      @show ctx[sup_exp, -inter_degs[i]]
      @show jj+1
      cech_map = map(ctx[sup_exp, -inter_degs[i]], jj+1)
      @show matrix(cech_map)
      # We would like to lift via the Cech map, but modulo 
      # the image of previous differentials. 
      # In principal, this can be done with the following code, 
      # but that computes a Groebner basis every single time. 
      e0 = _minimal_exponent_vector(ctx, -inter_degs[i])
      H = prev_page[ii, jj]
      F = ambient_free_module(H)
      pr = canonical_projection(F, i)
      coh_inc = cohomology_model_inclusion(ctx, -inter_degs[i], jj)
      mult = ctx[e0, sup_exp, -inter_degs[i]][jj]
      I, inc = sub(parent(v), [mult(coh_inc(pr(v))) for v in relations(H)])
      @show gens(I)
      M, pro = quo(parent(v), I)
      psi = compose(cech_map, pro)
      @show v
      @show pro(v)
      vv = preimage(psi, pro(v))
      @show vv
      @show cech_map(vv) - v
      # Instead, we would like to make use of the rather simple form 
      # of the Cech map (whose image is a primitive direct summand 
      # of its codomain) and lift with that via some brute force. 
      push!(lifted_buckets, (i, vv))
    end
    
    # The lifted elements still need to be mapped.
    mapped_buckets = Tuple{Int, FreeModElem}[]
    orig_map = map(graded_complex(cssp), ii)
    for (k, v) in lifted_buckets
      img_gen = images_of_generators(orig_map)[k]
      for (l, p) in coordinates(img_gen)
        phi = multiplication_map(ctx, p, sup_exp, -inter_degs[k], jj+1)
        vv = phi(v)
        is_zero(vv) && continue
        ind = findfirst(ll == l for (ll, _) in mapped_buckets)
        if isnothing(ind)
          push!(mapped_buckets, (l, vv))
        else
          _, ww = mapped_buckets[ind]
          mapped_buckets[ind] = l, ww + vv
        end
      end
    end

    # store the result
    push!(new_lifted_kernel_gens, (sup_exp, mapped_buckets))
  end
  for (v, (e, b)) in zip(u_next, new_lifted_kernel_gens)
    println("$(coordinates(v)) -> $e : $b")
  end
  return new_lifted_kernel_gens
end


function produce_map(cssp::CSSPage, i::Int, j::Int)
  @assert can_compute_index(graded_complex(cssp), i) "index out of bounds"
  @assert can_compute_index(graded_complex(cssp), i-page_number(cssp)) "index out of bounds"
  @assert j <= 0 "index out of bounds"
  @assert j + page_number(cssp) - 1 <= 0 "index out of bounds"
  #is_one(page_number(cssp)) && return produce_map_on_initial_page(cssp, i, j)
  #page_number(cssp) == 2 && return produce_map_on_second_page(cssp, i, j)

  p = page_number(cssp)
  css = spectral_sequence(cssp)
  ctx = pushforward_ctx(cssp)
  
  # the indices of the codomain
  i1 = i - p 
  j1 = j + p - 1

  new_lifted_kernel_gens = lifted_kernel_generators(cssp, i1, j1)

  # assemble the actual induced morphism
  dom = cssp[i, j]
  cod = cssp[i1, j1]
  F = ambient_free_module(cod)
  img_gens = elem_type(F)[]
  cod_degs = degrees_of_generators(graded_complex(cssp)[i1])
  for (e, buckets) in new_lifted_kernel_gens
    img_gen = zero(F)
    for (k, v) in buckets
      is_zero(v) && continue
      d1 = -cod_degs[k]
      e1 = _minimal_exponent_vector(ctx, d1)
      pr = ctx[e, e1, d1][j1]
      vv = pr(v)
      is_zero(vv) && continue
      coh_pr = cohomology_model_projection(ctx, d1, j1)
      img_gen += canonical_injection(F, k)(coh_pr(vv))
      #img_gen += cohomology_model_projection(ctx, d1, j1)(pr(v))
    end
    push!(img_gens, img_gen)
  end
  return hom(dom, cod, elem_type(cod)[cod(v) for v in img_gens]) #; check=false)

  #=
  css = spectral_sequence(cssp)
  p = page_number(cssp)
  cmplx = graded_complex(css)
  orig_dom = cmplx[i]
  orig_inter = cmplx[i-p+1]
  orig_map2 = map(cmplx, i-p+1)
  p1 = css[1]
  ctx = pushforward_ctx(css)
  S = graded_ring(ctx)
  dom = cssp[i, j]
  inter = p1[i-p+1, j]
  cod = cssp[i-p, j+p-1]
  if is_zero(dom) || is_zero(cod)
    return hom(dom, cod, elem_type(cod)[zero(cod) for _ in 1:ngens(dom)])
  end
  img_gens = elem_type(cod)[]
  map_on_previous_page = map(css[p-1], i, j)
  inter_injs = canonical_injections(inter)
  inter_coh_incs = [cohomology_model_inclusion(ctx, d1, j) for d1 in degrees_of_generators(orig_inter)]
  # We need to do a knights move 
  #
  #     p1[i-2, j+1] <-----  p1[i-1, j+1]
  #
  #                          p1[i-1, j]   <------ p1[i, j] 
  #
  # with a double complex of free modules behind it:
  #
  #         D        <--ψ---       C
  #                                ↓∂
  #                                B      <---ϕ----   A
  #
  # The horizontal maps are induced from those in the `graded_complex`.
  # The vertical maps are direct sums of the standard Cech-complexes 
  # for either of the summands in the `graded_complex`. 
  # We exploit this block structure and construct the induced map 
  # term by term.
  for (ii, v) in enumerate(gens(dom))
    @show v
    @show is_zero(v)
    if is_zero(v)
      push!(img_gens, zero(cod))
      continue
    end
    v1 = repres(v)::FreeModElem
    v1_blocks = typeof(v1)[pr(v1) for pr in canonical_projections(parent(v1))]
    exps = Vector{Int}[_minimal_exponent_vector(ctx, -degrees_of_generators(orig_dom)[k]) for (k, v1_block) in enumerate(v1_blocks) if !is_zero(v1_block)]
    @show exps
    sup_exp = [maximum([e[i] for e in exps]; init=0) for i in 1:ngens(S)]
    @show sup_exp

    @assert is_zero(map_on_previous_page(v1))
    @show v1
    F = parent(v1)
    # to hold the block-pieces of the image of v in B.
    img_buckets = [zero(ctx[_minimal_exponent_vector(ctx, -d1), -d1][j]) for d1 in degrees_of_generators(orig_inter)]
    # we can lift v1 blockwise
    for (k, v_block) in enumerate(v1_blocks)
      @show k
      d0 = -degree(orig_dom[k]) # the original shift 
      is_zero(v_block) && continue
      e0 = _minimal_exponent_vector(ctx, d0)
      coh_mod = ctx[e0, d0]
      inc = cohomology_model_inclusion(ctx, d0, j)
      v3 = inc(v_block)
      @show v3
      @show j 
      @assert codomain(inc) === coh_mod[j]
      for (l, p) in coordinates(images_of_generators(orig_map1)[k])
        @show l
        d1 = d0 + degree(p; check=false)
        @assert d1 == -degrees_of_generators(orig_inter)[l]
        e1 = _minimal_exponent_vector(ctx, d1)
        phi = multiplication_map(ctx, p, e0, d0, j)
        w1 = phi(v3)
        is_zero(w1) && continue
        @show w1
        str_pr = ctx[e0, e1, d1]
        @show codomain(str_pr[j]) 
        @show parent(img_buckets[l])
        @assert codomain(str_pr[j]) === parent(img_buckets[l])
        w2 = str_pr[j](w1)
        cech_map = map(ctx[e1, d1], j+1)
        @assert codomain(cech_map) === parent(w2)
        @show w2
        is_zero(w2) && continue
        img_buckets[l] += w2
      end
    end
    img_gen = zero(ambient_free_module(cod))
    for (l, w, dd1) in zip(1:length(img_buckets), img_buckets, degrees_of_generators(orig_inter))
      is_zero(w) && continue
      d1 = -dd1
      #is_zero(w) && continue
      e1 = _minimal_exponent_vector(ctx, d1)
      cech_map = map(ctx[e1, d1], j+1)
      @assert codomain(cech_map) === parent(w)
      ww = preimage(cech_map, w)
      for (r, q) in coordinates(images_of_generators(orig_map2)[l])
        d2 = d1 + degree(q; check=false)
        ext_cod_strand = ctx[e1, d2][j+1]
        mult_map = multiplication_map(ctx, q, e1, d1, j+1)
        @assert domain(mult_map) === domain(cech_map)
        @assert parent(ww) === domain(mult_map)
        ww2 = mult_map(ww)
        is_zero(ww2) && continue
        e2 = _minimal_exponent_vector(ctx, d2)
        ww3 = ctx[e1, e2, d2][j+1](ww2)
        is_zero(ww3) && continue
        ww4 = cohomology_model_projection(ctx, d2, j+1)(ww3)
        is_zero(ww4) && continue
        ww5 = canonical_injection(ambient_free_module(cod), r)(ww4)
        @show ww5
        @show i-2, j+1
        if i > 2
          bb = map(p1, i-2, j+1)
          @show ambient_free_module(domain(bb)) 
          @show parent(ww5)
          @assert ambient_free_module(domain(bb)) === parent(ww5)
          @show is_zero(bb(ww5))
        end
        img_gen += ww5
      end
    end
    push!(img_gens, cod(img_gen))
  end
  return hom(dom, cod, img_gens)
  =#
end

relations(F::FreeMod) = elem_type(F)[]

function multiplication_map(
    ctx::PushForwardCtx, 
    p::MPolyDecRingElem,
    e0::Vector{Int}, d0::FinGenAbGroupElem, 
    j::Int
  )
  return get!(ctx.mult_map_cache, (p, e0, d0, j)) do
    d1 = d0 + degree(p; check=false)
    dom_cplx = ctx[e0, d0]
    cod_cplx = ctx[e0, d1]
    dom = dom_cplx[j]
    cod = cod_cplx[j]
    dom_strand_inc = inclusion_map(dom_cplx)[j]
    cod_strand_pr = projection_map(cod_cplx)[j]
    img_gens = elem_type(cod)[cod_strand_pr(p*dom_strand_inc(v)) for v in gens(dom)]
    hom(dom, cod, img_gens)
  end
end

function produce_map_on_initial_page(cssp::CSSPage, i::Int, j::Int)
  dom = cssp[i, j]::FreeMod # the cohomology model
  cod = cssp[i-1, j]::FreeMod # the cohomology model
  c = graded_complex(cssp) # the original complex of graded modules
  orig_dom = c[i]::FreeMod # the graded module for this column
  orig_cod = c[i-1]::FreeMod # the graded module for the column of the codomain
  orig_map = map(c, i)::FreeModuleHom # the original map between the graded modules
  # initialize the result
  result = hom(dom, cod, elem_type(cod)[zero(cod) for _ in 1:ngens(dom)])
  # `ctx` holds precomputed truncated Cech-complexes, their inclusions, and their strands
  ctx = pushforward_ctx(cssp)
  # the cohomology model is a direct sum of modules for each generator of `orig_dom`.
  dom_proj = canonical_projections(dom) # projections for the direct sum
  cod_inj = canonical_injections(cod)   # inclusions for the direct sum

  # We assemble the induced maps on the cohomology model.
  # It comes in blocks for the generators of `orig_dom` and `orig_cod`.
  for (k, v) in enumerate(images_of_generators(orig_map))
    d0 = -degree(orig_dom[k]) # the relevant degree in the domain for this block
    e0 = _minimal_exponent_vector(ctx, d0) # determines the exponent for the 
                                           # truncated Cech-complex for the cohomology 
                                           # model in the domain for this block
    dom_complex = ctx[e0, d0] # the said truncated Cech-complex
    dom_summand = dom_complex[j] 
    pr_k = dom_proj[k] # from the cohomology module to the summand
    @assert domain(pr_k) === dom
    dom_coh_complex = cohomology_model(ctx, d0) # the simplified version of dom_complex
    dom_coh = dom_coh_complex[j]
    @assert codomain(pr_k) === dom_coh

    h_inj = cohomology_model_inclusion(ctx, d0, j) # from the summand to the minimal 
                                                   # truncated Cech complex
    @assert codomain(pr_k) === domain(h_inj)
    dom_strand = codomain(h_inj)

    dom_str_inc = inclusion_map(dom_complex)[j]
    @assert domain(dom_str_inc) === dom_strand

    for (l, p) in coordinates(v)
      # On the `(k, l)`-th block we have the map which is induced by multiplication 
      # with `p` as `-⋅p : S¹ → S¹`. 
      # We extract this map on the strands of the truncated Cech-complex for the 
      # exponent vector `e0` (for degree `-d0`). This will then be projected 
      # to the minimal degree for the `l`-th generator and subsequently to the 
      # cohomology model. 
      inc_l = cod_inj[l] # inclusion of the direct summand for the `l`-th generator
                         # of the codomain
      @assert codomain(inc_l) === cod

      # Assemble the morphism given by multiplication with p
      #
      #   x ↦ p⋅x : dom_strand → ext_cod_strand
      #
      deg_p = degree(p; check=false)
      d1 = d0 + deg_p
      @assert d1 == -degree(orig_cod[l])
      h_pr = cohomology_model_projection(ctx, d1, j)
      cod_coh = cohomology_model(ctx, d1)[j]
      @assert cod_coh === codomain(h_pr)
      @assert domain(inc_l) === codomain(h_pr)

      e1 = _minimal_exponent_vector(ctx, d1)
      cod_cplx_ext = ctx[e0, d1]
      cod_cplx = ctx[e1, d1]
      strand_reduction = ctx[e0, e1, d1]
      @assert domain(strand_reduction) === cod_cplx_ext
      @assert codomain(strand_reduction) === cod_cplx


      cod_str_pr = projection_map(cod_cplx_ext)[j] # the map from the graded module to the extended strand
      ext_cod_strand = cod_cplx_ext[j]
      @assert codomain(cod_str_pr) === ext_cod_strand

      #img_gens = elem_type(ext_cod_strand)[cod_str_pr(p*dom_str_inc(v)) for v in gens(dom_strand)]
      #induced_map = hom(dom_strand, ext_cod_strand, img_gens)
      induced_map = multiplication_map(ctx, p, e0, d0, j)
      @assert domain(induced_map) === dom_strand
      @assert codomain(induced_map) === ext_cod_strand

      h_pr = cohomology_model_projection(ctx, d1, j)
      red = strand_reduction[j]
      block_img_gens = elem_type(cod_coh)[h_pr(
                                            red(
                                              induced_map(
                                                h_inj(v)
                                          ))) for v in gens(dom_coh)]
      result += compose(compose(dom_proj[k], hom(dom_coh, cod_coh, block_img_gens)), cod_inj[l])
    end
  end
  return result
end

function produce_map_on_second_page(cssp::CSSPage, i::Int, j::Int)
  css = spectral_sequence(cssp)
  cmplx = graded_complex(css)
  orig_dom = cmplx[i]
  orig_inter = cmplx[i-1]
  orig_map1 = map(cmplx, i)
  orig_map2 = map(cmplx, i-1)
  p1 = css[1]
  ctx = pushforward_ctx(css)
  dom = cssp[i, j]
  inter = p1[i-1, j]
  cod = cssp[i-2, j+1]
  if is_zero(dom) || is_zero(cod)
    return hom(dom, cod, elem_type(cod)[zero(cod) for _ in 1:ngens(dom)])
  end
  img_gens = elem_type(cod)[]
  map_on_first_page = map(p1, i, j)
  inter_injs = canonical_injections(inter)
  inter_coh_incs = [cohomology_model_inclusion(ctx, d1, j) for d1 in degrees_of_generators(orig_inter)]
  # We need to do a knights move 
  #
  #     p1[i-2, j+1] <-----  p1[i-1, j+1]
  #
  #                          p1[i-1, j]   <------ p1[i, j] 
  #
  # with a double complex of free modules behind it:
  #
  #         D        <--ψ---       C
  #                                ↓∂
  #                                B      <---ϕ----   A
  #
  # The horizontal maps are induced from those in the `graded_complex`.
  # The vertical maps are direct sums of the standard Cech-complexes 
  # for either of the summands in the `graded_complex`. 
  # We exploit this block structure and construct the induced map 
  # term by term.
  for (ii, v) in enumerate(gens(dom))
    @show v
    @show is_zero(v)
    if is_zero(v)
      push!(img_gens, zero(cod))
      continue
    end
    v1 = repres(v)::FreeModElem
    @assert is_zero(map_on_first_page(v1))
    @show v1
    F = parent(v1)
    # to hold the block-pieces of the image of v in B.
    img_buckets = [zero(ctx[_minimal_exponent_vector(ctx, -d1), -d1][j]) for d1 in degrees_of_generators(orig_inter)]
    # we can lift v1 blockwise
    for (k, pr) in enumerate(canonical_projections(F))
      @show k
      d0 = -degree(orig_dom[k]) # the original shift 
      v2 = pr(v1)
      @show v2
      is_zero(v2) && continue
      e0 = _minimal_exponent_vector(ctx, d0)
      coh_mod = ctx[e0, d0]
      inc = cohomology_model_inclusion(ctx, d0, j)
      @assert codomain(pr) === domain(inc)
      v3 = inc(v2)
      @show v3
      @show j 
      @assert codomain(inc) === coh_mod[j]
      for (l, p) in coordinates(images_of_generators(orig_map1)[k])
        @show l
        d1 = d0 + degree(p; check=false)
        @assert d1 == -degrees_of_generators(orig_inter)[l]
        e1 = _minimal_exponent_vector(ctx, d1)
        phi = multiplication_map(ctx, p, e0, d0, j)
        w1 = phi(v3)
        is_zero(w1) && continue
        @show w1
        str_pr = ctx[e0, e1, d1]
        @show codomain(str_pr[j]) 
        @show parent(img_buckets[l])
        @assert codomain(str_pr[j]) === parent(img_buckets[l])
        w2 = str_pr[j](w1)
        cech_map = map(ctx[e1, d1], j+1)
        @assert codomain(cech_map) === parent(w2)
        @show w2
        is_zero(w2) && continue
        img_buckets[l] += w2
      end
    end
    img_gen = zero(ambient_free_module(cod))
    for (l, w, dd1) in zip(1:length(img_buckets), img_buckets, degrees_of_generators(orig_inter))
      is_zero(w) && continue
      d1 = -dd1
      #is_zero(w) && continue
      e1 = _minimal_exponent_vector(ctx, d1)
      cech_map = map(ctx[e1, d1], j+1)
      @assert codomain(cech_map) === parent(w)
      ww = preimage(cech_map, w)
      for (r, q) in coordinates(images_of_generators(orig_map2)[l])
        d2 = d1 + degree(q; check=false)
        ext_cod_strand = ctx[e1, d2][j+1]
        mult_map = multiplication_map(ctx, q, e1, d1, j+1)
        @assert domain(mult_map) === domain(cech_map)
        @assert parent(ww) === domain(mult_map)
        ww2 = mult_map(ww)
        is_zero(ww2) && continue
        e2 = _minimal_exponent_vector(ctx, d2)
        ww3 = ctx[e1, e2, d2][j+1](ww2)
        is_zero(ww3) && continue
        ww4 = cohomology_model_projection(ctx, d2, j+1)(ww3)
        is_zero(ww4) && continue
        ww5 = canonical_injection(ambient_free_module(cod), r)(ww4)
        @show ww5
        @show i-2, j+1
        if i > 2
          bb = map(p1, i-2, j+1)
          @show ambient_free_module(domain(bb)) 
          @show parent(ww5)
          @assert ambient_free_module(domain(bb)) === parent(ww5)
          @show is_zero(bb(ww5))
        end
        img_gen += ww5
      end
    end
    push!(img_gens, cod(img_gen))
  end
  return hom(dom, cod, img_gens)
end


function stable_index(css::CohomologySpectralSequence, i::Int, j::Int)
  # return a page number k from which onwards the cohomology 
  # for the (i, j)-th entry has stabilized.
  error("not implemented")
end

function spectral_sequence(cssp::CSSPage)
  return cssp.parent::CohomologySpectralSequence # TODO: Add proper type assertion
end

function graded_complex(css::CohomologySpectralSequence) 
  return css.graded_complex
end

function graded_complex(cssp::CSSPage)
  return graded_complex(spectral_sequence(cssp))
end

function page_number(cssp::CSSPage)
  return cssp.page_number
end

function Base.show(io::IO, css::CohomologySpectralSequence)
  print(io, "spectral sequence for the cohomology of graded complex with Betti table\n")
  print(io, betti_table(graded_complex(css)))
end

function Base.show(io::IO, cssp::CSSPage)
  print(io, "page $(page_number(cssp)) of $(spectral_sequence(cssp))")
end
