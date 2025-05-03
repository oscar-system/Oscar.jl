mutable struct CSSPage#{GradedRingType, CoeffRingType}
  parent::Any # The spectral sequence; type declaration below
  page_number::Int
  entries::Dict{Tuple{Int, Int}, ModuleFP}
  maps::Dict{Tuple{Int, Int}, ModuleFPHom}

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

function produce_entry(cssp::CSSPage, i::Int, j::Int)
  is_one(page_number(cssp)) && return produce_entry_on_initial_page(cssp, i, j)
  page_number(cssp) == 2 && return produce_entry_on_second_page(cssp, i, j)
  error("not implemented")
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
  return can_compute_index(cssp, i, j) && can_compute_index(cssp, i-p, j-p+1)
end

function produce_entry_on_second_page(cssp::CSSPage, i::Int, j::Int)
  css = spectral_sequence(cssp)
  first_page = css[1]
  F = first_page[i, j]
  Z, inc_Z = can_compute_map(first_page, i, j) ? kernel(map(first_page, i, j)) : sub(F, gens(F))
  B, inc_B = can_compute_map(first_page, i+1, j) ? image(map(first_page, i+1, j)) : sub(F, elem_type(F)[])
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

function produce_map(cssp::CSSPage, i::Int, j::Int)
  @assert has_index(graded_complex(cssp), i) "index out of bounds"
  @assert has_index(graded_complex(cssp), i-page_number(cssp)) "index out of bounds"
  @assert j <= 0 "index out of bounds"
  @assert j + page_number(cssp) - 1 <= 0 "index out of bounds"
  is_one(page_number(cssp)) && return produce_map_on_initial_page(cssp, i, j)
  page_number(cssp) == 2 && return produce_map_on_second_page(cssp, i, j)
  error("not implemented")
end

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
  for (i, v) in enumerate(gens(dom))
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
    img_gen = zero(cod)
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
        ww5 = cod(ww4)
        img_gen += ww5
      end
    end
    push!(img_gens, img_gen)
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
