########################################################################
# Data structures and methods for the lazy computation of spectral 
# sequences for ̌Cech cohomology. 
########################################################################

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
  pfctx::Union{PushForwardCtx, ToricCtx, ToricCtxWithParams}

@doc raw"""
    CohomologySpectralSequence(S::MPolyRing, comp::AbsHyperComplex)

Suppose `comp` is a left-bounded complex of graded `FreeMod`s over a (multi-)graded 
polynomial ring ``S = R[xᵢⱼ : i = 1,…,m, j = 0,…,rᵢ]``. We assume 
that ``S`` is the homogeneous coordinate ring of a product of projective 
spaces ``ℙ = ℙʳ¹ × … × ℙʳᵐ → Spec R`` over ``R``, i.e. the grading group 
is ``ℤᵐ = ⨁ ᵢ₌₁ᵐ ℤ ⋅eᵢ`` with ``deg(xᵢⱼ) = eᵢ``.

Then `comp` represents a complex of coherent sheaves ``ℱ *`` on ``ℙ`` and 
there is a spectral sequence 
```math
    Eᵢⱼ¹ = Rʲπ_* (ℱ ⁱ) ⇒ Rⁱ⁺ʲπ_* (ℱ *)
```
This is a constructor for that spectral sequence as a lazy object from 
which subsequent pages and their entries can be computed.
```jldoctest
julia> R, (t,) = polynomial_ring(QQ, [:t]); # ring for the base with parameter t

julia> S, (x, y, z, w) = graded_polynomial_ring(R, [:x, :y, :z, :w]); # IP^3 over R

julia> f1 = sum(x^4 for x in gens(S); init=zero(S)); # a smooth K3 surface

julia> f0 = prod(x for x in gens(S); init=one(S)); # a degenerate quartic

julia> f = f0 + t*f1; # a family with degeneration for t=0

julia> S1 = graded_free_module(S, [0]);

julia> I, _ = ideal(f)*S1;

julia> M, _ = quo(S1, I); # a module representing the structure sheaf of the family

julia> res, _ = free_resolution(Oscar.SimpleFreeResolution, M);

julia> css = Oscar.CohomologySpectralSequence(S, res) # the spectral sequence
spectral sequence for the cohomology of graded complex with Betti table
degree: 0
---------
     0: 1
---------
 total: 1

julia> p1 = css[1]; # the first page

julia> [p1[i, j] for j in 0:-1:-3, i in 0:2] # indices have reversed signs here for internal reasons!
4×3 Matrix{FreeMod{QQMPolyRingElem}}:
 R^1  0    0
 0    0    0
 0    0    0
 0    R^1  0

julia> [!is_zero(css[2, i, j]) for j in 0:-1:-3, i in 0:2] # entries on the second page
4×3 Matrix{Bool}:
 1  0  0
 0  0  0
 0  0  0
 0  1  0

julia> map(p1, 1, -3) # the outgoing map on the first page at index (1, -3)
Module homomorphism
  from R^1
  to 0

```
"""
  function CohomologySpectralSequence(S::MPolyRing, comp::AbsHyperComplex)
    @assert isone(dim(comp)) "complex must be 1-dimensional"
    A = coefficient_ring(S)
    return new{typeof(S), typeof(A)}(S, A, comp, Dict{Int, CSSPage}())
  end
  
  function CohomologySpectralSequence(X::NormalToricVariety, comp::AbsHyperComplex; algorithm::Symbol=:ext)
    @assert isone(dim(comp)) "complex must be 1-dimensional"
    S = cox_ring(X)
    ctx = ToricCtx(X; algorithm)
    A = coefficient_ring(S)
    return new{typeof(S), typeof(A)}(S, A, comp, Dict{Int, CSSPage}(), ctx)
  end
  
  function CohomologySpectralSequence(X::NormalToricVariety, transfer::Map, comp::AbsHyperComplex)
    @assert isone(dim(comp)) "complex must be 1-dimensional"
    ctx = ToricCtxWithParams(X, transfer)
    S = graded_ring(ctx)
    R = coefficient_ring(S)
    return new{typeof(S), typeof(R)}(S, R, comp, Dict{Int, CSSPage}(), ctx)
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
  #@assert can_compute_index(graded_complex(spectral_sequence(cssp)), i) "index out of bounds"
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
  B, inc_B = can_compute_map(previous_page, i+p-1, j-p+2) ? image(map(previous_page, i+p-1, j-p+2)) : sub(H, elem_type(H)[])
  Z, inc_Z = can_compute_map(previous_page, i, j) ? kernel(map(previous_page, i, j)) : sub(H, gens(H))
  # The code below relies on a particular generating set for the modulus
  # so we need to carefully build that here.
  F = ambient_free_module(H)
  new_B, _ = sub(F, vcat(relations(H), repres.(gens(B))))
  g = ambient_representatives_generators(Z)
  return SubquoModule(ambient_free_module(H), g, ambient_representatives_generators(new_B))
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
  rngs = get_attribute(dom, :ranges)::Vector{UnitRange{Int}}
  cur_rng_ind = 1
  summands = get_attribute(dom, :direct_product)::Vector{typeof(dom)}
  for (v_num, v) in enumerate(gens(dom))
    buckets = Tuple{Int, FreeModElem}[]
    e = Int[]
    ind, _ = only(coordinates(v))
    while !(ind in rngs[cur_rng_ind])
      cur_rng_ind+=1
    end
    # manual way to project (saves allocations)
    vv = FreeModElem(coordinates(v)[rngs[cur_rng_ind]], summands[cur_rng_ind]) #canonical_projection(dom, k)(v)
    @assert !is_zero(vv)
    d = -dom_degs[cur_rng_ind]
    vv2 = cohomology_model_inclusion(ctx, d, j)(vv)
    @assert !is_zero(vv2)
    e = _minimal_exponent_vector(ctx, d)
    for (l, p) in coordinates(images_of_generators(orig_map)[cur_rng_ind])
      phi = multiplication_map(ctx, p, e, d, j)
      vvv = phi(vv2)
      @assert !is_zero(vvv)
      @assert !any(ll == l for (ll, _) in buckets)
      push!(buckets, (l, vvv))
    end
    @assert !is_empty(e)
    push!(result, (e, buckets))
  end
  return result
end

function produce_lifted_kernel_generators(cssp::CSSPage, i::Int, j::Int)
  #@assert can_compute_index(graded_complex(cssp), i) "index out of bounds"
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

  v = lifted_kernel_generators(prev_page, ii, jj) # lifts of the kernel generators on the previous page

  # The generators of the kernel Zₚ = ker ∂ₚ₋₁ are extracted from this page.
  u_next = [repres(simplify!(g)) for g in gens(cssp[i0, j0])] # effectively discard zero generators

  orig_cplx = graded_complex(cssp)
  # The degrees of the generators in the original graded complex 
  # determine the exponents of the numerators in the Cech complex 
  # which are required to properly represent the involved cohomology classes. 
  dom_degs = degrees_of_generators(orig_cplx[i0])
  inter_degs = degrees_of_generators(orig_cplx[ii])
  cod_degs = degrees_of_generators(orig_cplx[i])

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

  hor_maps = Dict{Vector{Int}, Map}()
  vert_maps = Dict{Vector{Int}, Map}()
  other_vert_maps = Dict{Vector{Int}, Map}()

  next_hor_maps = Dict{Vector{Int}, Map}()
  next_vert_maps = Dict{Vector{Int}, Map}()

  for (u_num, u) in enumerate(u_next)
    # We write `u` in the generators of `Zₚ₋₁` on which the previous 
    # map ∂ₚ₋₁ is defined.
    u_coords = p == 2 ? coordinates(u) : coordinates(u, prev_page[i0, j0])

    # To bring everything to the same denominator in the Cech complex 
    # we need to find a common exponent.
    exps = Vector{Int}[v[i][1] for (i, c) in u_coords]
    exps = vcat(exps, Vector{Int}[_minimal_exponent_vector(ctx, -d) for d in unique(inter_degs)])
    exps = vcat(exps, Vector{Int}[_minimal_exponent_vector(ctx, -d) for d in unique(cod_degs)])
    sup_exp = Int[maximum([e[i] for e in exps]; init=0) for i in 1:length(first(exps))]

    buckets = Tuple{Int, FreeModElem}[]
    # ∂ₚ₋₁(u) = 0, so the linear combination for `u` in the generators of `Zₚ₋₁` 
    # applied to the `lifted_kernel_generators` for ψₚ₋₁ gives us elements 
    # which can be lifted further
    for (l, c) in u_coords
      e, g = v[l]
      # bring everything to the same denominators
      for (j, v) in (e != sup_exp ? Tuple{Int, <:FreeModElem}[(j, ctx[e, sup_exp, -inter_degs[j]][jj](v)) for (j, v) in g] : g)
        ind = findfirst(a -> a[1]==j, buckets)
        if isnothing(ind)
          push!(buckets, (j, c*v))
        else
          _, w = buckets[ind]
          buckets[ind] = (j, w + c*v)
        end
      end
    end
      

    # Now the buckets hold lifts of the kernel generators on page `p`. 
    # We need to lift them once more through the double complex.

    # First we need to prepare the representatives
    coh_rem = zero(css[1, ii, jj])
    for (i, v) in buckets
      is_zero(v) && continue
      e0 = _minimal_exponent_vector(ctx, -inter_degs[i])
      pr1 = ctx[sup_exp, e0, -inter_degs[i]][jj]
      @assert parent(v) === domain(pr1)
      pr2 = cohomology_model_projection(ctx, -inter_degs[i], jj)
      @assert domain(pr2) === codomain(pr1)
      v_img = pr2(pr1(v))
      coh_rem += canonical_injection(parent(coh_rem), i)(v_img)
    end

    next_sup_exp = copy(sup_exp)

    # if coh_rem is not zero, then the element represented by `buckets` is not liftable 
    # through the Cech map. The adjustment can be done modulo the images of the incoming 
    # previous differentials.
    if !is_zero(coh_rem)
      # in this case we need to do adjustments by the images of the 
      # previous incoming maps.
      H = prev_page[ii, jj]
      I = H.quo
      @assert coh_rem in I
      c = coordinates(coh_rem, I)
      offset = 0
      ranges = UnitRange[]
      for q in 1:p-2
        kernel_lifts = lifted_kernel_generators(css[q], ii, jj)
        n = length(kernel_lifts)
        c_sub = c[offset+1:offset+n]
        push!(ranges, offset+1:offset+n)
        offset += n
      end

      # We first find a common denominator for the updated images
      next_sup_exp = copy(sup_exp)
      for (q, rng) in enumerate(ranges)
        c_sub = c[rng]
        is_zero(c_sub) && continue
        kernel_lifts = lifted_kernel_generators(css[q], ii, jj)
        exps = Vector{Int}[kernel_lifts[i][1] for (i, _) in c_sub]
        next_sup_exp = Int[maximum([e[i] for e in exps]; init=next_sup_exp[i]) for i in 1:ngens(S)]
      end

      # We build the corresponding linear combination and update the buckets
      for (q, rng) in enumerate(ranges)
        c_sub = c[rng]
        is_zero(c_sub) && continue
        kernel_lifts = lifted_kernel_generators(css[q], ii, jj)
        for (i, c) in c_sub
          (e, prev_buckets) = kernel_lifts[i]
          for (l, v) in prev_buckets
            ind = findfirst(a->a[1]==l, buckets)
            vv = -c*(e == next_sup_exp ? v : ctx[e, next_sup_exp, -inter_degs[l]][jj](v))
            if isnothing(ind)
              push!(buckets, (l, vv))
            else
              _, w = buckets[ind]
              if sup_exp != next_sup_exp
                w = ctx[sup_exp, next_sup_exp, -inter_degs[l]][jj](w)
              end
              buckets[ind] = (l, w + vv)
            end
          end
        end
      end
    end

    # Perform the lifting.
    lifted_buckets = Tuple{Int, FreeModElem}[]
    for (i, v) in buckets
      is_zero(v) && continue
      cech_map = map(ctx[sup_exp, -inter_degs[i]], jj+1)
      push!(lifted_buckets, (i, preimage(cech_map, v)))
    end

    sup_exp = next_sup_exp
    
    # The lifted elements still need to be mapped.
    mapped_buckets = Tuple{Int, FreeModElem}[]
    orig_map = map(graded_complex(cssp), ii)
    for (k, v) in lifted_buckets
      img_gen = images_of_generators(orig_map)[k]
      for (l, p) in coordinates(img_gen)
        phi = multiplication_map(ctx, p, sup_exp, -inter_degs[k], jj+1)
        vv = phi(v)
        is_zero(vv) && continue
        ind = findfirst(a->a[1] == l, mapped_buckets)
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
  return new_lifted_kernel_gens
end

# Some internal methods for debugging. 
# They are for turning the sparsely stored intermediate data for the 
# `lifted_kernel_generators` into actual mathematical objects in the 
# numerous direct limits and get instances of the maps on these objects 
# to check coherence and validity of the data.
function vertical_map(
    css::CohomologySpectralSequence, exps::Vector{Int}, i::Int, j::Int;
    domain::FreeMod=begin
      S = graded_ring(css)
      G = grading_group(S)
      cplx = graded_complex(css)
      ctx = pushforward_ctx(css)
      dom_degs = degrees_of_generators(cplx[i])
      direct_sum([ctx[exps, -d][j] for d in dom_degs])[1]
    end,
    codomain::FreeMod=begin
      S = graded_ring(css)
      G = grading_group(S)
      cplx = graded_complex(css)
      ctx = pushforward_ctx(css)
      cod_degs = degrees_of_generators(cplx[i])
      direct_sum([ctx[exps, -d][j-1] for d in cod_degs])[1]
    end
  )
  S = graded_ring(css)
  G = grading_group(S)
  cplx = graded_complex(css)
  ctx = pushforward_ctx(css)
  dom_degs = degrees_of_generators(cplx[i])
  img_gens = elem_type(codomain)[]
  cech_maps = [map(ctx[exps, -d], j) for d in dom_degs]
  result = hom(domain, codomain, elem_type(codomain)[zero(codomain) for _ in 1:ngens(domain)])
  @assert length(canonical_projections(domain)) == ngens(cplx[i])
  @assert length(canonical_injections(codomain)) == ngens(cplx[i])
  @assert length(cech_maps) == ngens(cplx[i])
  for (pr, inc, cech_map) in zip(canonical_projections(domain), canonical_injections(codomain), cech_maps)
    @assert Oscar.codomain(pr) === Oscar.domain(cech_map)
    @assert Oscar.codomain(cech_map) === Oscar.domain(inc)
    result += compose(pr, compose(cech_map, inc))
  end
  return result
end


function horizontal_map(
    css::CohomologySpectralSequence, exps::Vector{Int}, i::Int, j::Int;
    domain::FreeMod=begin
      S = graded_ring(css)
      G = grading_group(S)
      cplx = graded_complex(css)
      ctx = pushforward_ctx(css)
      dom_degs = degrees_of_generators(cplx[i])
      direct_sum([ctx[exps, -d][j] for d in dom_degs])[1]
    end,
    codomain::FreeMod=begin
      S = graded_ring(css)
      G = grading_group(S)
      cplx = graded_complex(css)
      ctx = pushforward_ctx(css)
      cod_degs = degrees_of_generators(cplx[i-1])
      direct_sum([ctx[exps, -d][j] for d in cod_degs])[1]
    end
  )
  S = graded_ring(css)
  G = grading_group(S)
  cplx = graded_complex(css)
  ctx = pushforward_ctx(css)
  dom_degs = degrees_of_generators(cplx[i])
  cod_degs = degrees_of_generators(cplx[i-1])
  orig_map = map(cplx, i)
  img_gens = elem_type(codomain)[]
  for (k, pr) in enumerate(canonical_projections(domain))
    is_zero(Oscar.codomain(pr)) && continue
    for v in gens(Oscar.codomain(pr))
      img = zero(codomain)
      for (l, p) in coordinates(images_of_generators(orig_map)[k])
        mm = multiplication_map(ctx, p, exps, -dom_degs[k], j)
        vv = mm(v)
        inc = canonical_injection(codomain, l)
        vvv = inc(vv)
        img += vvv
      end
      push!(img_gens, img)
    end
  end
  return hom(domain, codomain, img_gens)
end

function cohomology_injection_map(
    css::CohomologySpectralSequence, exps::Vector{Int}, i::Int, j::Int;
    codomain::FreeMod=begin
      S = graded_ring(css)
      G = grading_group(S)
      cplx = graded_complex(css)
      ctx = pushforward_ctx(css)
      cod_degs = degrees_of_generators(cplx[i])
      direct_sum([ctx[exps, -d][j] for d in cod_degs])[1]
    end
  )
  S = graded_ring(css)
  G = grading_group(S)
  cplx = graded_complex(css)
  ctx = pushforward_ctx(css)
  degs = degrees_of_generators(cplx[i])
  dom = css[1, i, j]
  result = hom(dom, codomain, elem_type(codomain)[zero(codomain) for _ in 1:ngens(dom)])
  for (k, pr) in enumerate(canonical_projections(dom))
    is_zero(Oscar.codomain(pr)) && continue
    inc = cohomology_model_inclusion(ctx, -degs[k], j)
    @assert Oscar.codomain(pr) === domain(inc)
    e0 = _minimal_exponent_vector(ctx, -degs[k])
    phi = ctx[e0, exps, -degs[k]][j]
    @assert Oscar.codomain(inc) === domain(phi)
    @assert Oscar.codomain(phi) === domain(canonical_injection(codomain, k))
    result += compose(compose(pr, inc), compose(phi, canonical_injection(codomain, k)))
  end
  return result
end

function check_sanity(cssp::CSSPage, i::Int, j::Int; error_on_false::Bool=true)
  p = page_number(cssp)
  css = spectral_sequence(cssp)
  ctx = pushforward_ctx(css)
  cplx = graded_complex(cssp)
  S = graded_ring(css)
  G = grading_group(S)

  # check plausibility of the `lifted_kernel_generators`
  vert_maps = Dict{Vector{Int}, FreeModuleHom}()
  secondary_vert_maps = Dict{Vector{Int}, FreeModuleHom}()
  hor_maps = Dict{Vector{Int}, FreeModuleHom}()
  secondary_hor_maps = Dict{Vector{Int}, FreeModuleHom}()
  if j - p + 1 >= -ngens(S) + rank(G) # if there is an incoming map
    for (e, buckets) in lifted_kernel_generators(cssp, i, j)
      if i > 0
        # sanity check of the input
        hor_map = get!(hor_maps, e) do
          horizontal_map(css, e, i, j)
        end
        v0 = sum(canonical_injection(domain(hor_map), l)(v) for (l, v) in buckets if !is_zero(v); init=zero(domain(hor_map)))
        res = is_zero(hor_map(v0))
        !res && (error_on_false ? error("horizontal map does not annihilate lifted kernel generator") : return false)
      end
      if j > -ngens(S) + rank(grading_group(S))
        vert_map = get!(vert_maps, e) do
          vertical_map(css, e, i, j)
        end
        v0 = sum(canonical_injection(domain(vert_map), l)(v) for (l, v) in buckets if !is_zero(v); init=zero(domain(vert_map)))
        res = is_zero(vert_map(v0))
        !res && (error_on_false ? error("vertical map does not annihilate lifted kernel generator") : return false)
      end
    end
  end
  # check sanity of the generators
  exps = [_minimal_exponent_vector(ctx, -d) for d in degrees_of_generators(cplx[i])]
  e = Int[maximum([a[i] for a in exps]) for i in 1:nvars(S)]
  for (k, g) in enumerate(gens(cssp[i, j]))
    rep_g = repres(g)

    if i > 0 && j < 0 && p > 1
      # sanity check of the input
      hor_map = get!(hor_maps, e) do
        horizontal_map(css, e, i, j)
      end
      inc = cohomology_injection_map(css, e, i, j; codomain=domain(hor_map))
      v0 = inc(rep_g)
      res = hor_map(v0)
      # `res` needs to be in the image of the vertical map
      vert_map = get!(secondary_hor_maps, e) do
        vertical_map(css, e, i-1, j+1; codomain=codomain(hor_map))
      end
      res in image(vert_map)[1] || (error_on_false ? error("horizontal map does not take the generator to the image of the incoming vertical map") : return false)
    end
    if j > -ngens(S) + rank(grading_group(S))
      vert_map = get!(vert_maps, e) do
        vertical_map(css, e, i, j)
      end
      inc = cohomology_injection_map(css, e, i, j; codomain=domain(vert_map))
      v0 = inc(rep_g)
      res = vert_map(v0)
      if can_compute_index(cplx, (i+1,))
        hor_map = get!(secondary_hor_maps, e) do
          horizontal_map(css, e, i+1, j-1; codomain=codomain(vert_map))
        end
        if codomain(vert_map) !== codomain(hor_map)
          hor_map = horizontal_map(css, e, i+1, j-1; codomain=codomain(vert_map))
        end
        res in image(hor_map)[1] || (error_on_false ? error("vertical map does not take the generator to the image of the horizontal map") : return false)
      else
        is_zero(res) && (error_on_false ? error("vertical map does not annihilate lifted kernel generator") : return false)
      end
    end
  end
  return true
end

function produce_map(cssp::CSSPage, i::Int, j::Int)
  #@assert can_compute_index(graded_complex(cssp), i) "index out of bounds"
  #@assert can_compute_index(graded_complex(cssp), i-page_number(cssp)) "index out of bounds"
  @assert j <= 0 "index out of bounds"
  @assert j + page_number(cssp) - 1 <= 0 "index out of bounds"
  #is_one(page_number(cssp)) && return produce_map_on_initial_page(cssp, i, j)
  #page_number(cssp) == 2 && return produce_map_on_second_page(cssp, i, j)

  p = page_number(cssp)
  S = graded_ring(cssp)
  css = spectral_sequence(cssp)
  ctx = pushforward_ctx(cssp)
  
  # the indices of the codomain
  i1 = i - p 
  j1 = j + p - 1

  # assemble the actual induced morphism
  dom = cssp[i, j]
  cod = cssp[i1, j1]

  # take the chance for early return
  if (is_known(is_zero, dom) && is_zero(dom)) || (is_known(is_zero, cod) && is_zero(cod))
    return hom(dom, cod, elem_type(cod)[zero(cod) for _ in 1:ngens(dom)])
  end
  
  new_lifted_kernel_gens = lifted_kernel_generators(cssp, i1, j1)
  F = ambient_free_module(cod)
  img_gens = elem_type(F)[]
  cod_degs = degrees_of_generators(graded_complex(cssp)[i1])
  vert_maps = Dict{Vector{Int}, FreeModuleHom}()
  other_vert_maps = Dict{Vector{Int}, FreeModuleHom}()
  hor_maps = Dict{Vector{Int}, FreeModuleHom}()
  incs = Dict{Vector{Int}, FreeModuleHom}()
  for (e, buckets) in new_lifted_kernel_gens # iterate through the representatives 
                                             # of the images
    img_gen = zero(F) # initialize the result for this generator
    for (k, v) in buckets # go through the buckets
      is_zero(v) && continue # discard zero elements 
      d1 = -cod_degs[k] # find the shift for this component
      e1 = _minimal_exponent_vector(ctx, d1) # get the minimal exponent for this 
                                             # component's cohomology complex
      @assert all(a <= b for (a, b) in zip(e1, e)) # make sure we have not lost information
                                                   # along the way by taking to small denoms
      pr = ctx[e, e1, d1][j1] # project down to the minimal complex
      vv = pr(v)
      is_zero(vv) && continue
      coh_pr = cohomology_model_projection(ctx, d1, j1) # project to the cohomology model
      img_gen += canonical_injection(F, k)(coh_pr(vv)) # add the component to the result
    end
    push!(img_gens, img_gen)
  end
  return hom(dom, cod, elem_type(cod)[cod(v) for v in img_gens]; check=false)
end

relations(F::FreeMod) = elem_type(F)[]

function multiplication_map(
    ctx::Union{PushForwardCtx, ToricCtx},
    p::MPolyDecRingElem,
    e0::Vector{Int}, d0::FinGenAbGroupElem, 
    j::Int
  )
  cache = get!(ctx.mult_map_cache, (e0, d0, j)) do
    WeakKeyDict{typeof(p), Map}()
  end
  neg_res = get(cache, -p, nothing)
  !isnothing(neg_res) && return MapFromFunc(domain(neg_res), codomain(neg_res), v->-neg_res(v))
  return get!(cache, p) do
    d1 = d0 + degree(p; check=false)
    dom_cplx = ctx[e0, d0]
    cod_cplx = ctx[e0, d1]
    dom = dom_cplx[j]
    cod = cod_cplx[j]
    dom_strand_inc = inclusion_map(dom_cplx)[j]
    cod_strand_pr = projection_map(cod_cplx)[j]
    return MapFromFunc(dom, cod, v->cod_strand_pr(p*dom_strand_inc(v)))
  end
end

function multiplication_map(
    ctx::ToricCtxWithParams,
    p::MPolyDecRingElem,
    e0::Vector{Int}, d0::FinGenAbGroupElem, 
    j::Int
  )
  # we store the multiplication maps in weak dictionaries. 
  # The same polynomial tends to appear multiple times in the matrices of a 
  # given complex; for instance, the Koszul complex has the same `n` polynomials 
  # all over itself, but with different signs. Using weak dictionaries should 
  # prevent us from memory overflows. 
  cache = get!(ctx.mult_map_cache, (e0, d0, j)) do
    WeakKeyDict{typeof(p), Map}()
  end
  # If the map is already known for the negative of a given polynomial, use that.
  neg_res = get(cache, -p, nothing)
  !isnothing(neg_res) && return MapFromFunc(domain(neg_res), codomain(neg_res), v->-neg_res(v))
  # We want to use as little cache as possible while still providing reasonable speed. 
  # The bases of the free modules in question are usually quite big. Storing a list of 
  # images of generators is therefore rather costly, even though they will be sparse. 
  # Given the nature of the setup, the relevant information to compute the image of an 
  # element is already stored in a format which allows for fast evaluation by means 
  # different from taking a linear combination. 
  return get!(cache, p) do
    d1 = d0 + degree(p; check=false)
    dom_cplx = ctx[e0, d0]
    cod_cplx = ctx[e0, d1]
    dom = dom_cplx[j]
    cod = cod_cplx[j]

    S = cox_ring(toric_variety(ctx))
    R = base_ring(cod)
    mult_cache = Dict{Int, sparse_row_type(R)}() # store the image of the `i`-th generator
    return MapFromFunc(dom, cod, function(v)
                         res_coords = sparse_row(R)
                         for (i, q) in coordinates(v)
                           # `q` will be an element of `R`, the parameter ring.
                           img_gen = get!(mult_cache, i) do
                             new_row = sparse_row(R)
                             for (c, e) in zip(AbstractAlgebra.coefficients(p), 
                                               AbstractAlgebra.exponent_vectors(p))
                               # While the coefficient `c` is in `R`, the multiplication
                               # map ϕₑfor the monomial for `e` is determined in the ring 
                               # `S` over `kk`. We need to construct ϕₑ ⊗ R.
                               mon = S([one(QQ)], [e])
                               mon_mult = multiplication_map(ctx.pure_ctx, mon, e0, d0, j)
                               # `domain(mon_mult)` and `dom` have generating sets so that 
                               # `dom[i] = domain(mon_mult)[i] ⊗ R`. 
                               mon_gen = domain(mon_mult)[i]
                               row = map_entries(R, coordinates(mon_mult(mon_gen)))
                               new_row = Hecke.add_scaled_row!(row, new_row, c)
                             end
                             new_row
                           end
                           res_coords = Hecke.add_scaled_row!(img_gen, res_coords, q)
                         end
                         return cod(res_coords)
                       end
                      )
  end
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

function set_global_exponent_vector!(css::CohomologySpectralSequence, e::Vector{Int})
  @assert css.pfctx isa ToricCtx "manual setting of exponent vectors is only supported in the toric setup"
  return set_global_exponent_vector!(css.pfctx, e)
end

function set_global_exponent_vector!(css::CohomologySpectralSequence, k::Int)
  @assert css.pfctx isa ToricCtx "manual setting of exponent vectors is only supported in the toric setup"
  return set_global_exponent_vector!(css.pfctx, Int[k for _ in 1:length(cech_complex_generators(css.pfctx))])
end

