function _maximum(v::Vector{FinGenAbGroupElem})
  @assert !is_empty(v)
  G = parent(first(v))
  m = Vector{Int}()
  for i in 1:ngens(G)
    push!(m, maximum([w[i] for w in v]))
  end
  return sum(a*G[i] for (i, a) in enumerate(m); init=zero(G))
end

function _regularity_bound(F::FreeMod)
  @assert is_graded(F)
  is_zero(ngens(F)) && return zero(grading_group(F))
  degs = degrees_of_generators(F)
  result = _maximum(degs)
end

function _regularity_bound(comp::AbsHyperComplex, rng::UnitRange)
  @assert isone(dim(comp))
  m = [_regularity_bound(comp[i]) for i in rng]
  return _maximum(m)
end

function _derived_pushforward(M::FreeMod)
  S = base_ring(M)
  G = grading_group(M)
  r = rank(G)
  variables = [[x for x in gens(S) if degree(x) == G[i]] for i in 1:r]
  dims = [length(x)-1 for x in variables]

  d = _regularity_bound(M) # the degrees of the generators
  d = d - sum(n*G[i] for (i, n) in enumerate(dims); init=zero(G))
  d = sum((d[i] < 0 ? 0 : d[i])*G[i] for i in 1:r; init=zero(G))

  g = vcat([[x^(Int(d[i])) for x in v] for (i, v) in enumerate(variables)]...)
  kosz = [shift(Oscar.HomogKoszulComplex(S, [x^(Int(d[i])) for x in v])[1:length(v)], 1) for (i, v) in enumerate(variables)]
  K = simplify(total_complex(tensor_product(kosz)))

  KoM = hom(K, M)
  st = strand(KoM, zero(G))[1]
  return st
end

@doc raw"""
    _derived_pushforward(M::SubquoModule)

We consider a graded module `M` over a standard graded polynomial ring 
``S = A[x₀,…,xₙ]`` as a representative of a coherent sheaf ``ℱ`` on 
relative projective space ``ℙ ⁿ_A``. Then we compute ``Rπ_* ℱ`` as a
complex of ``A``-modules where ``π : ℙ ⁿ_A → Spec(A)`` is the projection 
to the base. 
"""
function _derived_pushforward(M::SubquoModule)
  S = base_ring(M)
  n = ngens(S)-1

  d = _regularity_bound(M) - n
  d = (d < 0 ? 0 : d)

  return _derived_pushforward(M, d)
end

# Method for the multigraded case
function _derived_pushforward(M::FreeMod{T}, seq::Vector{T}) where {T <: RingElem}
  return _derived_pushforward(ZeroDimensionalComplex(M), seq)
end

function _derived_pushforward(M::SubquoModule{T}, seq::Vector{T}) where {T <: RingElem}
  res, _ = free_resolution(Oscar.SimpleFreeResolution, M)
  return _derived_pushforward(res, seq)
  S = base_ring(M)
  @assert all(parent(x) === S for x in seq)
  n = length(seq)
  Sd = graded_free_module(S, [0 for i in 1:n])
  v = sum(seq[i]*Sd[i] for i in 1:n; init=zero(Sd))
  kosz = koszul_complex(Oscar.KoszulComplex, v)
  K = shift(Oscar.DegreeZeroComplex(kosz)[1:n+1], 1)

  res, _ = free_resolution(Oscar.SimpleFreeResolution, M)
  KoM = hom(K, res)
  tot = total_complex(KoM)
  tot_simp = simplify(tot)

  G = grading_group(M)
  st = strand(tot_simp, zero(G))
  return st[1]
end

function _derived_pushforward(comp::AbsHyperComplex, seq::Vector{T}) where {T <: RingElem}
  @assert dim(comp) <= 1
  @assert !isempty(seq)
  S = parent(first(seq))
  @assert all(parent(x) === S for x in seq)
  n = length(seq)
  Sd = graded_free_module(S, [0 for i in 1:n])
  v = sum(seq[i]*Sd[i] for i in 1:n; init=zero(Sd))
  kosz = koszul_complex(Oscar.KoszulComplex, v)
  K = shift(Oscar.DegreeZeroComplex(kosz)[1:n+1], 1)

  KoM = hom(K, comp)
  tot = total_complex(KoM)
  #tot_simp = simplify(tot)

  G = grading_group(S)
  st, _ = strand(tot, zero(G))
  return simplify(st)
end


# Method for the standard graded case
function _derived_pushforward(M::SubquoModule, bound::Int)
  S = base_ring(M)
  n = ngens(S)-1

  Sd = graded_free_module(S, [0 for i in 1:ngens(S)])
  v = sum(x^bound*Sd[i] for (i, x) in enumerate(gens(S)); init=zero(Sd))
  kosz = koszul_complex(Oscar.KoszulComplex, v)
  K = shift(Oscar.DegreeZeroComplex(kosz)[1:n+1], 1)

  res, _ = free_resolution(Oscar.SimpleFreeResolution, M)
  KoM = hom(K, res)
  tot = total_complex(KoM)
  tot_simp = simplify(tot)

  st = strand(tot_simp, 0)
  return st[1]
end

# The code below is supposedly faster, because it does not need to pass 
# through a double complex. 
#
# In practice, however, it is impossible to program the iterators for
# monomials in quotient modules in a sufficiently efficient way for this 
# to work. Hence, it's disabled for the time being.
#=
function _derived_pushforward(M::SubquoModule{T}) where {T <:MPolyRingElem{<:FieldElem}}
  S = base_ring(M)
  n = ngens(S)-1

  d = _regularity_bound(M) - n
  d = (d < 0 ? 0 : d)

  Sd = graded_free_module(S, [0 for i in 1:ngens(S)])
  v = sum(x^d*Sd[i] for (i, x) in enumerate(gens(S)); init=zero(Sd))
  kosz = koszul_complex(Oscar.KoszulComplex, v)
  K = shift(Oscar.DegreeZeroComplex(kosz)[1:n+1], 1)

  pres = presentation(M)
  N = cokernel(map(pres, 1))

  KoM = hom(K, N)
  st = strand(KoM, 0)
  return st[1]
end
=#

function rank(phi::FreeModuleHom{FreeMod{T}, FreeMod{T}, Nothing}) where {T<:FieldElem}
  return ngens(domain(phi)) - nrows(kernel(sparse_matrix(phi), side = :left))
end


@doc raw"""
    simplify(c::ComplexOfMorphisms{ChainType}) where {ChainType<:ModuleFP}

For a complex `c` of free modules over some `base_ring` `R` this looks for 
unit entries `u` in the representing matrices of the (co-)boundary maps and 
eliminates the corresponding generators in domain and codomain of that map. 

It returns a triple `(d, f, g)` where 
  * `d` is the simplified complex, 
  * `f` is a `Vector` of morphisms `f[i] : d[i] → c[i]` and 
  * `g` is a `Vector` of morphisms `g[i] : c[i] → d[i]` 
such that the composition `f ∘ g` is homotopy equivalent to the identity 
and `g ∘ f` is the identity on `d`.
"""
function simplify(c::ComplexOfMorphisms{ChainType}) where {ChainType<:ModuleFP}
  # the maps from the new to the old complex
  phi = morphism_type(ChainType)[identity_map(c[first(range(c))])] 
  # the maps from the old to the new complex
  psi = morphism_type(ChainType)[identity_map(c[first(range(c))])]
  # the boundary maps in the new complex
  new_maps = morphism_type(ChainType)[]
  for i in map_range(c)
    f = compose(last(phi), map(c, i))
    M = domain(f)
    N = codomain(f)
    A = sparse_matrix(f)
    Acopy = copy(A)
    S, Sinv, T, Tinv, ind = _simplify_matrix!(A)
    m = nrows(A)
    n = ncols(A)
    I = [i for (i, _) in ind]
    I = [i for i in 1:m if !(i in I)]
    J = [j for (_, j) in ind]
    J = [j for j in 1:n if !(j in J)]
    img_gens_dom = elem_type(M)[sum(c*M[j] for (j, c) in S[i]; init=zero(M)) for i in I]
    img_gens_cod = elem_type(N)[sum(c*N[i] for (i, c) in T[j]; init=zero(N)) for j in J]
    new_cod = _make_free_module(N, img_gens_cod)
    new_dom = _make_free_module(M, img_gens_dom)
    cod_map = hom(new_cod, N, img_gens_cod)
    dom_map = hom(new_dom, M, img_gens_dom)

    img_gens_dom = elem_type(new_dom)[]
    for i in 1:m
      w = Sinv[i]
      v = zero(new_dom)
      for j in 1:length(I)
        a = w[I[j]]
        !iszero(a) && (v += a*new_dom[j])
      end
      push!(img_gens_dom, v)
    end

    img_gens_cod = elem_type(new_cod)[]
    for i in 1:n
      w = Tinv[i]
      new_entries = Vector{Tuple{Int, elem_type(base_ring(w))}}()
      for (real_j, b) in w
        j = findfirst(==(real_j), J)
        j === nothing && continue
        push!(new_entries, (j, b))
      end
      w_new = sparse_row(base_ring(w), new_entries)
      push!(img_gens_cod, FreeModElem(w_new, new_cod))
    end

    cod_map_inv = hom(N, new_cod, img_gens_cod)
    dom_map_inv = hom(M, new_dom, img_gens_dom)

    v = gens(new_cod)
    img_gens = elem_type(new_cod)[]
    for k in 1:length(I)
      w = A[I[k]]
      push!(img_gens, sum(w[J[l]]*v[l] for l in 1:length(J); init=zero(new_cod)))
    end
    g = hom(new_dom, new_cod, img_gens)
    if !isempty(new_maps)
      new_maps[end] = compose(last(new_maps), dom_map_inv)
      @assert codomain(last(new_maps)) === domain(g)
    end
    push!(new_maps, g)
    phi[end] = compose(dom_map, phi[end])
    push!(phi, cod_map)
    psi[end] = compose(psi[end], dom_map_inv)
    push!(psi, cod_map_inv)
    # Enable for debugging
    @assert compose(g, cod_map) == compose(dom_map, f)
  end
  result = ComplexOfMorphisms(ChainType, new_maps, seed=last(range(c)), typ=typ(c), check=true)
  return result, phi, psi
end

function _vdim(M::SubquoModule)
  F = ambient_free_module(M)
  all(repres(x) in gens(F) for x in gens(M)) || return _vdim(presentation(M)[-1])

  return Singular.vdim(Singular.std(singular_generators(M.quo.gens)))
end

########################################################################
# A context object for computing the spectral sequence associated 
# to a the ̌Cech double complex for a complex of coherent sheaves. 
########################################################################
mutable struct PushForwardCtx
  S::MPolyRing
  variable_groups::Vector{Vector{<:MPolyRingElem}}
  var_group_indices::Vector{Vector{Int}}
  dims::Vector{Int}
  truncated_cech_complexes::Dict{Vector{Int}, AbsHyperComplex}
  inclusions::Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}
  projections::Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}
  strands::Dict{Vector{Int}, Dict}
  strand_inclusions::Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}
  strand_projections::Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}
  cohomology_models::Dict{FinGenAbGroupElem, AbsHyperComplex}
  cohomology_inclusions::Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}
  cohomology_projections::Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}
  # mult_map_cache::Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, Dict}
  mult_map_cache::Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, WeakKeyDict}
  S1::AbsHyperComplex

  function PushForwardCtx(S::MPolyRing)
    G = grading_group(S)
    var_grp_ind = Vector{Vector{Int}}()
    for i in 1:rank(G)
      push!(var_grp_ind, [j for j in 1:ngens(S) if degree(S[j]) == G[i]])
    end
    variable_groups = [gens(S)[ind] for ind in var_grp_ind]
    return new(S, 
               variable_groups,
               var_grp_ind,
               [length(v)-1 for v in variable_groups],
               Dict{Vector{Int}, AbsHyperComplex}(),
               Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}(),
               Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}(),
               Dict{Vector{Int}, Dict}(),
               Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}(),
               Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}(), 
               Dict{FinGenAbGroupElem, AbsHyperComplex}(),
               Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}(),
               Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}(),
               # Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, Dict}()
               Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, WeakKeyDict}()
              )
  end
end

graded_ring(ctx::PushForwardCtx) = ctx.S
variable_groups(ctx::PushForwardCtx) = ctx.variable_groups
variable_group_indices(ctx::PushForwardCtx, i::Int) = ctx.var_group_indices[i]
variable_group(ctx::PushForwardCtx, i::Int) = ctx.variable_groups[i]
number_of_factors(ctx::PushForwardCtx) = length(ctx.variable_groups)
dimensions(ctx::PushForwardCtx) = ctx.dims
dimension(ctx::PushForwardCtx, i::Int) = ctx.dims[i]

function ring_as_hypercomplex(ctx::PushForwardCtx)
  if !isdefined(ctx, :S1)
    S = graded_ring(ctx)
    ctx.S1 = ZeroDimensionalComplex(graded_free_module(S, [zero(grading_group(S))]))
  end
  return ctx.S1
end

function getindex(ctx::PushForwardCtx, alpha::Vector{Int})
  @assert all(>=(0), alpha)
  return get!(ctx.truncated_cech_complexes, alpha) do
    S = graded_ring(ctx)
    G = grading_group(S)
    cod = ring_as_hypercomplex(ctx)
    kosz = [shift(hom(Oscar.HomogKoszulComplex(S, elem_type(S)[S[i]^alpha[i] for i in variable_group_indices(ctx, j)]), cod)[-dimension(ctx, j)-1:-1], -1) for j in 1:number_of_factors(ctx)]
    K = total_complex(tensor_product(kosz))
    return K
  end
end

function getindex(ctx::PushForwardCtx, alpha::Vector{Int}, d::FinGenAbGroupElem)
  G = parent(d)
  S = graded_ring(ctx)
  @assert G === grading_group(S)
  strands = get!(ctx.strands, alpha) do 
    Dict{typeof(d), AbsHyperComplex}()
  end
  return get!(strands, d) do
    #offset = sum(a*degree(x) for (x, a) in zip(gens(S), alpha); init=zero(G))
    strand(ctx[alpha], d)[1]
  end
end

function cohomology_model(ctx::PushForwardCtx, d::FinGenAbGroupElem)
  get!(ctx.cohomology_models, d) do
    simplify(ctx[_minimal_exponent_vector(ctx, d), d])
  end
end

function cohomology_model_inclusion(ctx::PushForwardCtx, d::FinGenAbGroupElem, i::Int)
  h = cohomology_model(ctx, d)
  to_orig = map_to_original_complex(h)[i]
  @assert domain(to_orig) === h[i]
  return to_orig
end

function cohomology_model_projection(ctx::PushForwardCtx, d::FinGenAbGroupElem, i::Int)
  h = cohomology_model(ctx, d)
  from_orig = map_from_original_complex(h)
  @assert codomain(from_orig) === h
  return from_orig[i]
end

# return the minimal exponent vector `alpha` such that the whole 
# cohomology in degree `d` is contained in the truncated ̌Cech-complex for `alpha`
function _minimal_exponent_vector(ctx::PushForwardCtx, d::FinGenAbGroupElem)
  S = graded_ring(ctx)
  G = grading_group(S)
  @assert parent(d) === G
  result = [0 for _ in 1:ngens(S)]
  for i in 1:number_of_factors(ctx)
    inds = variable_group_indices(ctx, i)
    di = Int(d[i])
    di >= 0 && continue # Nothing to do in this case
    for j in inds
      result[j] = -di < dimension(ctx, i) ? 0 : -di - dimension(ctx, i)
    end
  end
  @assert length(result) == length(cech_complex_generators(ctx))
  return result
end

function getindex(ctx::PushForwardCtx, alpha::Vector{Int}, beta::Vector{Int})
  @assert all(a <= b for (a, b) in zip(alpha, beta))
  return get!(ctx.inclusions, (alpha, beta)) do
    S = graded_ring(ctx)
    c_alpha = ctx[alpha]::TotalComplex
    c_beta = ctx[beta]::TotalComplex
    c_alpha_orig = original_complex(c_alpha)::HCTensorProductComplex
    c_beta_orig = original_complex(c_beta)::HCTensorProductComplex
    facs = AbsHyperComplexMorphism[]
    for (a_fac, b_fac, d) in zip(factors(c_alpha_orig), factors(c_beta_orig), dimensions(ctx))
      # both a_fac and b_fac are shifted, truncated hom-complexes of Koszul complexes
      a_fac_unshift = original_complex(a_fac::ShiftedHyperComplex)
      b_fac_unshift = original_complex(b_fac::ShiftedHyperComplex)
      a_fac_untrunc = original_complex(a_fac_unshift::HyperComplexView)
      b_fac_untrunc = original_complex(b_fac_unshift::HyperComplexView)
      a_kosz = domain(a_fac_untrunc::HomComplex)
      b_kosz = domain(b_fac_untrunc::HomComplex)
      a_seq = sequence(a_kosz)
      b_seq = sequence(b_kosz)
      trans_mat = sparse_matrix(S, 0, d+1)
      for (i, p, q) in zip(1:d+1, a_seq, b_seq)
        push!(trans_mat, sparse_row(S, [(i, divexact(q, p))]))
      end
      ind_kosz = Oscar.InducedKoszulMorphism(b_kosz, a_kosz; transition_matrix=trans_mat)
      ind_hom = hom(ind_kosz, codomain(a_fac_untrunc); domain=a_fac_untrunc, codomain=b_fac_untrunc)
      ind_trunc = ind_hom[-d-1:-1, domain=a_fac_unshift, codomain=b_fac_unshift]
      ind_shift = shift(ind_trunc, [-1]; domain=a_fac, codomain=b_fac)
      push!(facs, ind_shift)
    end
    tensor_map = tensor_product(facs; domain=c_alpha_orig, codomain=c_beta_orig)
    tot_map = total_complex(tensor_map; domain=c_alpha, codomain=c_beta)
    tot_map
  end
end

function getindex(ctx::PushForwardCtx, alpha::Vector{Int}, beta::Vector{Int}, d::FinGenAbGroupElem)
  if all(a <= b for (a, b) in zip(alpha, beta))
    return get!(ctx.strand_inclusions, (alpha, beta, d)) do 
      strand(ctx[alpha, beta], d; domain=ctx[alpha, d], codomain=ctx[beta, d])
    end
  elseif all(a >= b for (a, b) in zip(alpha, beta))
    return get!(ctx.strand_projections, (alpha, beta, d)) do
      SummandProjection(ctx[beta, alpha, d])
    end
  end
  error("neither sector is fully contained in the other")
end

########################################################################
# A context object for computing the spectral sequence associated 
# to a the ̌Cech double complex on a toric variety
########################################################################
mutable struct ToricCtx
  X::NormalToricVariety
  S::MPolyRing
  algorithm::Symbol
  truncated_cech_complexes::Dict{Vector{Int}, AbsHyperComplex}
  inclusions::Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}
  projections::Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}
  strands::Dict{Vector{Int}, Dict}
  simplified_strands::IdDict{AbsHyperComplex, AbsHyperComplex}
  strand_inclusions::Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}
  strand_projections::Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}
  induced_cohomology_maps::Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem, Int}, Map}
  cohomology_models::Dict{FinGenAbGroupElem, AbsHyperComplex}
  cohomology_inclusions::Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}
  cohomology_projections::Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}
  # mult_map_cache::Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, Dict}
  mult_map_cache::Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, WeakKeyDict}
  proportionality_factors::Union{Tuple{Int, Int}, Nothing}
  exp_vec_cache::Dict{FinGenAbGroupElem, Vector{Int}} # Caching the _minimal_exponent_vector s
  S1::AbsHyperComplex
  cech_gens::Vector{<:MPolyRingElem}
  fixed_exponent_vector::Vector{Int} # If this field is set, then it is used for all computations
  D::Dict # An auxiliary cache; see `_minimal_exponent_vector` for more info.

  function ToricCtx(X::NormalToricVariety; algorithm::Symbol=:ext)
    S = cox_ring(X)
    G = grading_group(S)
    return new(X, S, algorithm,
               Dict{Vector{Int}, AbsHyperComplex}(),
               Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}(),
               Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}(),
               Dict{Vector{Int}, Dict}(),
               IdDict{AbsHyperComplex, AbsHyperComplex}(),
               Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}(),
               Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}(), 
               Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem, Int}, Map}(), 
               Dict{FinGenAbGroupElem, AbsHyperComplex}(),
               Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}(),
               Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}(),
               # Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, Dict}()
               Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, WeakKeyDict}(), 
               nothing,
               Dict{FinGenAbGroupElem, Vector{Int}}()
              )
  end
end

toric_variety(ctx::ToricCtx) = ctx.X
graded_ring(ctx::ToricCtx) = ctx.S


function ring_as_hypercomplex(ctx::ToricCtx)
  if !isdefined(ctx, :S1)
    S = graded_ring(ctx)
    ctx.S1 = ZeroDimensionalComplex(graded_free_module(S, [zero(grading_group(S))]))
  end
  return ctx.S1
end

function cech_complex_generators(ctx::ToricCtx)
  if !isdefined(ctx, :cech_gens)
    ctx.cech_gens = gens(irrelevant_ideal(toric_variety(ctx)))
  end
  return ctx.cech_gens::Vector{elem_type(graded_ring(ctx))}
end

function getindex(ctx::ToricCtx, alpha::Vector{Int})
  @assert all(>=(0), alpha)
  return get!(ctx.truncated_cech_complexes, alpha) do
    S = graded_ring(ctx)
    G = grading_group(S)
    cod = ring_as_hypercomplex(ctx)
    cech_gens = cech_complex_generators(ctx)
    n = length(cech_gens)
    if ctx.algorithm == :ext
      K = hom(free_resolution(SimpleFreeResolution, 
                              ideal(S, elem_type(S)[x^i for (x, i) in zip(cech_gens, alpha)]))[1],
              cod)
      return K
    elseif ctx.algorithm == :cech
      kosz = hom(shift(HomogKoszulComplex(S, elem_type(S)[x^i for (x, i) in zip(cech_gens, alpha)])[1:n], 1), cod)
      return kosz
    else
      error("algorithm not recognized")
    end
  end
end

function getindex(ctx::ToricCtx, alpha::Vector{Int}, d::FinGenAbGroupElem)
  G = parent(d)
  S = graded_ring(ctx)
  @assert G === grading_group(S)
  strands = get!(ctx.strands, alpha) do 
    Dict{typeof(d), AbsHyperComplex}()
  end
  return get!(strands, d) do
    strand(ctx[alpha], d)[1]
  end
end

function simplified_strand(ctx::ToricCtx, alpha::Vector{Int}, d::FinGenAbGroupElem)
  str = ctx[alpha, d]
  return get!(ctx.simplified_strands, str) do
    simplify(str; with_homotopy_maps=true)
  end
end

function induced_cohomology_map(
    ctx::ToricCtx, e0::Vector{Int},
    e1::Vector{Int}, d::FinGenAbGroupElem,
    i::Int
  )
  return get!(ctx.induced_cohomology_maps, (e0, e1, d, i)) do 
    if all(a <= b for (a, b) in zip(e0, e1))
      return compose(map_to_original_complex(simplified_strand(ctx, e0, d))[i],
                     compose(
                             ctx[e0, e1, d][i],
                             map_from_original_complex(simplified_strand(ctx, e1, d))[i]
                            )
                    )
    elseif all(a >= b for (a, b) in zip(e0, e1))
      return inv(induced_cohomology_map(ctx, e1, e0, d, i))
    else
      error("not implemented")
    end
  end
end
function simplified_strand_homotopy(
    ctx::ToricCtx, alpha::Vector{Int}, d::FinGenAbGroupElem, p::Int
  )
  return homotopy_map(simplified_strand(ctx, alpha, d), p)
end

function simplified_strand_inclusion(
    ctx::ToricCtx, alpha::Vector{Int}, d::FinGenAbGroupElem, p::Int
  )
  return map_to_original_complex(simplified_strand(ctx, alpha, d))[p]
end

function simplified_strand_projection(
    ctx::ToricCtx, alpha::Vector{Int}, d::FinGenAbGroupElem, p::Int
  )
  return map_from_original_complex(simplified_strand(ctx, alpha, d))[p]
end


function cohomology_model(ctx::ToricCtx, d::FinGenAbGroupElem)
  get!(ctx.cohomology_models, d) do
    simplified_strand(ctx, _minimal_exponent_vector(ctx, d), d)
  end
end

function cohomology_model_inclusion(ctx::ToricCtx, d::FinGenAbGroupElem, i::Int)
  h = cohomology_model(ctx, d)
  c = ctx[_minimal_exponent_vector(ctx, d), d]
  to_orig = map_to_original_complex(h)[i]
  @assert domain(to_orig) === h[i]
  return to_orig
end

function cohomology_model_projection(ctx::ToricCtx, d::FinGenAbGroupElem, i::Int)
  h = cohomology_model(ctx, d)
  c = ctx[_minimal_exponent_vector(ctx, d), d]
  from_orig = map_from_original_complex(h)
  @assert codomain(from_orig) === h
  return from_orig[i]
end

function get_lattice_points(A::Any, Q::BitVector, m::Vector{Int}, n::Int)
  minus_indices = findall(x -> x != 0, Q)
  W = copy(A); W[:, minus_indices] = -W[:,minus_indices]
  P = polyhedron((-Matrix{Int}(identity_matrix(ZZ,n)),zeros(Int,n)),(W,m + A*Q))#First tuples is an affine halfspace that specifies that all exponents are non-negative
  #Second tuple is a hyperplane that specifies that the degree of the corresponding rationome is m
  points = Vector{Vector{Int}}()
  list = try 
    collect(lattice_points(P))#It is sometimes possible that there will be inf solutions which prob means that the corresponding rationome will contribute 0(secondary complex should be 0)
  catch
    Vector{Vector{Int}}()
  end
  for point in list
    point = Vector{Int}(point)
    point[minus_indices] = -(point[minus_indices] .+ 1)
    push!(points, point)
  end
  return points
end

function traverse_SR_ideal(v::NormalToricVariety)
#This traverses the whole powerset of the Stanley-Reisner ideal
  D=Dict{BitVector, Vector{Int}}()
  SR = [BitVector(collect(exponents(p))[1]) for p in gens(stanley_reisner_ideal(v))]
  n_coords = length(SR[1])
  n_sr = length(SR)
  for k in 0:length(SR)
    P_k = subsets(SR,k)#Subsets of size k
    for S in P_k
      !is_empty(S) ? Q = reduce(.|, S) : Q = falses(n_rays(v))
      if !haskey(D,Q)
        D[Q] = zeros(Int, n_sr + n_coords)
      end
      N = sum(Q) - k #This element will contribute to the N-th cohomology group

      D[Q][N+1+n_sr] += 1
    end
  end
  return D
end

function cohomology_support(v::NormalToricVariety, m::Vector{Int}; D = Dict())#returns vectors of exponents of rationoms that generate the cohomologies of the line bundle O(m) of v
  #m should be in the class group of v

  Classes=[divisor_class(toric_divisor_class(x)) for x in torusinvariant_prime_divisors(v)]
  Classes=[Int.(getindex(x,collect(1:length(m)))) for x in Classes]#Using this one needs to be sure that m is valid

  A = hcat(Classes...)
  m_lift = solve(matrix(ZZ,A), ZZ.(m), side=:right)
  m_dual = Vector{Int}(matrix(ZZ,A)*(-m_lift .- 1))

  if isempty(D)
    D = traverse_SR_ideal(v)
  end
  list = Dict{BitVector, Vector{Vector{Int}}}()

  for (Q,c_degrees) in D
    list[Q] = get_lattice_points(A, Q, m, n_rays(v))
  end

  for (Q,c_degrees) in D
    if !isempty(list[Q])

      Q_rest = .~Q
      if length(findall(x -> x!=0,c_degrees)) > 1
        if !haskey(D, Q_rest)
          list[Q] = []#The dual subset wasn't in the list so we omit any lattice points
          continue
        end

        #This could omit a small amount of lattice points, which could improve the optimal_k bound very slightly, but could be outcommented for some performance
        ########
        list_rest = get_lattice_points(A, Q_rest, m_dual, n_rays(v))
        if length(list_rest)==0
          list[Q] = []#The dual subset was in the list, but had 0 or inf dual lattice points, so for Serre duality to make sense, this subset of SR cannot contribute
          continue
        end
        ########
      end
    end
  end

  filter!(x -> !isempty(x[2]), list)
  return collect(Iterators.flatten(values(list))), D
end

function set_global_exponent_vector!(ctx::ToricCtx, v::Vector{Int})
  @assert length(v) == length(cech_complex_generators(ctx)) "exponent vector needs to have the same length as the generators for the `irrelevant_ideal`"
  ctx.fixed_exponent_vector = v
end

function set_global_exponent_vector!(ctx::ToricCtx, k::Int)
  return set_global_exponent_vector!(ctx, Int[k for _ in 1:length(cech_complex_generators(ctx))])
end

# return the minimal exponent vector `alpha` such that the whole 
# cohomology in degree `m` is contained in the truncated ̌Cech-complex for `alpha`
function _minimal_exponent_vector(ctx::ToricCtx, m::FinGenAbGroupElem)
  if isdefined(ctx, :fixed_exponent_vector)
    return ctx.fixed_exponent_vector
  end
  if ctx.algorithm == :cech
    error("dynamic computation of the minimal exponent vector is not implemented for cech cohomology")
  end
  # The following is based on the cohomCalg algorithm(See [BJRR10, BJRR10*1](@cite)), 
  # but only part of the algorithm is executed and explicit lattice points are 
  # calculated for each polyhedron.
  return get!(ctx.exp_vec_cache, m) do
    rationoms, ctx.D = cohomology_support(toric_variety(ctx), Int[m[i] for i in 1:rank(parent(m))]; D=get_chamber_dict(ctx))
    k = -minimum(hcat(rationoms...); init=-1)
    n = ngens(irrelevant_ideal(toric_variety(ctx)))
    return Int[k for i in 1:n]
  end
  # We use [CLS11](@cite), Lemma 9.5.8 and Theorem 9.5.10 for this.
  # TODO: Remove this?
  p, q = _proportionality_factors(ctx)
  G = grading_group(graded_ring(ctx))
  k0, rem = divrem(p*maximum([Int(abs(d[i])) for i in 1:ngens(G)]), q)
  if !is_zero(rem)
    k0 += 1
  end
  return Int[k0 for _ in 1:length(cech_complex_generators(ctx))]
end

function get_chamber_dict(ctx::ToricCtx)
  if !isdefined(ctx, :D)
    ctx.D = Dict()
  end
  return ctx.D
end

# See [CLS11](@cite), Theorem 9.5.10
# TODO: Remove this?
function _proportionality_factors(ctx::ToricCtx)
  if isnothing(ctx.proportionality_factors)
    S = graded_ring(ctx)
    G = grading_group(S)
    X = toric_variety(ctx)
    A = matrix(ZZ, rays(X))
    n = ncols(A)
    q_n = minimum([abs(a) for a in minors(A, n) if !is_zero(a)])
    Q_1 = maximum([abs(a) for a in A])
    Q_n1 = maximum([abs(a) for a in minors(A, n-1)])
    p = n^2*Q_1*Q_n1
    q = q_n
    ctx.proportionality_factors = (Int(p), Int(q))
  end
  return ctx.proportionality_factors::Tuple{Int, Int}
end

function getindex(ctx::ToricCtx, alpha::Vector{Int}, beta::Vector{Int})
  @assert length(alpha) == length(beta) == length(cech_complex_generators(ctx))
  @assert all(a <= b for (a, b) in zip(alpha, beta))
  if ctx.algorithm == :ext
    return get!(ctx.inclusions, (alpha, beta)) do
      S = graded_ring(ctx)
      c_alpha = ctx[alpha]::HomComplex
      c_beta = ctx[beta]::HomComplex
      a_dom = domain(c_alpha)
      b_dom = domain(c_beta)
      init_map = hom(b_dom[0], a_dom[0], [x^(beta[k]-alpha[k])*a for (k, x, a) in zip(1:length(alpha), cech_complex_generators(ctx), gens(a_dom[0]))])
      lift = lift_map(b_dom, a_dom, init_map)
      return hom(lift, codomain(c_alpha); domain=c_alpha, codomain=c_beta)
    end
  elseif ctx.algorithm == :cech
    return get!(ctx.inclusions, (alpha, beta)) do
      S = graded_ring(ctx)
      c_alpha = ctx[alpha]::HomComplex
      c_beta = ctx[beta]::HomComplex
      a_dom = domain(c_alpha)
      b_dom = domain(c_beta)
      a_unshift = original_complex(a_dom::ShiftedHyperComplex)
      b_unshift = original_complex(b_dom::ShiftedHyperComplex)
      a_kosz = original_complex(a_unshift::HyperComplexView)
      b_kosz = original_complex(b_unshift::HyperComplexView)
      a_seq = sequence(a_kosz)
      b_seq = sequence(b_kosz)
      n = length(cech_complex_generators(ctx))
      trans_mat = sparse_matrix(S, 0, n)
      for (i, (p, q)) in enumerate(zip(a_seq, b_seq))
        push!(trans_mat, sparse_row(S, [(i, divexact(q, p))]))
      end
      ind_kosz = Oscar.InducedKoszulMorphism(b_kosz, a_kosz; transition_matrix=trans_mat)
      ind_kosz_trunc = ind_kosz[1:n, domain=b_unshift, codomain=a_unshift]
      ind_kosz_shift = shift(ind_kosz_trunc, [1]; domain=b_dom, codomain=a_dom)
      hom(ind_kosz_shift, codomain(c_alpha); domain=c_alpha, codomain=c_beta)
    end
  end
  error("algorithm not recognized")
end

function getindex(ctx::ToricCtx, alpha::Vector{Int}, beta::Vector{Int}, d::FinGenAbGroupElem)
  @assert length(alpha) == length(beta) == length(cech_complex_generators(ctx))
  if all(a <= b for (a, b) in zip(alpha, beta))
    return get!(ctx.strand_inclusions, (alpha, beta, d)) do 
      strand(ctx[alpha, beta], d; domain=ctx[alpha, d], codomain=ctx[beta, d])
    end
  elseif all(a >= b for (a, b) in zip(alpha, beta)) 
    if ctx.algorithm == :cech
      return get!(ctx.strand_projections, (alpha, beta, d)) do
        SummandProjection(ctx[beta, alpha, d])
      end
    elseif ctx.algorithm == :ext
      return get!(ctx.strand_projections, (alpha, beta, d)) do
        HomologyPseudoInverse(ctx[beta, alpha, d])
      end
    else
      error("algorithm not recognized")
    end
  elseif ctx.algorithm == :cech
    c = Int[max(a, b) for (a, b) in zip(alpha, beta)]
    return compose(ctx[alpha, c, d], ctx[c, beta, d])
  end
  error("the given constellation of exponent vectors can not be handled with the chosen algorithm")
end


########################################################################
# A context object for computing the spectral sequence associated 
# to a the ̌Cech double complex on a toric variety with parameters
#
# Let X be a toric variety with (multi-)graded Cox Ring S over ℚ. 
# For a ℚ-algebra R, say R = ℚ[a₁,…, aₙ], one can form the ring 
# S' = S ⊗ R (tensor product taken over ℚ). The natural map R → S' 
# then corresponds to the projection X × Spec R → Spec R. 
#
# The `ToricCtxWithParams` is a context object to allow for 
# computation of cohomology of pullbacks of toric line bundles 
# ρ* ℒ along the projection ρ : X × Spec R → X and induced 
# morphisms 
# 
#    ϕ : ρ* ℒ → ρ* ℒ'
#
# defined on X × Spec R (rather than just X). It is clear that the 
# cohomology computations for ρ* ℒ should be obtained by base 
# change from ℒ. The induced maps in cohomology, however, require 
# coefficients in R. This leads to a quite convoluted caching and 
# transformation pattern. 
########################################################################
mutable struct ToricCtxWithParams
  pure_ctx::ToricCtx
  R::Ring
  truncated_cech_complexes::Dict{Vector{Int}, AbsHyperComplex}
  inclusions::Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}
  projections::Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}
  strands::Dict{Vector{Int}, Dict}
  simplified_strands::IdDict{AbsHyperComplex, AbsHyperComplex}
  simplified_strands_to_orig::IdDict{Map, Map}
  simplified_strands_from_orig::IdDict{Map, Map}
  simplified_strands_homotopy::IdDict{Map, Map}
  strand_inclusions::Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}
  strand_projections::Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}
  induced_cohomology_maps::Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem, Int}, Map}
  cohomology_models::Dict{FinGenAbGroupElem, AbsHyperComplex}
  cohomology_inclusions::Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}
  cohomology_projections::Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}
  # mult_map_cache::Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, Dict}
  mult_map_cache::Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, WeakKeyDict}
  transfer::Map

  function ToricCtxWithParams(X::NormalToricVariety, transfer::Map; algorithm::Symbol=:ext)
    pure_ctx = ToricCtx(X; algorithm)
    @assert domain(transfer) === cox_ring(X)
    R = coefficient_ring(codomain(transfer))
    return new(pure_ctx, R, 
               Dict{Vector{Int}, AbsHyperComplex}(),
               Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}(),
               Dict{Tuple{Vector{Int}, Vector{Int}}, AbsHyperComplexMorphism}(),
               Dict{Vector{Int}, Dict}(),
               IdDict{AbsHyperComplex, AbsHyperComplex}(),
               IdDict{Map, Map}(),
               IdDict{Map, Map}(),
               IdDict{Map, Map}(),
               Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}(),
               Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem}, AbsHyperComplexMorphism}(), 
               Dict{Tuple{Vector{Int}, Vector{Int}, FinGenAbGroupElem, Int}, Map}(), 
               Dict{FinGenAbGroupElem, AbsHyperComplex}(),
               Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}(),
               Dict{Tuple{FinGenAbGroupElem, Vector{Int}}, AbsHyperComplexMorphism}(),
               # Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, Dict}()
               Dict{Tuple{Vector{Int}, FinGenAbGroupElem, Int}, WeakKeyDict}(), 
               transfer
              )
  end
end

toric_variety(ctx::ToricCtxWithParams) = toric_variety(ctx.pure_ctx)
graded_ring(ctx::ToricCtxWithParams) = codomain(ctx.transfer)

function getindex(ctx::ToricCtxWithParams, alpha::Vector{Int})
  return get!(ctx.truncated_cech_complexes, alpha) do
    res, tr = change_base_ring(ctx.transfer, ctx.pure_ctx[alpha])
    return res
  end
end

function getindex(ctx::ToricCtxWithParams, alpha::Vector{Int}, d::FinGenAbGroupElem)
  strands = get!(ctx.strands, alpha) do 
    Dict{typeof(d), AbsHyperComplex}()
  end
  return get!(strands, d) do
    res, tr = change_base_ring(ctx.R, ctx.pure_ctx[alpha, d])
    res
  end
end

function simplified_strand(ctx::ToricCtxWithParams, 
    alpha::Vector{Int}, d::FinGenAbGroupElem
  )
  str = simplified_strand(ctx.pure_ctx, alpha, d)
  return get!(ctx.simplified_strands, str) do
    change_base_ring(ctx.R, str)[1]
  end
end

function simplified_strand_homotopy(
    ctx::ToricCtxWithParams, alpha::Vector{Int}, d::FinGenAbGroupElem, p::Int
  )
  h = simplified_strand_homotopy(ctx.pure_ctx, alpha, d, p)
  return get!(ctx.simplified_strands_homotopy, h) do
    change_base_ring(ctx.R, h; 
                     domain=ctx[alpha, d][p], 
                     codomain=ctx[alpha, d][p+1]
                    )[1]
  end
end

function simplified_strand_inclusion(
    ctx::ToricCtxWithParams, alpha::Vector{Int}, d::FinGenAbGroupElem, p::Int
  )
  inc = simplified_strand_inclusion(ctx.pure_ctx, alpha, d, p)
  res = get!(ctx.simplified_strands_to_orig, inc) do
    change_base_ring(ctx.R, inc; 
                     domain=simplified_strand(ctx, alpha, d)[p],
                     codomain=ctx[alpha, d][p]
                    )[1]
  end
  @assert domain(res) === simplified_strand(ctx, alpha, d)[p]
  @assert codomain(res) === ctx[alpha, d][p]
  return res
end

function simplified_strand_projection(
    ctx::ToricCtxWithParams, alpha::Vector{Int}, d::FinGenAbGroupElem, p::Int
  )
  h = simplified_strand_projection(ctx.pure_ctx, alpha, d, p)
  return get!(ctx.simplified_strands_from_orig, h) do
    change_base_ring(ctx.R, h;
                     domain=ctx[alpha, d][p],
                     codomain=simplified_strand(ctx, alpha, d)[p]
                    )[1]
  end
end

function induced_cohomology_map(
    ctx::ToricCtxWithParams, e0::Vector{Int},
    e1::Vector{Int}, d::FinGenAbGroupElem, 
    i::Int
  )
  return get!(ctx.induced_cohomology_maps, (e0, e1, d, i)) do 
    change_base_ring(ctx.R, 
                     induced_cohomology_map(ctx.pure_ctx, e0, e1, d, i);
                     domain=simplified_strand(ctx, e0, d)[i],
                     codomain=simplified_strand(ctx, e1, d)[i]
                    )[1]
  end
end

function cohomology_model(ctx::ToricCtxWithParams, d::FinGenAbGroupElem)
  get!(ctx.cohomology_models, d) do
    simplified_strand(ctx, _minimal_exponent_vector(ctx, d), d)
  end
end

function cohomology_model_inclusion(ctx::ToricCtxWithParams, d::FinGenAbGroupElem, i::Int)
  h = cohomology_model(ctx, d)
  c = ctx[_minimal_exponent_vector(ctx.pure_ctx, d), d]
  to_orig = map_to_original_complex(cohomology_model(ctx.pure_ctx, d))[i]
  res, _, _ = change_base_ring(ctx.R, to_orig; domain=h[i], codomain=c[i])
  return res
end

function cohomology_model_projection(ctx::ToricCtxWithParams, d::FinGenAbGroupElem, i::Int)
  h = cohomology_model(ctx, d)
  c = ctx[_minimal_exponent_vector(ctx.pure_ctx, d), d]
  from_orig = map_from_original_complex(cohomology_model(ctx.pure_ctx, d))[i]
  res, _, _ = change_base_ring(ctx.R, from_orig; domain=c[i], codomain=h[i])
  return res
end

function getindex(ctx::ToricCtxWithParams, alpha::Vector{Int}, beta::Vector{Int})
  return get!(ctx.inclusions, (alpha, beta)) do
    pure = ctx.pure_ctx[alpha, beta]
    res, _, _ = change_base_ring(ctx.transfer, pure; domain=ctx[alpha], codomain=ctx[beta])
    res
  end
end

function getindex(ctx::ToricCtxWithParams, alpha::Vector{Int}, beta::Vector{Int}, d::FinGenAbGroupElem)
  if all(a <= b for (a, b) in zip(alpha, beta))
    return get!(ctx.strand_inclusions, (alpha, beta, d)) do 
      res, _, _ = change_base_ring(ctx.R, ctx.pure_ctx[alpha, beta, d]; domain=ctx[alpha, d], codomain=ctx[beta, d])
      res
    end
  elseif all(a >= b for (a, b) in zip(alpha, beta)) 
    # TODO: Make the stuff below work with base change, too.
    if ctx.pure_ctx.algorithm == :cech
      return get!(ctx.strand_projections, (alpha, beta, d)) do
        SummandProjection(ctx[beta, alpha, d])
      end
    elseif ctx.pure_ctx.algorithm == :ext
      return get!(ctx.strand_projections, (alpha, beta, d)) do
        HomologyPseudoInverse(ctx[beta, alpha, d])
      end
    else
      error("algorithm not recognized")
    end
  elseif ctx.pure_ctx.algorithm == :cech
    c = Int[max(a, b) for (a, b) in zip(alpha, beta)]
    return compose(ctx[alpha, c, d], ctx[c, beta, d])
  end
  error("the given constellation of exponent vectors can not be handled with the chosen algorithm")
end

function change_base_ring(bc::Any, f::ModuleFPHom{DT, CT, Nothing}; 
    domain::ModuleFP=change_base_ring(bc, Oscar.domain(f))[1],
    codomain::ModuleFP=change_base_ring(bc, Oscar.codomain(f))[1]
  ) where {DT, CT}
  img_gens = elem_type(codomain)
  bc_dom = hom(Oscar.domain(f), domain, gens(domain), bc)
  bc_cod = hom(Oscar.codomain(f), codomain, gens(codomain), bc)
  res = hom(domain, codomain, bc_cod.(images_of_generators(f)))
  return res, bc_dom, bc_cod
end

function _minimal_exponent_vector(ctx::ToricCtxWithParams, m::FinGenAbGroupElem)
  return _minimal_exponent_vector(ctx.pure_ctx, m)
end

# outer constructor
@attr ToricCtx function local_cohomology_context_object(X::NormalToricVariety; algorithm::Symbol=:ext)
  return ToricCtx(X; algorithm)
end

