function _derived_pushforward(M::FreeMod)
  S = base_ring(M)
  n = ngens(S)-1

  d = _regularity_bound(M) - n
  d = (d < 0 ? 0 : d)

  Sd = graded_free_module(S, [0 for i in 1:ngens(S)])
  v = sum(x^d*Sd[i] for (i, x) in enumerate(gens(S)); init=zero(Sd))
  kosz = koszul_complex(Oscar.KoszulComplex, v)
  K = shift(Oscar.DegreeZeroComplex(kosz)[1:n+1], 1)

  KoM = hom(K, M)
  #KoM_simp, _, _ = simplify(KoM)
  st = strand(KoM, 0)[1]
  return st
end

@doc raw"""
    _derived_pushforward(M::SubquoModule)

We consider a graded module `M` over a standard graded polynomial ring 
``S = A[x₀,…,xₙ]`` as a representative of a coherent sheaf ``ℱ`` on 
relative projective space ``ℙ ⁿ_A``. Then we compute ``Rπ_* ℱ``` as a
complex of ``A``-modules where ``π : ℙ ⁿ_A → Spec(A)`` is the projection 
to the base. 
"""
function _derived_pushforward(M::SubquoModule)
  S = base_ring(M)
  n = ngens(S)-1

  d = _regularity_bound(M) - n
  d = (d < 0 ? 0 : d)

  Sd = graded_free_module(S, [0 for i in 1:ngens(S)])
  v = sum(x^d*Sd[i] for (i, x) in enumerate(gens(S)); init=zero(Sd))
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


_regularity_bound(F::FreeMod) = maximum(Int(degree(a; check=false)[1]) for a in gens(F))

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
        j = findfirst(k->k==real_j, J)
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

