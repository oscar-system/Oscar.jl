@doc raw"""
    FPGroup(W::WeylGroup) -> FPGroup
    fp_group(W::WeylGroup) -> FPGroup

Construct a group of type `FPGroup` that is isomorphic to `W`.

If one needs the isomorphism then [`isomorphism(::Type{FPGroup}, W::WeylGroup)`](@ref)
can be used instead.
"""
function FPGroup(W::WeylGroup)
  return codomain(isomorphism(FPGroup, W))
end

@doc raw"""
    PermGroup(W::WeylGroup) -> PermGroup
    permutation_group(W::WeylGroup) -> PermGroup

Construct a group of type `PermGroup` that is isomorphic to `W`.

If one needs the isomorphism then [`isomorphism(::Type{PermGroup}, W::WeylGroup)`](@ref)
can be used instead.
"""
function PermGroup(W::WeylGroup)
  return codomain(isomorphism(PermGroup, W))
end

fp_group(W::WeylGroup) = FPGroup(W)
permutation_group(W::WeylGroup) = PermGroup(W)

@doc raw"""
    isomorphism(::Type{FPGroup}, W::WeylGroup) -> Map{WeylGroup, FPGroup}

Return an isomorphism from `W` to a group `H` of type `FPGroup`.

`H` will be the quotient of a free group with the same rank as `W`,
where we have the natural 1-to-1 correspondence of generators, modulo the Coxeter relations of `W`.

Isomorphisms are cached in `W`, subsequent calls of `isomorphism(FPGroup, W)` yield identical results.

If only the image of such an isomorphism is needed, use `fp_group(W)`.
"""
function isomorphism(T::Type{FPGroup}, W::WeylGroup; on_gens::Bool=true)
  on_gens = true # we ignore the on_gens flag, the iso will *always* map gens onto gens
  isos =
    get_attribute!(Dict{Tuple{Type,Bool},Any}, W, :isomorphisms)::Dict{Tuple{Type,Bool},Any}
  return get!(isos, (T, on_gens)) do
    G = _isomorphic_group_on_gens(T, W)

    # help GAP a bit
    set_is_finite(G, is_finite(W))
    is_finite(W) && set_order(G, order(W))

    iso = function (w::WeylGroupElem)
      return G(syllables(w)) # TODO: change to letters once G supports that input
    end

    isoinv = function (g::FPGroupElem)
      return W(abs.(letters(g)))
    end

    return MapFromFunc(W, G, iso, isoinv)
  end::MapFromFunc{WeylGroup,T}
end

# Constructs the same object as `fp_group(W) === codomain(isomorphism(FPGroup, W))`
# but without setting GAP attributes, creating the iso functions, and caching.
# We use this function during testing as setting GAP attributes may skip
# some computations and thus make the tests less meaningful.
function _isomorphic_group_on_gens(::Type{FPGroup}, W::WeylGroup)
  R = root_system(W)
  F = free_group(rank(R))

  gcm = cartan_matrix(R)
  rels = [
    (gen(F, i) * gen(F, j))^coxeter_matrix_entry_from_cartan_matrix(gcm, i, j) for
    i in 1:rank(R) for j in i:rank(R)
  ]

  G, _ = quo(F, rels)
  return G
end

@doc raw"""
    isomorphism(::Type{PermGroup}, W::WeylGroup) -> Map{WeylGroup, PermGroup}

Return an isomorphism from `W` to a group `H` of type `PermGroup`.
An exception is thrown if no such isomorphism exists.

The generators of `H` are in the natural 1-1 correspondence with the generators of `W`.

If the type of `W` is irreducible and not $E_6$ or $E_7$, then the degree of `H` is optimal.
See [Sau14](@cite) for the optimal permutation degrees of Weyl groups.

Isomorphisms are cached in `W`, subsequent calls of `isomorphism(PermGroup, W)` yield identical results.

If only the image of such an isomorphism is needed, use `permutation_group(W)`.
"""
function isomorphism(T::Type{PermGroup}, W::WeylGroup; on_gens::Bool=true)
  on_gens = true # we ignore the on_gens flag, the iso will *always* map gens onto gens
  isos =
    get_attribute!(Dict{Tuple{Type,Bool},Any}, W, :isomorphisms)::Dict{Tuple{Type,Bool},Any}
  return get!(isos, (T, on_gens)) do
    G = _isomorphic_group_on_gens(T, W)

    # help GAP a bit
    set_order(G, order(W))

    # Create a homomorphism mapping gens(W) to gen_G.
    # This is a workaround until hom works for Weyl groups.
    epi = epimorphism_from_free_group(G)

    iso = function (w::WeylGroupElem)
      map_word(w, gens(G); init=one(G))
    end

    isoinv = function (p::PermGroupElem)
      rep_word = abs.(word(preimage(epi, p))) # `abs` may be used as all gens of W are self-inverse
      return W(rep_word)
    end

    return MapFromFunc(W, G, iso, isoinv)
  end::MapFromFunc{WeylGroup,T}
end

# Constructs the same object as `permutation_group(W) === codomain(isomorphism(PermGroup, W))`
# but without setting GAP attributes, creating the iso functions, and caching.
# We use this function during testing as setting GAP attributes may skip
# some computations and thus make the tests less meaningful.
function _isomorphic_group_on_gens(::Type{PermGroup}, W::WeylGroup)
  @req is_finite(W) "Weyl group is not finite"
  R = root_system(W)
  type, ordering = root_system_type_with_ordering(R)

  if length(type) != 1
    # Apply recursion to the irreducible factors
    # Note that the gens of each irreducible factor are in canonical ordering
    factors = irreducible_factors(W; morphisms=false)
    factors_as_perm = [_isomorphic_group_on_gens(PermGroup, factor) for factor in factors]
    G = inner_direct_product(factors_as_perm; morphisms=false)
    # gens(G) corresponds to gens(W)[ordering], so
    # gens(G)[sortperm(ordering)] corresponds to gens(W)
    G_sorted, _ = sub(G, gens(G)[sortperm(ordering)])
    return G_sorted
  end

  # Compute generators of the permutation group to which the simple reflections are mapped.
  # These generators correspond to the sorted ordering. They will later be reordered.
  coxeter_type, n = only(type)
  if coxeter_type == :A
    Sym = symmetric_group(n + 1)
    gen_G = [cperm(Sym, [i, i + 1]) for i in 1:n]
  elseif coxeter_type == :B || coxeter_type == :C
    Sym = symmetric_group(2n)
    gen_G = vcat(
      [cperm(Sym, [i, i + 1], [i + n, i + 1 + n]) for i in 1:(n - 1)], cperm(Sym, [n, 2n])
    )
  elseif coxeter_type == :D
    Sym = symmetric_group(2n)
    gen_G = vcat(
      [cperm(Sym, [i, i + 1], [i + n, i + 1 + n]) for i in 1:(n - 1)],
      cperm(Sym, [n - 1, 2n], [n, 2n - 1]),
    )
  elseif coxeter_type == :E
    # Permutation representation on the root system. The permutation degree is not optimal for E_6 and E_7.
    m = number_of_roots(R)
    Sym = symmetric_group(m)
    gen_G = [
      perm(Sym, [is_root_with_index(reflect(root(R, j), i))[2] for j in 1:m]) for i in 1:n
    ]
  elseif coxeter_type == :F
    Sym = symmetric_group(24)
    #! format: off
    gen_G = [ # Computed by hand
      perm(Sym, [1, 2, 5, 4, 3, 7, 6, 9, 8, 10, 11, 12, 13, 14, 17, 16, 15, 19, 18, 21, 20, 22, 23, 24]),
      perm(Sym, [3, 2, 1, 6, 5, 4, 7, 8, 10, 9, 11, 12, 15, 14, 13, 18, 17, 16, 19, 20, 22, 21, 23, 24]),
      perm(Sym, [13, 4, 3, 2, 5, 8, 9, 6, 7, 11, 10, 12, 1, 16, 15, 14, 17, 20, 21, 18, 19, 23, 22, 24]),
      perm(Sym, [4, 14, 6, 1, 7, 3, 5, 8, 9, 10, 12, 11, 16, 2, 18, 13, 19, 15, 17, 20, 21, 22, 24, 23]),
    ]
    #! format: on
  elseif coxeter_type == :G
    Sym = symmetric_group(5)
    gen_G = [cperm(Sym, [1, 2], [3, 5]), cperm(Sym, [4, 5])] # gen_G[1]*gen_G[2] = cperm([1,2], [3,4,5])
  end

  # Reorder generators
  # (Details: simple_roots(R)[ordering] is in canonical ordering.
  # s = sortperm(ordering) is the inverse of the corresponding permutation.
  # Hence (gen_G[s])[i] is the image of gen(W, i).)
  gen_G = gen_G[sortperm(ordering)]
  G, _ = sub(Sym, gen_G)
  return G
end

@doc raw"""
    parabolic_subgroup(W::WeylGroup, vec::Vector{<:Integer}, w::WeylGroupElem=one(W)) -> WeylGroup, Map{WeylGroup, WeylGroup}

Return a Weyl group `P` and an embedding $f:P\to W$ such that $f(P)$ is
the subgroup `U` of `W` generated by `[inv(w)*u*w for u in gens(W)[vec]]`.
Further, `f` maps `gen(P, i)` to `inv(w)*gen(W, vec[i])*w`.
The elements of `vec` must be pairwise distinct integers in
`1:number_of_generators(W)` and `vec` must be non-empty.

# Examples
```jldoctest
julia> W = weyl_group(:B, 3)
Weyl group
  of root system of rank 3
    of type B3

julia> P1, f1 = parabolic_subgroup(W, [1, 2])
(Weyl group of root system of type A2, Map: P1 -> W)

julia> f1(P1[1] * P1[2]) == W[1] * W[2]
true

julia> P2, f2 = parabolic_subgroup(W, [1, 2], W[1])
(Weyl group of root system of type A2, Map: P2 -> W)

julia> f2(P2[1]) == W[1] && f2(P2[2]) == W[1] * W[2] * W[1]
true

julia> P3, f3 = parabolic_subgroup(W, [1,3,2])
(Weyl group of root system of type B3 (non-canonical ordering), Map: P3 -> W)

julia> f3(P3[2]) == W[3]
true
```
"""
function parabolic_subgroup(W::WeylGroup, vec::Vector{<:Integer}, w::WeylGroupElem=one(W))
  @req allunique(vec) "Elements of vector are not pairwise distinct"
  @req all(i -> 1 <= i <= number_of_generators(W), vec) "Invalid indices"
  cm = cartan_matrix(W)[vec, vec]
  para = weyl_group(cm)
  genimgs = [conj(W[i], w) for i in vec]
  emb = function (u::WeylGroupElem)
    return map_word(u, genimgs)
  end
  return para, MapFromFunc(para, W, emb)
end

@doc raw"""
    parabolic_subgroup_with_projection(W::WeylGroup, vec::Vector{<:Integer}; check::Bool=true) -> WeylGroup, Map{WeylGroup, WeylGroup}, Map{WeylGroup, WeylGroup}

Return a triple `(P, emb, proj)` that describes a factor of `W`, that is,
a product of irreducible factors.
Here `P, emb = `[`parabolic_subgroup`](@ref)`(W, vec)`
and `proj` is the projection map from `W` onto `P`,
which is a left-inverse of `emb`.

If `check = true`, then it is checked whether `vec` actually describes
a union of irreducible components of the Dynkin diagram.

# Examples
```jldoctest
julia> W = weyl_group([(:A, 3), (:B, 3)])
Weyl group
  of root system of rank 6
    of type A3 x B3

julia> P1, f1, p1 = parabolic_subgroup_with_projection(W, [1,2,3])
(Weyl group of root system of type A3, Map: P1 -> W, Map: W -> P1)

julia> p1(W[1]*W[4]*W[2]*W[6]) == P1[1] * P1[2]
true

julia> P2, f2, p2 = parabolic_subgroup_with_projection(W, [4,6,5])
(Weyl group of root system of type B3 (non-canonical ordering), Map: P2 -> W, Map: W -> P2)

julia> p2(W[5]) == P2[3]
true
```
"""
function parabolic_subgroup_with_projection(
  W::WeylGroup, vec::Vector{<:Integer}; check::Bool=true
)
  if check
    # Check that every generator in gens(W)[vec] commutes with every other generator.
    # In other words, vec describes a union of irreducible components of the Coxeter diagram.
    cm = cartan_matrix(W)
    for i in setdiff(1:number_of_generators(W), vec)
      for j in vec
        @req is_zero_entry(cm, i, j) begin
          "Input vector must describe a direct factor of the Weyl group"
        end
      end
    end
  end

  factor, emb = parabolic_subgroup(W, vec)
  # Generators of W are mapped to the corresponding generators of factor,
  # or to 1 if there is no corresponding generator
  proj_gen_imgs = fill(one(factor), ngens(W))
  for i in 1:length(vec)
    proj_gen_imgs[vec[i]] = gen(factor, i)
  end
  proj = function (w::WeylGroupElem)
    return map_word(w, proj_gen_imgs)
  end
  return factor, emb, MapFromFunc(W, factor, proj)
end

@doc raw"""
    irreducible_factors(W::WeylGroup; morphisms::Bool=false)

If `morphisms` is `true`, return a triple `(U, emb, proj)` that describes the irreducible factors of `W`.
That is, `U`, `emb`, `proj` are vectors whose length is
the number of irreducible components of the Dynkin diagram of `W`, and for each `i`,
`(U[i], emb[i], proj[i])` describes the `i`-th factor of `W` as described in
[`parabolic_subgroup_with_projection`](@ref).

If `morphisms` is `false`, return only `U` .

The order of the irreducible factors is the one given by [`cartan_type_with_ordering`](@ref).

See also [`inner_direct_product(::AbstractVector{WeylGroup})`](@ref).

# Examples
```jldoctest
julia> W = weyl_group([(:A, 3), (:B, 4)]);

julia> irreducible_factors(W)
2-element Vector{WeylGroup}:
 Weyl group of root system of type A3
 Weyl group of root system of type B4

julia> U, emb, proj = irreducible_factors(W; morphisms=true)
(WeylGroup[Weyl group of root system of type A3, Weyl group of root system of type B4], Map{WeylGroup, WeylGroup}[Map: Weyl group -> W, Map: Weyl group -> W], Map{WeylGroup, WeylGroup}[Map: W -> Weyl group, Map: W -> Weyl group])

julia> emb[1](U[1][1]) == W[1]
true

julia> proj[2](W[4]) == U[2][1]
true
```
"""
function irreducible_factors(W::WeylGroup; morphisms::Bool=false)
  types, ordering = root_system_type_with_ordering(root_system(W))
  num_factors = length(types) # Number of irreducible factors of W

  U = Vector{WeylGroup}(undef, num_factors)
  emb = Vector{Map{WeylGroup,WeylGroup}}(undef, num_factors)
  proj = Vector{Map{WeylGroup,WeylGroup}}(undef, num_factors)

  start_index = 1
  for (i, type) in enumerate(types)
    rk = type[2] # rank of the i-th irreducible factor
    end_index = start_index + rk - 1
    U[i], emb[i], proj[i] = parabolic_subgroup_with_projection(
      W, ordering[start_index:end_index]; check=false # No check needed by construction
    )
    start_index = end_index + 1
  end

  if morphisms
    return U, emb, proj
  else
    return U
  end
end

@doc raw"""
    inner_direct_product(L::AbstractVector{WeylGroup}; morphisms::Bool=false)
    inner_direct_product(L::WeylGroup...; morphisms::Bool=false)

If `morphisms` is `false`, then return a Weyl group `W` that is isomorphic to the
direct product of the Weyl groups in `L`.
The generators of `W` are in natural 1-1 correspondence with the generators of the groups in `L`
in the given order.

If `morphisms` is `true`, return a triple `(W, emb, proj)` where
`W` is as above and `emb`, `proj` are vectors such that
`emb[i]` (resp., `proj[i]`) is the embedding of `L[i]` into `W`
(resp., the projection of `W` onto `L[i]`).

See also [`inner_direct_product(::AbstractVector{T}) where {T<:Union{PcGroup, SubPcGroup, FPGroup, SubFPGroup}}`](@ref).

# Examples
```jldoctest
julia> W1 = weyl_group(:A, 2); W2 = weyl_group(:B, 3);

julia> W, emb, proj = inner_direct_product(W1, W2; morphisms=true)
(Weyl group of root system of type A2 x B3, Map{WeylGroup, WeylGroup}[Map: W1 -> W, Map: W2 -> W], Map{WeylGroup, WeylGroup}[Map: W -> W1, Map: W -> W2])

julia> proj[2](W[3]) == W2[1]
true

julia> gens(W) == vcat(emb[1].(gens(W1)), emb[2].(gens(W2)))
true
```
"""
function inner_direct_product(L::AbstractVector{WeylGroup}; morphisms::Bool=false)
  if length(L) > 0
    cm_blocks = [cartan_matrix(factor) for factor in L]
    cm = block_diagonal_matrix(cm_blocks)
  else
    cm = zero_matrix(ZZ, 0, 0)
  end
  product = weyl_group(cm)

  if !morphisms
    return product
  end

  emb = Vector{Map{WeylGroup,WeylGroup}}(undef, length(L))
  proj = Vector{Map{WeylGroup,WeylGroup}}(undef, length(L))
  start_index = 1
  for (i, factor) in enumerate(L)
    end_index = start_index + ngens(factor) - 1
    # Embedding
    emb_gen_imgs = [product[j] for j in start_index:end_index]
    emb[i] = MapFromFunc(
      factor, product, w::WeylGroupElem -> map_word(w, emb_gen_imgs; init=one(product))
    )
    # Projection
    proj_gen_imgs = fill(one(factor), ngens(product))
    proj_gen_imgs[start_index:end_index] = gens(factor)
    proj[i] = MapFromFunc(
      product, factor, w::WeylGroupElem -> map_word(w, proj_gen_imgs; init=one(factor))
    )
    start_index = end_index + 1
  end

  return product, emb, proj
end

function inner_direct_product(L::WeylGroup, Ls::WeylGroup...; morphisms::Bool=false)
  return inner_direct_product([L, Ls...]; morphisms=morphisms)
end
