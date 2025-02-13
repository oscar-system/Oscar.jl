function coxeter_matrix(W::WeylGroup)
  return cartan_to_coxeter_matrix(cartan_matrix(root_system(W)))
end

@doc raw"""
    fp_group(W::WeylGroup) -> FPGroup

Construct a group of type `FPGroup` that is isomorphic to `W`.

The `FPGroup` will be the quotient of a free group with the same rank as `W`,
where we have the natural 1-to-1 correspondence of generators, modulo the Coxeter relations of `W`.

Also see: [`isomorphism(::Type{FPGroup}, ::WeylGroup)`](@ref).
"""
function fp_group(W::WeylGroup; set_properties::Bool=true)
  return codomain(isomorphism(FPGroup, W; set_properties))
end

@doc raw"""
    isomorphism(::Type{FPGroup}, W::WeylGroup) -> Map{WeylGroup, FPGroup}

Construct an isomorphism between `W` and a group of type `FPGroup`.

The properties of the codomain group and the isomorphism are described in [`fp_group(::WeylGroup)`](@ref).
"""
function isomorphism(::Type{FPGroup}, W::WeylGroup; set_properties::Bool=true)
  R = root_system(W)
  F = free_group(rank(R))

  gcm = cartan_matrix(R)
  rels = [
    (gen(F, i) * gen(F, j))^coxeter_matrix_entry_from_cartan_matrix(gcm, i, j) for
    i in 1:rank(R) for j in i:rank(R)
  ]

  G, _ = quo(F, rels)

  if set_properties
    set_is_finite(G, is_finite(W))
    is_finite(W) && set_order(G, order(W))
  end

  iso = function (w::WeylGroupElem)
    return G([i => 1 for i in word(w)])
  end

  isoinv = function (g::FPGroupElem)
    return W(abs.(letters(g))) # TODO: check if normalize=false can be added here (probably not)
  end

  return MapFromFunc(W, G, iso, isoinv)
end

@doc raw"""
    permutation_group(W::WeylGroup) -> PermGroup

Construct a group of type `PermGroup` that is isomorphic to `W`. The generators of the `PermGroup` are in
the natural 1-1 correspondence with the generators of `W`.

If the type of `W` is irreducible and not $E_6$ or $E_7$, then the degree of the constructed `PermGroup` is optimal.
See [Sau14](@cite) for the optimal permutation degrees of Weyl groups.

Also see: [`isomorphism(::Type{PermGroup}, ::WeylGroup)`](@ref).
"""
function permutation_group(W::WeylGroup; set_properties::Bool=true)
  return codomain(isomorphism(PermGroup, W; set_properties))
end

@doc raw"""
    isomorphism(::Type{PermGroup}, W::WeylGroup) -> Map{WeylGroup, PermGroup}

Construct an isomorphism between `W` and a group of type `PermGroup`.

The properties of the codomain group and the isomorphism are described in [`permutation_group(::WeylGroup)`](@ref).
"""
function isomorphism(::Type{PermGroup}, W::WeylGroup; set_properties::Bool=true)
  @req is_finite(W) "Weyl group is not finite"
  R = root_system(W)
  type, ordering = root_system_type_with_ordering(R)

  if length(type) != 1
    error("Not implemented (yet)")
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
    gen_G = [ # Computed by hand
      perm(
        Sym,
        [
          1,
          2,
          5,
          4,
          3,
          7,
          6,
          9,
          8,
          10,
          11,
          12,
          13,
          14,
          17,
          16,
          15,
          19,
          18,
          21,
          20,
          22,
          23,
          24,
        ],
      ),
      perm(
        Sym,
        [
          3,
          2,
          1,
          6,
          5,
          4,
          7,
          8,
          10,
          9,
          11,
          12,
          15,
          14,
          13,
          18,
          17,
          16,
          19,
          20,
          22,
          21,
          23,
          24,
        ],
      ),
      perm(
        Sym,
        [
          13,
          4,
          3,
          2,
          5,
          8,
          9,
          6,
          7,
          11,
          10,
          12,
          1,
          16,
          15,
          14,
          17,
          20,
          21,
          18,
          19,
          23,
          22,
          24,
        ],
      ),
      perm(
        Sym,
        [
          4,
          14,
          6,
          1,
          7,
          3,
          5,
          8,
          9,
          10,
          12,
          11,
          16,
          2,
          18,
          13,
          19,
          15,
          17,
          20,
          21,
          22,
          24,
          23,
        ],
      ),
    ]
  elseif coxeter_type == :G
    Sym = symmetric_group(5)
    gen_G = [cperm(Sym, [1, 2], [3, 5]), cperm(Sym, [4, 5])] # gen_G[1]*gen_G[2] = cperm([1,2], [3,4,5])
  end

  # Reorder generators
  # (Details: simple_roots(R)[ordering] is in canonical ordering.
  # s = sortperm(ordering) is the inverse of the corresponding permutation.
  # Hence (gen_G[s])[i] is the image of gen(W, i).)
  gen_G = gen_G[sortperm(ordering)]

  # Create a homomorphism mapping gens(W) to gen_G.
  # This is a workaround until hom works for Weyl groups.
  G, _ = sub(Sym, gen_G)
  epi = epimorphism_from_free_group(G)

  if set_properties
    set_order(G, order(W))
  end

  iso = function (w::WeylGroupElem)
    map_word(w, gens(G); init=one(G))
  end

  isoinv = function (p::PermGroupElem)
    rep_word = abs.(word(preimage(epi, p))) # `abs` may be used as all gens of W are self-inverse
    return W(rep_word)
  end

  return MapFromFunc(W, G, iso, isoinv)
end

@doc raw"""
    map_word(w::WeylGroupElem, genimgs::Vector; genimgs_inv::Vector = genimgs, init = nothing)

If `init` is `nothing` and `word(w) = [`$i_1$`, ..., `$i_n$`]`,
then return the product $R_1 R_2 \cdots R_n$ with $R_j =$ `genimgs[`$i_j$`]`.
Otherwise return the product $xR_1 R_2 \cdots R_n$ where $x =$ `init`.

The length of `genimgs` must be equal to the rank of the parent of `w`.
If `w` is the trivial element, then `init` is returned if it is different
from `nothing`, and otherwise `one(genimgs[1])` is returned if `genimgs` is non-empty.
If `w` is trivial, `init` is nothing and `genimgs` is empty, an error occurs.

See also: [`map_word(::Union{FPGroupElem, SubFPGroupElem}, ::Vector)`](@ref),
[`map_word(::Union{PcGroupElem, SubPcGroupElem}, ::Vector)`](@ref).
Note that `map_word(::WeylGroupElem)` accepts the `genimgs_inv` keyword argument
for consistency with other `map_word` methods, but ignores it because the
generators of a Weyl group are always self-inverse.

# Examples
```jldoctest
julia> W = weyl_group(:B, 3); imgs = [2, 3, 5];

julia> map_word(one(W), imgs)
1

julia> map_word(W([1]), imgs)
2

julia> map_word(W([1]), imgs; init=7)
14

julia> map_word(W([1,2,1]), imgs)
12

julia> map_word(W([2,1,2]), imgs) # W([2,1,2]) == W([1,2,1])
12

julia> map_word(W([3, 2, 1, 3, 2, 3]), imgs)
2250
```
"""
function map_word(
  w::WeylGroupElem, genimgs::Vector; genimgs_inv::Vector=genimgs, init=nothing
)
  @req length(genimgs) == number_of_generators(parent(w)) begin
    "Length of vector of images does not equal rank of Weyl group"
  end
  return map_word(Int.(word(w)), genimgs; init=init)
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
  cm = cartan_matrix(root_system(W))[vec, vec]
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
    cm = cartan_matrix(root_system(W))
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
