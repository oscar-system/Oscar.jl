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
    reduce(*, (gen(G, Int(i)) for i in word(w)); init=one(G))
  end

  isoinv = function (p::PermGroupElem)
    rep_word = abs.(word(preimage(epi, p))) # `abs` may be used as all gens of W are self-inverse
    return W(rep_word)
  end

  return MapFromFunc(W, G, iso, isoinv)
end
