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

function permutation_group(W::WeylGroup; set_properties::Bool=true)
  return codomain(isomorphism(PermGroup, W; set_properties))
end

function isomorphism(::Type{PermGroup}, W::WeylGroup; set_properties::Bool=true)
  @req is_finite(W) "Weyl group is not finite"
  R = root_system(W)
  type, ordering = root_system_type_with_ordering(R)

  if length(type) != 1
    error("Not implemented (yet)")
  end
  if !issorted(ordering)
    error("Not implemented (yet)")
  end
  coxeter_type, n = only(type)
  if coxeter_type == :A
    G = symmetric_group(n + 1)

    iso = function (w::WeylGroupElem)
      reduce(*, [cperm(G, [i, i + 1]) for i in word(w)]; init=cperm(G))
    end

    isoinv = function (p::PermGroupElem)
      word = UInt8[]
      for cycle in cycles(p)
        transpositions = [
          sort([c, cycle[i + 1]]) for (i, c) in enumerate(cycle) if i < length(cycle)
        ]
        for t in transpositions
          word = reduce(
            vcat,
            [
              [i for i in t[1]:(t[2] - 1)],
              [i for i in reverse(t[1]:(t[2] - 2))],
              word,
            ],
          )
        end
      end
      return W(word) # TODO: check if normalize=false can be added here
    end
  else
    error("Not implemented (yet)")
  end

  if set_properties
    set_order(G, order(W))
  end

  return MapFromFunc(W, G, iso, isoinv)
end
