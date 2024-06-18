
# methods working with representations

# We assume that `G` is abs. irreducible and
# respects a nondegenerate quadratic form.
# In even characteristic (0 or 2), we assume that `G` is written over its
# character field.
# In odd characteristic, we assume that `G` is written over an odd degree
# extension of its character field.
#
# In characteristic zero, we have to interpret the result computed in
# `base_ring(G)` as an element of the character field.
# Thus we have to specify `embG` and `embF`, where
# `embG` is the embedding of `base_ring(G)` into the abelian closure of `QQ`
# and `embF` is the embedding of the character field of `G` into the abelian
# closure of `QQ`.
function od_from_generators(G::MatrixGroup, embG = nothing, embF = nothing)
  FG = base_ring(G)
  p = characteristic(FG)
  if p != 2
    # Compute the determinant of the invariant form as the determinant
    # of a skew-symmetric endomorphism w.r.t. the invariant form,
    # of full rank.
    # (In particular, do not compute the invariant form.)
    # The idea is that for each matrix `M` in `G`,
    # `Y = M - inv(M)` is such an endomorphism, i.e., satisfies
    # `B Y^{tr} inv(B) == -Y`, where `B` is the matrix of the invariant
    # bilinear form for `G`, i.e., `M B M^{tr} == B` holds for all `M` in `G`.
    # Thus we try to find a linear combination of such matrices `M - inv(M)`
    # that is invertible.
    # (In practice, very often three summands suffice.)
    M = matrix(rand_pseudo(G))
#T initialization is expensive if the matrices are large
    E = M - inv(M)
    for i in 1:100
      d = det(E)
      if is_zero(d)
        # Add/subtract another random element
        M = matrix(rand_pseudo(G))
        c = rand([1, -1])
        E = E + c * (M - inv(M))
      else
        # We have found an invertible skew-symmetric endomorphism.
        if p == 0
          # Switch from determinant to discriminant.
          if mod(degree(G), 4) == 2
            d = - d
          end
          
          # Coerce `d` into the character field of `G`.
          od = preimage(embF, embG(d))

          # Reduce the representative `od` mod obvious squares
          # in the character field.
          od = reduce_mod_squares(od)

          # Embed this value into the alg. closure.
          od = embF(od)

          # Turn the value into a string (using Atlas notation).
          return true, atlas_description(od)
        else
          # Decide if `d` is a square in the char. field.
          # By our assumptions, we can simply test whether `d` is a square
          # (in its parent, i.e., `base_ring(G)`).
          str = is_square(d) ? "O+" : "O-"
        end
      end
    end
  end

  # We known no better method than the MeatAxe approach.
  sgn = orthogonal_sign(G)
  if sgn == 1
    return true, "O+"
  elseif sgn == -1
    return true, "O-"
  end

  # This should not happen.
  error("assumptions not satisfied")
end


@doc raw"""
    od_from_atlas_group(chi::GAPGroupClassFunction)

Return `(flag, val)` where `flag` is `true` if the function can compute
the orthogonal discriminant of `chi`.
In this case, `val` is a string that describes the
orthogonal discriminant of `chi`.

Otherwise `flag` is `false` and `val` is `""`.

The former case happens if `chi` is absolutely irreducible,
has even degree and indicator `+`,
and the Atlas of Group Representations contains a representation
affording `chi`, over the character field of `chi` or (in odd characteristic)
a suitable extension field.

If the characteristic of `chi` is $2$ then we know no better method than
calling [`orthogonal_sign`](@ref).
Otherwise we try to find an invertible skew-symmetric endomorphism w.r.t.
the invariant bilinear form of the representation, with random methods.

```jldoctest
julia> t = character_table("A6");

julia> Oscar.OrthogonalDiscriminants.od_from_atlas_group(t[7])
(true, "-1")

julia> Oscar.OrthogonalDiscriminants.od_from_atlas_group(mod(t, 3)[4])
(true, "O-")

julia> t = character_table("L2(27)");

julia> Oscar.OrthogonalDiscriminants.od_from_atlas_group(t[4])
(false, "")
```
"""
function od_from_atlas_group(chi::GAPGroupClassFunction)
  tbl = parent(chi)
  chipos = findfirst(isequal(chi), tbl)
  chipos == nothing && return false, ""
  p = characteristic(chi)
  if p == 0
    id = identifier(tbl)
  else
    id = identifier(ordinary_table(tbl))
  end

  # Collect candidates that may afford `chi`.
  # Prefer an entry that knows to belong to `chi` or a member from its
  # orbit under group automorphisms,
  # but admit also other entries that do not know their characters,
  # provided that characteristic and degree determine the representation
  # up to group automorphisms.
#TODO: Improve by ruling out candidates via their character fields.
  d = degree(Int, chi)
  infos1 = []
  infos2 = []
  orbs = orbits_group_automorphisms(tbl)
  orb = orbs[findfirst(x -> chipos in x, orbs)]
  unique_orbit = all(x -> chipos in x || degree(Int, tbl[x[1]]) != d, orbs)
  for info in all_atlas_group_infos(id, characteristic => p, dim => d)
    if haskey(info, :constituents)
      # We trust in the precomputed data.
      if length(info[:constituents]) == 1 && info[:constituents][1] in orb
        push!(infos1, info)
      end
    elseif unique_orbit
      # Any abs. irred. character of degree `d` will be good.
      push!(infos2, info)
    end
  end

  # Run over the candidates (usually just one).
  for info in vcat(infos1, infos2)
    G = atlas_group(info)
    G == nothing && continue
    F = base_ring(G)
    K, embK = character_field(chi)
    if degree(F) != degree(K)
      if p == 2 || p == 0
        # For even `p`, we need representations over the character field.
        # Thus try to write the matrices over the character field.
        V = free_module(F, degree(G))
        M = gmodule(G, [hom(V, V, x) for x in [matrix(x) for x in gens(G)]])
        Mmin = Oscar.GModuleFromGap.gmodule_minimal_field(M)
#TODO: For p == 0, we need also the (!) embedding of `base_ring(Mmin)`
#      into `F == base_ring(M)`, because we know the embedding of `F`
#      into a cyclotomic field, and we will need this embedding in the end.
        F = base_ring(Mmin)
        (degree(F) != degree(K)) && continue
        G = matrix_group([matrix(x) for x in Mmin.ac])
      elseif is_even(degree(F) // degree(K))
        # For odd p, it is sufficient that `F` is an odd degree extension
        # of the character field.
        V = free_module(F, degree(G))
        M = gmodule(G, [hom(V, V, x) for x in [matrix(x) for x in gens(G)]])
        Mmin = Oscar.GModuleFromGap.gmodule_minimal_field(M)
        F = base_ring(Mmin)
        is_even(degree(F) // degree(K)) && continue
        G = matrix_group([matrix(x) for x in Mmin.ac])
      end
    end

    if ! haskey(info, :constituents)
      # We do not yet know that `G` has character `chi`.
      # Check whether the representation is abs. irred.
      if p == 0
        # We can check abs. irred. if we can reduce the repres.
        # modulo a prime l that does not divide the group order.
#TODO: Better create the natural module for G,
#T     and ask it whether it is abs. irreducible.
        l = next_prime(max([x[1] for x in factor(order(tbl))]...))
        Gbar, _ = Oscar.isomorphic_group_over_finite_field(G, min_char = Int(l))
        GG = GapObj(Gbar)
      else
        GG = GapObj(G)
      end
      M = GAP.Globals.GModuleByMats(GAP.Globals.GeneratorsOfGroup(GG),
                                    GAP.Globals.FieldOfMatrixGroup(GG))
      GAP.Globals.MTX.IsAbsolutelyIrreducible(M) || continue
    end

    # Now we are sure that `G` affords `chi`.
    if p == 0
#TODO: Add the relevant metadata to the Atlas of Group Representations.
if !(degree(F) <= 2 || Hecke.is_cyclotomic_type(F)[1])
  println( "missing field embedding info for ", info)
  continue
end
      KK, _ = abelian_closure(QQ)
      embG = Oscar.AbelianClosure._embedding(F, KK, KK(gen(F)))
      return od_from_generators(G, embG, embK)
    else
      return od_from_generators(G)
    end
  end

  return false, ""
end
