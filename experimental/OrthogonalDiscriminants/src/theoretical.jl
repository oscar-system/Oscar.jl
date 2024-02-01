
# character-theoretical methods


@doc raw"""
    od_from_order(chi::GAPGroupClassFunction)

Return `(flag, val)` where `flag` is `true` if the order of the group
of `chi` divides only one of the orders of the two orthogonal groups
[`omega_group`](@ref)`(+/-1, d, q)`, where `d` is the degree of `chi`
and `q` is the order of the field of definition of `chi`.

In this case, `val` is `"O+"` or `"O-"`.

```jldoctest
julia> t = character_table("L3(2)");

julia> Oscar.OrthogonalDiscriminants.od_from_order(mod(t, 3)[4])
(true, "O-")

julia> Oscar.OrthogonalDiscriminants.od_from_order(mod(t, 2)[4])
(false, "")
```
"""
function od_from_order(chi::GAPGroupClassFunction)
  characteristic(chi) == 0 && return (false, "")
  d = numerator(degree(chi))
  q = order_field_of_definition(chi)
  tbl = ordinary_table(parent(chi))
  ord = order(ZZRingElem, tbl)

  # Compute the order of the subgroup that shall embed into
  # the perfect group `omega_group(epsilon, d, q)`.
  n = sum(class_lengths(tbl)[class_positions_of_solvable_residuum(tbl)])
  flag1, flag2 = order_omega_mod_N(d, q, n)
  if flag1
    return flag2 ? (false, "") : (true, "O+")
  else
    return flag2 ? (true, "O-") : (false, "")
  end
end


@doc raw"""
    od_from_eigenvalues(chi::GAPGroupClassFunction)

Return `(flag, val)` where `flag` is `true` if there is a conjugacy class
on which representing matrices for `chi` have no eigenvalue $\pm 1$.
In this case, if `chi` is orthogonally stable (this is not checked here)
then `val` is a string that describes the orthogonal discriminant of `chi`.

If `flag` is `false` then `val` is equal to `""`.

This criterion works only if the characteristic of `chi` is not $2$,
`(false, "")` is returned if the characteristic is $2$.

# Examples
```jldoctest
julia> t = character_table("A5");

julia> Oscar.OrthogonalDiscriminants.od_from_eigenvalues(t[4])
(true, "5")

julia> Oscar.OrthogonalDiscriminants.od_from_eigenvalues(mod(t, 3)[4])
(true, "O-")

julia> Oscar.OrthogonalDiscriminants.od_from_eigenvalues(mod(t, 2)[4])
(false, "")
```
"""
function od_from_eigenvalues(chi::GAPGroupClassFunction)
  p = characteristic(chi)
  p == 2 && return (false, "")

  tbl = parent(chi)
  ord = orders_class_representatives(tbl)
  for i in 2:length(chi)
    n = ord[i]
    ev = multiplicities_eigenvalues(chi, i)
    if ev[end] != 0 || (iseven(n) && ev[divexact(n, 2)] != 0)
      continue
    end

    F, z = cyclotomic_field(n)
    od = prod(x -> x[1]^x[2], [(z^i-z^-i, ev[i]) for i in 1:n])
    if mod(degree(chi), 4) == 2
      od = -od
    end

    K, _ = abelian_closure(QQ)
    if p == 0
      # Coerce `od` into the character field of `chi`.
      F, emb = character_field(chi)
      od = preimage(emb, K(od))

      # Reduce the representative `od` mod obvious squares
      # in the character field.
      od = reduce_mod_squares(od)

      # Embed this value into the alg. closure.
      od = emb(od)

      # Turn the value into a string (using Atlas notation).
      str = atlas_description(od)
    else
      # Decide if the reduction mod `p` is a square in the char. field.
      str = is_square(reduce(K(od), character_field(chi)[1])) ? "O+" : "O-"
    end

    return true, str
  end

  return false, ""
end

@doc raw"""
    od_for_specht_module(chi::GAPGroupClassFunction)

Return `(flag, val)` where `flag` is `true` if `chi` is an ordinary
irreducible character of a symmetric group or of an alternating group
such that `chi` extends to the corresponding symmetric group.
In this case, if `chi` is orthogonally stable (this is not checked here)
then `val` is a string that describes the orthogonal discriminant of `chi`;
the discriminant is computed using the Jantzen-Schaper formula,
via [`gram_determinant_specht_module`](@ref).

`(false, "")` is returned in all cases where this criterion is not applicable.

# Examples
```jldoctest
julia> t = character_table("A5");

julia> Oscar.OrthogonalDiscriminants.od_for_specht_module(t[4])
(true, "5")

julia> Oscar.OrthogonalDiscriminants.od_for_specht_module(mod(t, 3)[4])
(false, "")
```
"""
function od_for_specht_module(chi::GAPGroupClassFunction)
  characteristic(chi) == 0 || return (false, "")

  # Find out to which alternating or symmetric group `chi` belongs.
  tbl = parent(chi)
  name = identifier(tbl)
  startswith(name, "A") || return (false, "")
  pos = findfirst('.', name)
  if pos === nothing
    n = parse(Int, name[2:end])
  else
    n = parse(Int, name[2:(pos-1)])
  end
  n === nothing && return (false, "")
  name == "A$n" || name == "A$n.2" || name == "A6.2_1" || return (false, "")

  chipos = findfirst(isequal(chi), tbl)
  chipos === nothing && return (false, "")
  para = character_parameters(tbl)[chipos]
  isa(para, Vector{Int}) || return (false, "")

  # Now we know that `chi` belongs to Sym(n) or extends to Sym(n)
  gramdet = gram_determinant_specht_module(partition(para))
  res = ZZRingElem(1)
  for pair in gramdet
    if is_odd(pair[2])
      res = res * pair[1]
    end
  end
  if mod(degree(ZZRingElem, chi), 4) == 2
    res = - res
  end

  return true, string(res)
end
