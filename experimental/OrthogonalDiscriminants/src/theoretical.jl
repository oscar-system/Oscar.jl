
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
      red = reduce(K(od), character_field(chi)[1])
      is_zero(red) && error("reduction mod $p is zero?")
      str = is_square(red) ? "O+" : "O-"
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


@doc raw"""
    od_from_p_subgroup(chi::GAPGroupClassFunction, p::Int[, pi::GAPGroupClassFunction])

Let `chi` be an irreducible (ordinary or Brauer) character of the group $G$,
and `p` be an odd prime integer.

Return `(flag, val)` where `flag` is `true` if and only if
enough information can be computed to prove that the restriction of `chi`
to a Sylow `p`-subgroup $P$ is orthogonally stable.
In this case,
`val` is a string that describes the orthogonal discriminant of `chi`.

If `flag` is `false` then `val` is equal to `""`.

If a character is given as the optional argument `pi` then it is assumed
that `pi` is the permutation character of $G$ induced from $P$.
Otherwise the function tries to compute the possible permutation characters.

# Examples
```jldoctest
julia> t = character_table("A5");

julia> Oscar.OrthogonalDiscriminants.od_from_p_subgroup(t[4], 5)
(true, "5")

julia> Oscar.OrthogonalDiscriminants.od_from_p_subgroup(mod(t, 3)[4], 5)
(true, "O-")

julia> Oscar.OrthogonalDiscriminants.od_from_p_subgroup(mod(t, 2)[4], 5)
(true, "O-")
```
"""
function od_from_p_subgroup(chi::GAPGroupClassFunction, p::Int)
  cand = possible_permutation_characters_from_sylow_subgroup(parent(chi), p)
  if cand != nothing
    res = collect(Set([od_from_p_subgroup(chi, p, pi) for pi in cand]))
    length(res) == 1 && return res[1]
  end
  return false, ""
end

function od_from_p_subgroup(chi::GAPGroupClassFunction, p::Int,
                            pi::GAPGroupClassFunction)
  tbl = parent(chi)
  Porder = div(order(tbl), degree(ZZRingElem, pi))
  flag, _, l = is_prime_power_with_data(Porder)
  @req (flag && l == p) "pi is not induced from a p-subgroup"
  p == 2 && return (false, "")

  l = characteristic(chi)
  l == p && return (false, "")

  # If `chi` is a Brauer character and `pi` is an ordinary character
  # then restrict `pi` to the `l`-regular classes.
  if characteristic(pi) == 0 && characteristic(chi) != 0
    pi = restrict(pi, tbl)
  end
  if scalar_product(chi, pi) != 0
    # The restriction of `chi` to P is not orthogonally stable.
    return false, ""
  end

  fus = filter(i -> pi[i] != 0, 1:number_of_conjugacy_classes(tbl))
  res = Oscar.GAPWrap.ELMS_LIST(GapObj(chi), GapObj(fus))

  if l != 0
    # possible shortcut:
    # `chi` itself can have a larger character field
    # than its restriction to P.
    # If this field contains a quadratic extension of the field of 'rest'
    # then we conclude that its OD is "O+",
    # no matter what the OD of the restriction to P is.
    Fchi = order_field_of_definition(chi)
    Fres = Oscar.GAPWrap.SizeOfFieldOfDefinition(res, l)
# We cannot form a class function object because we have no character table
# of the subgroup P (and we are happy that we do not need this table).
# Thus we cannot call `order_field_of_definition`.
# However, we know that the restriction to the non-vanishing classes of `pi`
# is closed under Galois conjugation,
# thus we would like to call `GAPWrap.SizeOfFieldOfDefinition` with this
# information.
#TODO: Let GAP support this.
    _, expchi, _ = is_prime_power_with_data(Fchi)
    _, expres, _ = is_prime_power_with_data(Fres)
    if mod(expchi, 2*expres) == 0
      return true, "O+"
    end
  end

  # We have to work.
  # First compute the ordinary OD of the restriction to P.
  # Compute the degree of the field `K` generated by the restriction,
  # and the degree of the enveloping cyclotomic field of `K`
  # or of the `p`-th cyclotomic field if `K` is the rational field.
  K = GAP.Globals.Field(GAP.Globals.Rationals, res)::GapObj
  pf = Oscar.GAPWrap.Conductor(K)
  if pf == 1
    pf = p
  end
  a = div(euler_phi(pf), GAP.Globals.Dimension(K))
  deg = degree(ZZRingElem, chi)
  @assert mod(deg, a) == 0
  KK, z = abelian_closure(QQ)
  if is_even(div(deg, a))
    # The value is a square in 'K'.
    od = KK(1)
  else
    # Compute \delta_K, as a `QQAbFieldElem`.
    # We need the norm of `z(p) + z(p)^-1 -2` w.r.t. the field extension
    # given by the real subfield of `cyclotomic_field(pf)` over `K`.
    # We compute this in the extension of `cyclotomic_field(p)` over
    # its intersection with `K`.
    # Note that this intersection is uniquely determined by `a`.
    # We compute the Galois group of the field extension.
    z = z(p)
    zinv = inv(z)
    od = z + zinv - 2
    galgen = GAP.Globals.PrimitiveRootMod(p)::Int
    galgen = powermod(galgen, Int(div(p-1, a)), p)
    e = galgen
    for k in 2:(a//2)
      od = od * (z^e + zinv^e - 2)
      e = mod(e * galgen, p)
    end
  end

  if l == 0
    # Turn the value into a string (using Atlas notation).
    return true, atlas_description(od)
  end

  # Reduce the ord. OD mod `l`.
  if conductor(od) == 1
    intod = ZZ(od)
    if is_odd(l)
      F = GF(l)
      return is_square(F(intod)) ? (true, "O+") : (true, "O-")
    elseif mod(intod, 8) == 1
      return (true, "O+")
    elseif mod(intod, 8) == 5
      return (true, "O-")
    else
      error("wrong congruence for OD of rational character")
    end
  end
  return (false, "")
end
