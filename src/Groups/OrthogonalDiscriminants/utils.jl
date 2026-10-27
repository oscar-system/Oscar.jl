@doc raw"""
    order_omega_mod_N(d::IntegerUnion, q::IntegerUnion, N::IntegerUnion) -> Pair{Bool, Bool}

Return `(flag_plus, flag_minus)` where `flag_plus` and `flag_minus`
are `true` or `false`, depending on whether `N` divides the order
of the orthogonal groups $\Omega^+(d, q)$ and $\Omega^-(d, q)$.

# Examples
```jldoctest
julia> Oscar.OrthogonalDiscriminants.order_omega_mod_N(4, 2, 60)
(false, true)

julia> Oscar.OrthogonalDiscriminants.order_omega_mod_N(4, 5, 60)
(true, true)
```
"""
function order_omega_mod_N(d::IntegerUnion, q::IntegerUnion, N::IntegerUnion)
  @req is_even(d) "d must be even"
  m = div(d, 2)
  exp, N = remove(N, q)
  facts = collect(factor(q))
  p = facts[1][1]
  if mod(N, p) == 0
    exp = exp + 1
    _, N = remove(N, p)
  end
  if m*(m-1) < exp
    # A group of order `N` does not embed in any candidate.
    return (false, false)
  end

  q2 = ZZ(q)^2
  q2i = ZZ(1)
  for i in  1:(m-1)
    q2i = q2 * q2i
    if i == 1 && is_odd(q)
      g = gcd(N, div(q2i-1, 2))
    else
      g = gcd(N, q2i-1)
    end
    N = div(N, g)
    if N == 1
      # A group of order N may embed in both candidates.
      return (true, true)
    end
  end

  # embeds in + type?, embeds in - type?
  return (mod(q^m-1, N) == 0, mod(q^m+1, N) == 0)
end


@doc raw"""
    reduce_mod_squares(val::AbsSimpleNumFieldElem)

Return an element of `F = parent(val)` that is equal to `val`
modulo squares in `F`.

If `val` describes an integer then the result corresponds to the
squarefree part of this integer.
Otherwise the coefficients of the result have a squarefree g.c.d.

# Examples
```jldoctest
julia> F, z = cyclotomic_field(4);

julia> Oscar.OrthogonalDiscriminants.reduce_mod_squares(4*z^0)
1

julia> Oscar.OrthogonalDiscriminants.reduce_mod_squares(-8*z^0)
-2
```
"""
function reduce_mod_squares(val::AbsSimpleNumFieldElem)
  is_zero(val) && return val
  d = denominator(val)
  if ! isone(d)
    val = val * d^2
  end
  F = parent(val)
  if is_integer(val)
    intval = ZZ(val)
    sgn = sign(intval)
    good = [x[1] for x in factor(intval) if is_odd(x[2])]
    return F(prod(good, init = sgn))
  elseif is_square(val)
    return F(1)
  end
  # Just get rid of the square part of the gcd of the coefficients.
  c = map(numerator, coefficients(val))
  s = 1
  for (p, e) in factor(gcd(c))
    if iseven(e)
      s = s * p^e
    elseif e > 1
      s = s * p^(e-1)
    end
  end
  return val//s
end


@doc raw"""
    possible_permutation_characters_from_sylow_subgroup(tbl::Oscar.GAPGroupCharacterTable, p::Int)

Return either `nothing` or the vector of all those characters of `tbl` that
may be equal to the permutation character $1_P^G$,
where $G$ is the group of `tbl` and $P$ is a Sylow `p`-subgroup of $G$.

# Examples
```jldoctest
julia> tbl = character_table("A5");

julia> l = Oscar.OrthogonalDiscriminants.possible_permutation_characters_from_sylow_subgroup(tbl, 2);

julia> length(l)
1

julia> println(coordinates(Int, l[1]))
[1, 0, 0, 1, 2]
```
"""
function possible_permutation_characters_from_sylow_subgroup(tbl::Oscar.GAPGroupCharacterTable, p::Int)
  @req is_prime(p) "p must be a prime integer"
  l = characteristic(tbl)
  if l == p
    # The permutation character is zero on all
    # nonidentity `p`-regular classes.
    pi = fill(0, nrows(tbl))
    # index of the Sylow `p`-subgroup
    _, n = remove(order(tbl), p)
    pi[1] = n
    return [Oscar.class_function(tbl, GapObj(pi, true))]
  elseif l != 0
    # Compute the corresponding ordinary characters, and restrict them.
    pi = possible_permutation_characters_from_sylow_subgroup(
             ordinary_table(tbl), p)
    if pi != nothing
      pi = [restrict(x, tbl) for x in pi]
    end
    return pi
  end

  # index of the Sylow 'p'-subgroup
  q, n = remove(order(tbl), p)
  q = p^q

  # Perhaps the Sylow `p`-subgroup is cyclic.
  pos = findfirst(==(q), orders_class_representatives(tbl))
  pos != nothing && return [induced_cyclic(tbl, [pos])[1]]

  # Do we have the table of marks?
  tom = GAP.Globals.TableOfMarks(GapObj(tbl))
  if tom != GAP.Globals.fail
    pos = findfirst(==(q), Vector{Int}(GAP.Globals.OrdersTom(tom)))
    return [Oscar.class_function(tbl, GAP.Globals.PermCharsTom(GapObj(tbl), tom)[pos])]
  end

  # Does the table have a nontrivial 'p'-core
  # such that the table of the factor is a library table?
  # If yes and if we can compute the result for the factor
  # then inflate it to 'tbl'.
  npos = class_positions_of_pcore(tbl, p)
  if length(npos) != 1
    for d in known_class_fusions(tbl)
      if class_positions_of_kernel(d[2]) == npos
        s = character_table(d[1])
        if s != nothing && characteristic(s) == 0
          # Try to recurse.
          pi = possible_permutation_characters_from_sylow_subgroup(s, p)
          if pi != nothing
            return [restrict(x, tbl) for x in pi]
          end
        end
      end
    end
  end

  candlist = nothing
  for name in names_of_fusion_sources(tbl)
    s = character_table(name)
    if s != nothing && characteristic(s) == 0
      flag, fus = known_class_fusion(s, tbl)
      if flag && length(class_positions_of_kernel(fus)) == 1
        if mod(order(s), q) == 0
          # Run over the known character tables of subgroups
          # that contain Sylow subgroups of `tbl`.
          npos = class_positions_of_pcore(s, p)
          if sum(class_lengths(s)[npos]) == q
            # The Sylow `p`-subgroup of `tbl` is normal in `s`.
            pi = fill(0, nrows(s))
            index = div(order(s), q)
            for i in npos
              pi[i] = index
            end
            pi = GAP.Globals.ClassFunction(GapObj(s), GapObj(pi, true))
            return [Oscar.class_function(s, pi)^tbl]
          else
            # Try to recurse.
            pi = possible_permutation_characters_from_sylow_subgroup(s, p)
            if pi != nothing
              pi = collect(Set([x^tbl for x in pi]))
              if candlist == nothing
                candlist = pi
              else
                intersect!(candlist, pi)
              end
              length(candlist) == 1 && return candlist
            end
          end
        elseif mod(div(order(tbl), order(s)), p) == 0 &&
               is_prime_power_with_data(div(order(tbl), order(s)))[1]
          # Compute the perm. char. of the Sylow p-subgroup in the subgroup,
          # and then try to extend the character uniquely to G.
          cand = possible_permutation_characters_from_sylow_subgroup(s, p)
          if cand != nothing
            extcand = []
            for pi in cand
              pi = GAP.Globals.CompositionMaps(GapObj(pi),
                       GAP.Globals.InverseMap(GapObj(fus)))
              union!(extcand, [Oscar.class_function(tbl, x) for x in
                                GAP.Globals.PermChars(GapObj(tbl),
                                  GapObj(Dict(:torso => pi), true))])
            end
            if candlist == nothing
              candlist = extcand
            else
              intersect!(candlist, extcand)
            end
            length(candlist) == 1 && return candlist
          end
        end
      end
    end
  end

  isdefined(tbl, :group) && return [trivial_character(sylow_subgroup(group(tbl), p)[1])^tbl]

  # If we arrive here,
  # we can try to recurse to subgroups for which fusions must be composed,
  # or to compute candidates.
  # (Eventually a function that returns all candidates might be useful then.)
  return candlist
end


# helper function:
# array of identifiers of available character tables of (non-dashed)
# automorphic extensions of a simple character table.
function _identifiers_of_almost_simple_tables(simplename::String)
  gapname = GapObj(simplename)
  nams = GAP.Globals.CTblLib.DisplayAtlasMap_ComputePortions(
           gapname, 0).identifiers[1]
  return filter(x -> x == gapname || x[end] != '\'', Vector{String}(nams))
end

# helper function:
# Let `M` be a vector of vectors representing a *regular* matrix
# (this is not checked),
# and let `column_orbits` be the set of orbits of some group `G`
# of matrix automorphisms of `M` on the positions of columns.
# Return the set of orbits of the corresponding action of `G`
# on the rows of `M`.
function _row_orbits_from_column_orbits(M::Vector, column_orbits::Vector{Vector{Int}})
  # Compute the auxiliary matrix of orbit sums.
  n = length(M)
  MM = [[] for _ in 1:n]
  for omega in column_orbits
    for i in 1:n
      push!(MM[i], sum(M[i][omega]))
    end
  end

  # Compute the orbits on rows.
  row_orbits = Vector{Int}[]
  reps = []
  poss = []
  for i in 1:n
    pos = findfirst(isequal(MM[i]), reps)
    if pos == nothing
      push!(row_orbits, [i])
      push!(reps, MM[i])
      push!(poss, length(row_orbits))
    else
      push!(row_orbits[poss[pos]], i)
    end
  end

  return row_orbits
end

# helper function:
# Compute the orbits obtained by joining the orbits in `orbs1`
# with those in `orbs2`.
function _common_orbits(orbs1::Vector{Vector{Int}}, orbs2::Vector{Vector{Int}})
    local orbs, orb, pos;

    orbs = orbs1
    for orb in orbs2
      inter = Int[]
      l = length(orbs)
      for i in 1:l
        is_empty(intersect(orb, orbs[i])) && continue
        push!(inter, i)
      end
      orbs = union([union(orbs[inter]...)], orbs[setdiff(1:l, inter)])
    end
    return sort!(map(sort!, orbs))
end

# helper function:
# Compute the inverse of a class fusion `fus`,
# that is, the vector of length `n` that contains at position `i`
# the (perhaps empty) vector of all those indices `j` in `fus`
# such that `fus[j] == i` holds.
# We assume that `n` is not smaller than the maximum of `fus`.
function _inverse_fusion(fus::Vector{Int}, n::Int)
  inv = [Int[] for i in 1:n]
  for i in 1:length(fus)
    push!(inv[fus[i]], i)
  end
  return inv
end

@doc raw"""
    orbits_group_automorphisms(tbl::Oscar.GAPGroupCharacterTable)

Return an array of arrays, each entry listing an orbit under the action of
those group automorphisms which are described by class fusions
that are stored on `tbl`, such that the group of `tbl` is a normal subgroup
of the group given by the image of the fusion.

The result describes the action of the full automorphism group only if
the relevant tables of automorphic extensions and the fusions to them
are available.

We assume that `tbl` is almost simple.

# Examples
```jldoctest
julia> Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
         character_table("L2(8)"))
5-element Vector{Vector{Int64}}:
 [1]
 [2]
 [3, 4, 5]
 [6]
 [7, 8, 9]
```
"""
function orbits_group_automorphisms(tbl::Oscar.GAPGroupCharacterTable)
  # 'tbl' is an upwards extension of the table of a simple group.
  # We compute the almost simple tables for this simple group,
  # and then take only those that correspond to automorphisms of 'tbl'.
  p = characteristic(tbl)
  if p == 0
    ordtbl = tbl
  else
    ordtbl = ordinary_table(tbl)
  end
  name = identifier(ordtbl)
  pos = findfirst('.', name)
  if pos != nothing
    simpname = name[1:(pos-1)]
  else
    simpname = name
  end
  n = nrows(tbl)
  orbs = [[i] for i in 1:n]
  ordauttbls = [character_table(x) for x in _identifiers_of_almost_simple_tables(simpname)]
  for ordauttbl in ordauttbls
    (order(ordauttbl) == order(tbl) || mod(order(ordauttbl), order(tbl)) != 0) && continue
    if p == 0
      auttbl = ordauttbl
    else
      auttbl = mod(ordauttbl, p)
      auttbl == nothing && continue
    end

    fl, fus = known_class_fusion(tbl, auttbl)
    if ! fl
      # If 'auttbl' does not count for 'ordtbl' then this is o.k.
      if length(possible_class_fusions(ordtbl, ordauttbl)) != 0
        # This criterion is sufficient for our almost simple tables.
        error("class fusion from $tbl to $auttbl missing?")
      end
      continue
    end

    colorbs = _inverse_fusion(fus, max(fus...))
    colorbs = filter(x -> length(x) != 0, colorbs)
    roworbs = _row_orbits_from_column_orbits(map(collect, collect(tbl)), colorbs)
    orbs = _common_orbits(orbs, roworbs)
  end

  return orbs
end


@doc raw"""
    is_orthogonally_stable(chi::GAPGroupClassFunction; check::Bool = true)

Return `nothing` if the indicator of some irreducible constituent of `chi`
is not known; this can happen only if `chi` has characteristic 2.

Otherwise return `true` if `chi` is orthogonally stable,
and `false` otherwise.

A character is called orthogonally stable if
- `chi` is orthogonal, that is, `chi` is real,
  and all its absolutely irreducible constituents of indicator `-`
  have even multiplicity and
- all its absolutely irreducible constituents of indicator `+`
  have even degree.

If we know that `chi` is orthogonal then we can set `check` to `false`;
in this case, some `nothing` results can be avoided.

# Examples
```jldoctest
julia> t = character_table("A6");

julia> println(map(is_orthogonally_stable, t))
Bool[0, 0, 0, 1, 1, 0, 1]

julia> println(map(is_orthogonally_stable, mod(t, 3)))
Bool[0, 0, 0, 1, 0]
```
"""
function is_orthogonally_stable(chi::GAPGroupClassFunction; check::Bool = true)
  chi != conj(chi) && return false

  # Now we know that `chi` is real.
  tbl = parent(chi)
  dec = coordinates(ZZRingElem, chi)
  ind = [indicator(chi) for chi in tbl]
  deg = [degree(ZZRingElem, chi) for chi in tbl]

  for i in 1:length(ind)
    if ind[i] == -1
      if is_odd(dec[i])
        # `chi` is not orthogonal.
        return false
      end
    elseif ind[i] == 1
      if dec[i] != 0 && is_odd(deg[i])
        # `chi` is not orthogonally stable.
        return false
      end
    elseif ind[i] != 0
      # The indicator is not known.
      # The value can be 1 or -1, so test both conditions.
      if is_odd(dec[i])
        check && return nothing
        # If we know that `chi` is the reduction of an orthogonal
        # character then each indicator `-` constituent occurs with
        # even multiplicity, thus this constituent has indicator `+`.
        # Check whether it has even degree.
        is_odd(deg[i]) && return false
      elseif dec[i] != 0 && is_odd(deg[i])
        # There is a constituent of odd degree, which may or may not
        # have indicator `+`.
        return nothing
      end
    end
  end
  return true
end


#   _orthogonal_discriminant_indicator0(chi::GAPGroupClassFunction)
#
# Return an string that describes the orthogonal discriminant of
# `chi + conj(chi)`, where `chi` must be absolutely irreducible and non-real.
#
function _orthogonal_discriminant_indicator0(chi::GAPGroupClassFunction)
  @req chi != conj(chi) "chi must be non-real"

  tbl = parent(chi)
  p = characteristic(tbl)

  # For Brauer characters, it may happen that the character fields of
  # `chi` and `chi + conj(chi)` are equal,
  # and then Theorem 10 yields "O+".
  # Otherwise, apply Prop. 2.9.
  if p > 0
    F, _ = character_field(chi)
    K, _ = character_field(chi + conj(chi))
    if degree(F) == degree(K)
      # happens for example for L3(2), degree 3, lives over GF(2)
      return "O+";
    else
      # We have "O+" iff the degree of `chi` is even,
      # also for odd characteristic.
      return is_even(degree(ZZRingElem, chi)) ? "O+" : "O-"
    end
  end

  # If the degree of `chi` is even then the discr. is 1.
  if is_even(degree(ZZRingElem, chi))
    return "1"
  end

  # Let `F` be the character field (not real), and `K` the real subfield.
  F, embF = character_field(chi)

  # Find delta in `K` such that `F` = `K(sqrt(delta))`.
  # Then the discr. is `delta`.
  if degree(F) == 2
    # We find a negative integer (nice).
    return string(Oscar.AbelianClosure.quadratic_irrationality_info(embF(gen(F)))[3])
  end

  # If `F` contains a quadratic field `Q` that is not real
  # then choose the square root from `Q`, again get an integral `delta`.
  for (Q, Qemb) in subfields(F, degree = 2)
    if ! is_real(complex_embeddings(Q, conjugates = false)[1])
      return string(AbelianClosure.quadratic_irrationality_info(embF(embQ(gen(Q))))[3])
    end
  end

  # Reduce from `F` to a subfield of 2-power degree,
  # we can choose `delta` in an index 2 subfield of this subfield.
  e, _ = remove(degree(F), 2)
  for (S, embS) in subfields(F, degree = 2^e)
    # Construct a nonreal element in `F` with real square.
    # (Usually this is not a nice element.)
    x = embF(embS(gen(S)))
    return atlas_description((x - conj(x))^2)
  end
end


@doc raw"""
    show_with_ODs(tbl::Oscar.GAPGroupCharacterTable)

Show `tbl` with 2nd indicators, known ODs, and degrees of character fields.
(See [`Base.show(io::IO, ::MIME"text/plain", tbl::Oscar.GAPGroupCharacterTable)`](@ref)
for ways to modify what is shown.)

# Examples
```jldoctest
julia> t = character_table("A5");

julia> Oscar.OrthogonalDiscriminants.show_with_ODs(t)
A5

          2  2  2  .  .  .
          3  1  .  1  .  .
          5  1  .  .  1  1
                          
            1a 2a 3a 5a 5b
         2P 1a 1a 3a 5b 5a
         3P 1a 2a 1a 5b 5a
         5P 1a 2a 3a 1a 1a
    d OD  2               
X_1 1     +  1  1  1  1  1
X_2 2     +  3 -1  .  A A*
X_3 2     +  3 -1  . A*  A
X_4 1  5  +  4  .  1 -1 -1
X_5 1     +  5  1 -1  .  .

A = z_5^3 + z_5^2 + 1
A* = -z_5^3 - z_5^2
```
"""
function show_with_ODs(tbl::Oscar.GAPGroupCharacterTable, io::IO = stdout)
   iob = IOBuffer()
   show(IOContext(iob, :indicator => true,
                       :OD => true,
                       :character_field => true,
                       :with_legend => true), MIME("text/plain"), tbl)
   print(io, String(take!(iob)))
   return
end


@doc raw"""
    show_OD_info(tbl::Oscar.GAPGroupCharacterTable)
    show_OD_info(name::String)

Show an overview of known information about the ordinary and modular
orthogonal discriminants for `tbl` or for the character table with
identifier `name`.

# Examples
```jldoctest
julia> show_OD_info("A5")
A5:  2^2*3*5
------------

i|chi|K|disc| 2| 3|       5
-+---+-+----+--+--+--------
4| 4a|Q|   5|4a|4a|(def. 1)
 |   | |    |O-|O-|        
```
"""
function show_OD_info(tbl::Oscar.GAPGroupCharacterTable, io::IO = stdout)
  id = identifier(tbl)
  haskey(OD_simple_names, id) || return
  simp = OD_simple_names[id]
  data = OD_data[simp]
  haskey(data, id) || return
  data = data[id]

  positions = [entry[2] for entry in data["0"]]
  ord = order(tbl)
  primes = prime_divisors(ord)
  factord = sort([pair for pair in factor(ord)])
  primes = [pair[1] for pair in factord]
  modtbls = [mod(tbl, p) for p in primes]
  header = "$(identifier(tbl)):  $(_string_factored(factord))"
  if is_unicode_allowed()
    header = [header, repeat("─", length(header)), ""]
  else
    header = [header, repeat("-", length(header)), ""]
  end

  labels_col = ["chi", "K", "disc"]
  append!(labels_col, [string(p) for p in primes])

  result = []
  labels_row = String[]

  for i in 1:length(positions)
    chipos = positions[i]
    push!(labels_row, string(chipos))
    push!(labels_row, "")
    chi = tbl[positions[i]]
    resulti = [[_character_name(tbl, chipos),
                _string_character_field(chi),
                data["0"][i][4]],
               ["", "", ""]]
    if mod(degree(ZZRingElem, chi), 4) == 2 && resulti[1][3] == "?"
      resulti[1][3] = "-?"
    end
    for j in 1:length(primes)
      p = primes[j]
      if modtbls[j] == nothing
        # We do not know the `p`-Brauer table.
        push!(resulti[1], "?")
        push!(resulti[2], "")
      else
        red = restrict(chi, modtbls[j])
        if is_orthogonally_stable(red) != false
          # (The result may be `nothing`, meaning that
          # the indicator for at least one constituent is not known.)
          dec = coordinates(ZZRingElem, red)
          res1 = String[]
          res2 = String[]
          for k in filter(x -> is_odd(dec[x]), 1:length(dec))
            ppos = filter(r -> r[2] == k, data[string(p)])
            if length(ppos) > 0
              ppos = ppos[1]
            end
            if length(ppos) != 0 && !("(indicator unknown)" in ppos[end])
              # must be indicator `+`
              push!(res1, _character_name(modtbls[j], k))
              push!(res2, ppos[4])
            else
              cc = findfirst(is_equal(conj(modtbls[j][k])), modtbls[j])
              if cc > k
                # indicator `o`
                push!(res1, _character_name(modtbls[j], k) *
                            filter(!isdigit, _character_name(modtbls[j], cc)))
                push!(res2, _orthogonal_discriminant_indicator0(modtbls[j][k]))
              elseif cc == k
                # indicator is *unknown*;
                # note that indicator `-` (disc is "O+") cannot occur here
                push!(res1, _character_name(modtbls[j], k) * "(?)")
                if length(ppos) == 0
                  push!(res2, "(?)")
                else
                  push!(res2, ppos[4] * "(?)")
                end
              end
            end
          end
          push!(resulti[1], join(res1, "+"))
          push!(resulti[2], join(res2, ", "))
        else
          bl = block_distribution(tbl, p)
          if bl[:defect][bl[:block][chipos]] == 1
            push!(resulti[1], "(def. 1)")
          else
            push!(resulti[1], "")
          end
          push!(resulti[2], "")
        end
      end
    end
    append!(result, resulti)
  end

  ioc = IOContext(io,
         :header => header,
         :corner => ["i"],
         :labels_row => labels_row,
         :labels_col => labels_col,
         :separators_row => 0:2:(length(result)-2),
         :separators_col => 0:(length(result[1])-1),
        )

  labeled_matrix_formatted(ioc, permutedims(reduce(hcat, result)))
  println(io, "")
end

show_OD_info(name::String, io::IO = stdout) = show_OD_info(character_table(name), io)

function _string_factored(pairs::Vector{Pair{ZZRingElem, Int}})
  l = map(pair -> pair[2] == 1 ? "$(pair[1])" : "$(pair[1])^$(pair[2])", pairs)
  return join(l, "*")
end

function _character_name(tbl::Oscar.GAPGroupCharacterTable, i::Int)
  deg = degree(ZZRingElem, tbl[i])
  n = count(j -> degree(ZZRingElem, tbl[j]) == deg, 1:i)
  alp = ["$x" for x in "abcdefghijklmnopqrstuvwxyz"]
  while n > length(alp)
    append!(alp, ["($x')" for x in alp])
  end

  return string(deg) * alp[n]
end

function _string_character_field(chi::Oscar.GAPGroupClassFunction)
  F, emb = character_field(chi)
  degree(F) == 1 && return "Q"
  names = map(atlas_description, [emb(x) for x in gens(F)])
  return "Q(" * join(names, ", " ) * ")"
end
