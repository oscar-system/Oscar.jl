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
    reduce_mod_squares(val::nf_elem)

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
function reduce_mod_squares(val::nf_elem)
  is_zero(val) && return val
  d = denominator(val)
  if ! isone(d)
    val = val * d^2
  end
  F = parent(val)
  if is_integer(val)
    intval = ZZ(val)
    sgn = sign(intval)
    good = [x[1] for x in collect(factor(intval)) if is_odd(x[2])]
    return F(prod(good, init = sgn))
  elseif is_square(val)
    return F(1)
  end
  # Just get rid of the square part of the gcd of the coefficients.
  c = map(numerator, coefficients(val))
  s = 1
  for (p, e) in collect(factor(gcd(c)))
    if iseven(e)
      s = s * p^e
    elseif e > 1
      s = s * p^(e-1)
    end
  end
  return val//s
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
    show_with_ODs(tbl::Oscar.GAPGroupCharacterTable)

Show `tbl` with 2nd indicators, known ODs, and degrees of character fields.
(See [`Base.show(io::IO, ::MIME"text/plain", tbl::GAPGroupCharacterTable)`](@ref)
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
