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
  if is_integer(val)
    intval = ZZ(val)
    sgn = sign(intval)
    good = [x[1] for x in collect(factor(intval)) if is_odd(x[2])]
    F = parent(val)
    return F(prod(good, init = sgn))
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
            if length(ppos) != 0 && !contains(ppos[end], "(indicator unknown)")
              # must be indicator `+`
              push!(res1, _character_name(modtbls[j], k))
              push!(res2, ppos[4])
            else
              cc = findfirst(is_equal(conj(modtbls[j][k])), modtbls[j])
              if cc > k
                # indicator `o`
                push!(res1, _character_name(modtbls[j], k) *
                            filter(!isdigit, _character_name(modtbls[j], cc)))
                push!(res2, orthogonal_discriminant_indicator0(modtbls[j][k]))
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

  labelled_matrix_formatted(ioc, permutedims(reduce(hcat, result)))
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
