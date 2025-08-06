
# functions for accessing the precomputed data

# the precomputed data for Atlas groups
const OD_data = JSON.parsefile(joinpath(@__DIR__, "../data/odresults.json"))

const OD_simple_names = Dict{String, String}()

@doc raw"""
    OD_split_string(str::String, sep::String)

Return a vector of strings obtained by splitting `str` at characters in `sep`,
but only at those positions where the brackets `(` and `)` are balanced.
"""
function OD_split_string(str::String, sep::String)
  open = 0
  res = String[]
  start = 0
  for i in 1:length(str)
    if str[i] == '('
      open = open + 1
    elseif str[i] == ')'
      open = open - 1
    elseif str[i] in sep && open == 0
      push!(res, str[(start+1):(i-1)])
      start = i
    end
  end
  push!(res, str[(start+1):length(str)])

  return res
end


for simpnam in OD_data["names"]
  for x in OD_data[simpnam]["names"]
    OD_simple_names[x] = simpnam
    for p in keys(OD_data[simpnam][x])
      for v in OD_data[simpnam][x][p]
        v[end] = OD_split_string(v[end], ",")
      end
    end
  end
end


@doc raw"""
    orthogonal_discriminant(chi::Oscar.GAPGroupClassFunction)

Return a string that describes the orthogonal discriminant of `chi`:
- `"?"` if the value is unknown,
- one of `"O+"`, `"O-"` in positive characteristic,
- something that can be evaluated with [`atlas_irrationality`](@ref)
  in characteristic `0`, and
- `""` if `chi is not irreducible, not orthogonal, or has odd degree.

# Examples
```jldoctest
julia> t = character_table("A6");

julia> orthogonal_discriminant(t[4])
"1"

julia> t2 = t % 2;

julia> orthogonal_discriminant(t2[4])
"O+"
```
"""
function orthogonal_discriminant(chi::Oscar.GAPGroupClassFunction)
  tbl = parent(chi)
  pos = findfirst(x -> x === chi, tbl)
  pos === nothing && return ""
  p = characteristic(tbl)
  if p == 0
    id = identifier(tbl)
  else
    id = identifier(ordinary_table(tbl))
  end
  if haskey(OD_simple_names, id)
    simp = OD_simple_names[id]
    data = OD_data[simp]
    data = data[id]
    p = string(p)
    if haskey(data, p)
      data = data[p]
      for l in data
        l[2] == pos && return l[4]
      end
      return ""
    else
      (indicator(chi) == 1 && is_even(numerator(degree(chi)))) || return ""
      return "?"
    end
  end

  # `tbl` is outside the scope of the database.
  (indicator(chi) == 1 && is_even(numerator(degree(chi)))) || return ""
  return "?"
end


@doc raw"""
    orthogonal_discriminants(tbl::Oscar.GAPGroupCharacterTable)

Return a vector of strings that describe the orthogonal discriminants of
the orthogonal irreducible characters of `tbl` of even degree.

The length of this vector is the number of irreducible characters of `tbl`,
the $i$-th entry is an empty string if the $i$-th character is not orthogonal
or has odd degree,
and the $i$-th entry is equal to `"?"` if the orthogonal discriminant is
unknown.

# Examples
```jldoctest
julia> t = character_table("A6");

julia> println(orthogonal_discriminants(t))
["", "", "", "1", "1", "", "-1"]

julia> println(orthogonal_discriminants(t % 3))
["", "", "", "O-", ""]
```
"""
function orthogonal_discriminants(tbl::Oscar.GAPGroupCharacterTable)
  p = characteristic(tbl)
  if p == 0
    ordtbl = tbl
  else
    ordtbl = ordinary_table(tbl)
  end
  id = identifier(ordtbl)
  res = fill("", number_of_conjugacy_classes(tbl))
  if haskey(OD_simple_names, id)
    simp = OD_simple_names[id]
    data = OD_data[simp]
    if haskey(data, id)
      data = data[id]
      p = string(p)
      if haskey(data, p)
        data = data[p]
        for l in data
          res[l[2]] = l[4]
        end
        return res
      else
#TODO: compute the reductions mod p (not dividing the group order)
      end
    end
  else
    # Perhaps a factor of `tbl` belongs to the database.
    facttbl = nothing
    for r in known_class_fusions(ordtbl)
      if 1 < length(class_positions_of_kernel(r[2]))
        id = r[1]
        if haskey(OD_simple_names, id)
          simp = OD_simple_names[id]
          data = OD_data[simp]
          if haskey(data, id)
            data = data[id]
            facttbl = character_table(id)
            if haskey(data, string(p))
              if p != 0
                facttbl = mod(facttbl, p)
              end
              facttbl != nothing && break
            end
          end
        end
      end
    end
    if facttbl != nothing
      # Set the known values for the factor.
      mp = [findfirst(is_equal(restrict(x, tbl)), tbl) for x in facttbl]
      data = data[string(p)]
      for l in data
        res[mp[l[2]]] = l[4]
      end
      # Set values for faithful characters if applicable.
      perf = (count(chi -> degree(chi) == 1, tbl) == 1)
      if perf && order(tbl) == 2*order(facttbl) && p != 2
        # deal with the faithful characters of a perfect central extension 2.G,
        # the image of the spinor norm is trivial
        for i in 1:length(tbl)
          chi = tbl[i]
          deg = degree(chi)
          ind = indicator(chi)
          if ind == 1 && length(class_positions_of_kernel(chi)) == 1
            if mod(deg, 4) == 0
              if p == 0
                res[i] = "1"
              else
                res[i] = "O+"
              end
            elseif mod(deg, 2) == 0
              if p == 0
                res[i] = "-1"
              else
                # Check whether -1 is a square in the character field.
                if mod(p-1, 4) == 0 || mod(degree_of_character_field(chi), 2) == 0
                  res[i] = "O+"
                else
                  res[i] = "O-"
                end
              end
            end
          end
        end
      end
    end
  end

  # Fill missing values.
  for i in 1:length(tbl)
    if res[i] == "" && indicator(tbl[i]) == 1 && is_even(numerator(degree(tbl[i])))
      res[i] = "?"
    end
  end
  return res
end


# Return the character described by `d`.
function character_of_entry(d::Dict)
  @req haskey(d, :groupname) "the dictionary has no :groupname"
  tbl = character_table(d[:groupname])
  @req haskey(d, :characteristic) "the dictionary has no :characteristic"
  p = d[:characteristic]
  if p != 0
    tbl = mod(tbl, p)
  end
  @req haskey(d, :charpos) "the dictionary has no :charpos"
  return tbl[d[:charpos]]
end


# Return `true` if `str` occurs as an entry in `d[:comment]`
function comment_matches(d::Dict, str::String)
  @req haskey(d, :comment) "the dictionary has no :comment"
  return str in d[:comment]
end


# Compare two character fields, by comparing their embeddings.
function is_equal_field(emb1, emb2)
  dom1 = domain(emb1)
  dom2 = domain(emb2)
  p = characteristic(dom1)
  p == characteristic(dom2) || return false
  degree(dom1) == degree(dom2) || return false
  p == 0 || return order(dom1) == order(dom2)
  return has_preimage_with_preimage(emb2, emb1(gen(dom1)))[1]
end


##############################################################################

const _OD_filter_attrs = Dict{Any,Tuple{Type, Any, Any}}()

function __init_OD()
  empty!(_OD_filter_attrs)

  props = [
    is_simple,
    is_sporadic_simple,
  ]

  for k in props
    _OD_filter_attrs[k] = (Bool, k, true)
    _OD_filter_attrs[!k] = (Bool, k, false)
  end

  _OD_filter_attrs[characteristic] =
    (Union{Oscar._IntOrIntVec, Function},
    characteristic, nothing)
  _OD_filter_attrs[character_field] =
    (Union{Map, Vector{Map}, Function},
    character_field, nothing)
  _OD_filter_attrs[degree] =
    (Union{Oscar._IntOrIntVec, Function},
    degree, nothing)
  _OD_filter_attrs[identifier] =
    (Union{String, Vector{String}, Function},
    identifier, nothing)
  _OD_filter_attrs[dim] =
    (Union{Oscar._IntOrIntVec, Function},
    dim, nothing)
  _OD_filter_attrs[orthogonal_discriminant] =
    (Union{String, Vector{String}, Function},
    orthogonal_discriminant, nothing)
  _OD_filter_attrs[comment_matches] =
    (Union{String, Vector{String}, Function},
    comment_matches, nothing)
end

__init_OD()


function _od_info(groupname, p, v)
  return Dict{Symbol, Any}(
           :groupname => groupname,
           :characteristic => p,
           :charname => v[1],
           :charpos => v[2],
           :degree => v[3],
           :valuestring => v[4],
           :comment => v[end])
end


@doc raw"""
    all_od_infos(L...)

Return the array of all those entries of the known OD data (see `OD_data`)
that satisfy the conditions in `L`.

The following conditions are supported.

- `is_simple` with value `true` or `false`,
  meaning entries only for simple or non-simple groups, respectively,

- `is_sporadic_simple` with value `true` or `false`,
  meaning entries only for sporadic simple or not sporadic simple groups,
  respectively,

- `characteristic`, with value `0` or a prime integer,
  meaning entries only for this characteristic,

- `character_field`, with value either a map or a vector of maps,
  meaning entries only for characters whose character fields
  (finite fields if the characteristic is positive, and subfields of
  cyclotomic fields in characteristic zero) have the given map(s) as
  embeddings into the algebraic closure or abelian closure, respectively,

- `degree`, with value a positive integer,
  meaning entries only for characters whose character field has the given
  degree over its prime field,

- `identifier`, with value a string denoting the name of an Atlas group,
  or a vector of such strings,
  meaning entries only for these groups,

- `dim`, with value a positive integer,
  or a vector of such integers,
  meaning entries only for characters of these degrees,

- `orthogonal_discriminant`, with value a string (`"O+"`, `"O-"`, or a string
  that encodes an algebraic integer),
  meaning entries only with this orthogonal discriminant,

- `comment_matches`, with value a string (one of `"ev"`, `"specht"`, ...),
  or a vector of such strings,
  meaning entries whose comment contains these values.

For all conditions except the boolean valued ones, also a function can be
given as value, meaning that all those entries satisfy this condition
for which the function returns `true` when applied to the stored value.
For example, the condition `characteristic => is_odd` matches all entries
for characteristics different from `0` and `2`,
and the condition `character_field => (emb -> degree(domain(emb)) == 1)`
matches all entries for which the character field is the field of rationals.

# Examples
```jldoctest
julia> length(all_od_infos(identifier => "A6"))
8

julia> length(all_od_infos(identifier => "A6", characteristic => 0))
3

julia> length(all_od_infos(identifier => "A6", characteristic => 2:5))
5
```
"""
function all_od_infos(L...)
  iso = nothing

  # scan the given conditions
  conditions = IdDict()

  res = Dict{Symbol, Any}[]

  gapargs = Any[]
  for arg in L
    if arg isa Pair
      # handle e.g.
      # `dimension => 4`, `dimension => [4, 8]`, `dimension => iseven`
      func = arg[1]
      data = arg[2]
      @req haskey(_OD_filter_attrs, func) "Function not supported"
      expected_type, func1, _ = _OD_filter_attrs[func]
      @req data isa expected_type "bad argument $(data) for function $(func)"
      conditions[func1] = data
    elseif arg isa Function
      # handle e.g. `is_simple` or `! is_simple`
      func = arg
      @req haskey(_OD_filter_attrs, func) "Function not supported"
      expected_type, func1, default = _OD_filter_attrs[func]
      @req default !== nothing "missing argument for function $(func)"
      conditions[func1] = default
    else
      throw(ArgumentError("expected a function or a pair, got $arg"))
    end
  end

  # Evaluate the conditions:
  # `is_simple` and `is_sporadic_simple`
  names = nothing
  simp = nothing
  if haskey(conditions, is_simple)
    simp = conditions[is_simple]
    if simp
      names = [(x, x) for x in OD_data["names"]]
    else
      names = Tuple{String, String}[]
      for simpnam in OD_data["names"]
        info = OD_data[simpnam]
        append!(names, [(simpnam, x) for x in info["names"][2:end]])
      end
    end
  end

  if haskey(conditions, is_sporadic_simple)
    sporsimp = conditions[is_sporadic_simple]
    spor_names = ["M11", "M12", "J1", "M22", "J2", "M23", "HS", "J3", "M24",
                  "McL", "He", "Ru", "Suz", "ON", "Co3", "Co2", "Fi22", "HN",
                  "Ly", "Th", "Fi23", "Co1", "J4", "Fi24'", "B", "M"]
    spor_names = filter(x -> haskey(OD_simple_names, x), spor_names)
    if sporsimp
      simp === false && return res
      names = [(x, x) for x in spor_names]
    elseif simp !== nothing
      if simp
        names = filter(x -> !(x[1] in spor_simp), names)
      end
    else
      names = Tuple{String, String}[]
      for simpnam in OD_data["names"]
        info = OD_data[simpnam]
        if simpnam in spor_names
          append!(names, [(simpnam, x) for x in info["names"][2:end]])
        else
          append!(names, [(simpnam, x) for x in info["names"]])
        end
      end
    end
  end

  if haskey(conditions, identifier)
    ids = conditions[identifier]
    if ids isa String
      haskey(OD_simple_names, ids) || return res
      ids = (OD_simple_names[ids], ids)
      names !== nothing && !(ids in names) && return []
      names = [ids]
    elseif ids isa Vector{String}
      if names !== nothing
        names = filter(in(names), [(OD_simple_names[x], x) for x in ids])
      else
        names = [(OD_simple_names[x], x) for x in ids]
      end
    elseif ids isa Function
      if names !== nothing
        names = filter(x -> ids(x[2]), names)
      else
        names = Tuple{String, String}[]
        for simpnam in OD_data["names"]
          info = OD_data[simpnam]
          append!(names, [(simpnam, x) for x in filter(ids, info["names"])])
        end
      end
    end
  elseif names === nothing
    names = Tuple{String, String}[]
    for simpnam in OD_data["names"]
      info = OD_data[simpnam]
      append!(names, [(simpnam, x) for x in info["names"]])
    end
  end

  for (simpnam, nam) in names
    D = OD_data[simpnam][nam]
    for char in sort!([parse(Int, x) for x in keys(D)])
      good_char = true
      if haskey(conditions, AbstractAlgebra.characteristic)
        characteristic = conditions[AbstractAlgebra.characteristic]
        if characteristic isa Function
          (!characteristic(char)) && (good_char = false)
        elseif characteristic isa IntegerUnion
          (characteristic != char) && (good_char = false)
        elseif !(char in characteristic)
          good_char = false
        end
      end
      good_char || continue
      for entry in D[string(char)]
        good_dim = true
        if haskey(conditions, Hecke.dim)
          dim = parse(Int, filter(isdigit, entry[1]))
          dimension = conditions[Hecke.dim]
          if dimension isa Function
            (!dimension(dim)) && (good_dim = false)
          elseif dimension isa IntegerUnion
            (dimension != dim) && (good_dim = false)
          elseif !(dim in dimension)
            good_dim = false
          end
        end
        good_dim || continue

        good_degree = true
        if haskey(conditions, degree)
          d = entry[3]
          deg = conditions[degree]
          if deg isa Function
            (!deg(d)) && (good_degree = false)
          elseif deg isa IntegerUnion
            (deg != d) && (good_degree = false)
          elseif !(d in deg)
            good_degree = false
          end
        end
        good_degree || continue

        good_od = true
        if haskey(conditions, Oscar.OrthogonalDiscriminants.orthogonal_discriminant)
          od = entry[4]
          orthogonal_discriminant = conditions[Oscar.OrthogonalDiscriminants.orthogonal_discriminant]
          if orthogonal_discriminant isa Function
            (!orthogonal_discriminant(od)) && (good_od = false)
          elseif orthogonal_discriminant isa String
            (orthogonal_discriminant != od) && (good_od = false)
          elseif !(od in orthogonal_discriminant)
            good_od = false
          end
        end
        good_od || continue

        good_comment = true
        if haskey(conditions, Oscar.OrthogonalDiscriminants.comment_matches)
          comment = entry[end]
          comment_matches = conditions[Oscar.OrthogonalDiscriminants.comment_matches]
          if comment_matches isa Function
            (!comment_matches(comment)) && (good_comment = false)
          elseif comment_matches isa String
            (!(comment_matches in comment)) && (good_comment = false)
          elseif is_empty(intersect(comment, comment_matches))
            good_comment = false
          end
        end
        good_comment || continue

        good_field = true
        if haskey(conditions, character_field)
          f = conditions[character_field]
          ordt = character_table(nam)
          t = (char == 0 ? ordt : mod(ordt, char))
          F, emb = character_field(t[entry[2]])
          if f isa Function
            (!f(emb)) && (good_field = false)
          elseif f isa Map
            (!is_equal_field(f, emb)) && (good_field = false)
          elseif all(x -> ! is_equal_field(x, emb), f)
            good_field = false
          end
        end
        good_field || continue

        push!(res, _od_info(nam, char, entry))
      end
    end
  end

  return res
end
