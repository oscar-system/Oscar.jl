##  This is a first attempt to implement group characters in Oscar.
##  
##  The idea is that the available GAP objects (groups, character tables,
##  class functions) are used in a first step, and that access to character
##  values yields `QQAbElem` objects.
##  
##  Once we agree on the functionality and the integration into Oscar,
##  this setup can in a second step be replaced by one that uses
##  native Julia objects for representing class functions,
##  but character tables and groups still have some counterpart in GAP.
##  
##  In a third step, we replace the character table objects by native Julia
##  objects.

# character values are elements from QQAbField


#############################################################################
##
##  Atlas irrationalities
##
"""
    atlas_irrationality([F::AnticNumberField, ]description::String)

Return the value encoded by `description`.
If `F` is given and is a cyclotomic field that contains the value then
the result is in `F`,
if `F` is not given then the result has type `QQAbElem`.

`description` is assumed to have the format defined in
[CCNPW85](@cite), Chapter 6, Section 10.

```jldoctest
julia> Oscar.with_unicode() do
         show(atlas_irrationality("r5"))
       end;
-2*ζ(5)^3 - 2*ζ(5)^2 - 1

julia> atlas_irrationality(CyclotomicField(5)[1], "r5")
-2*z_5^3 - 2*z_5^2 - 1

julia> Oscar.with_unicode() do
         show(atlas_irrationality("i"))
       end;
ζ(4)

julia> Oscar.with_unicode() do
         show(atlas_irrationality("b7*3"))
       end;
-ζ(7)^4 - ζ(7)^2 - ζ(7) - 1

julia> Oscar.with_unicode() do
         show(atlas_irrationality("3y'''24*13-2&5"))
       end;
-5*ζ(24)^7 - 2*ζ(24)^5 + 2*ζ(24)^3 - 3*ζ(24)
```
"""
function atlas_irrationality(F::AnticNumberField, description::String)
    return F(GAP.Globals.AtlasIrrationality(GapObj(description))::GAP.Obj)
end

function atlas_irrationality(description::String)
    F = abelian_closure(QQ)[1]
    return F(GAP.Globals.AtlasIrrationality(GapObj(description))::GAP.Obj)
end


#############################################################################
##
##  character tables
##
abstract type GroupCharacterTable end

"""
    GAPGroupCharacterTable <: GroupCharacterTable

This is the type of (ordinary or Brauer) character tables that can delegate
tasks to an underlying character table object in the GAP system
(field `GAPTable`).

The value of the field `characteristic` determines whether the table
is an ordinary one (value `0`) or a `p`-modular one (value `p`).

A group can (but need not) be stored in the field `group`.
If it is available then also the field `isomorphism` is available,
its value is a bijective map from the `group` value to a group in GAP.

Objects of type `GAPGroupCharacterTable` support [`get_attribute`](@ref),
for example in order to store the already computed `p`-modular tables
in an ordinary table, and to store the corresponding ordinary table
in a `p`-modular table.
"""
@attributes mutable struct GAPGroupCharacterTable <: GroupCharacterTable
    GAPTable::GapObj  # the GAP character table object
    characteristic::Int
    group::Union{GAPGroup, GrpAbFinGen}    # the underlying group, if any
    isomorphism::Map  # isomorphism from `group` to a group in GAP

    function GAPGroupCharacterTable(G::Union{GAPGroup, GrpAbFinGen}, tab::GapObj, iso::Map, char::Int)
      return new(tab, char, G, iso)
    end

    function GAPGroupCharacterTable(tab::GapObj, char::Int)
      #group is left undefined
      return new(tab, char)
    end
end

# access to field values via functions
GAPTable(tbl::GAPGroupCharacterTable) = tbl.GAPTable
group(tbl::GAPGroupCharacterTable) = tbl.group
characteristic(tbl::GAPGroupCharacterTable) = tbl.characteristic

# backwards compatibility:
# return `tbl.group` when one accesses `tbl.GAPGroup`
function Base.getproperty(tbl::GAPGroupCharacterTable, sym::Symbol)
   sym === :GAPGroup && return tbl.group
   return getfield(tbl, sym)
end


# If the character table stores an Oscar group `G` then
# an isomorphism from `G` to a GAP group is stored in the field `isomorphism`.
# Different kinds of isomorphisms can be stored,
# depending on the type of the stored group.
#
# - For `GAPGroup`, we compute images and preimages by unwrapping and
#   wrapping the groups elements, respectively.
# - For `GrpAbFinGen`, we compose this unwrapping/wrapping with
#   an isomorphism to a `PcGroup` .
#
function isomorphism_to_GAP_group(G::GAPGroup)
    f = function(x) return x.X; end
    finv = function(x::GAP.Obj) return group_element(G, x); end
    return MapFromFunc(f, finv, G, G.X)
end

function isomorphism_to_GAP_group(G::GrpAbFinGen)
    @req isfinite(G) "the group is not finite"
    iso = isomorphism(PcGroup, G)
    C = codomain(iso)
    @assert C isa GAPGroup
    f = function(x) return iso(x).X; end
    finv = function(x::GAP.Obj) return preimage(iso, group_element(C, x)); end
    return MapFromFunc(f, finv, G, C.X)
end

function isomorphism_to_GAP_group(tbl::GAPGroupCharacterTable)
    @req isdefined(tbl, :isomorphism) "character table stores no isomorphism"
    return tbl.isomorphism
end


# If the character table stores an Oscar group `G` then
# the columns correspond to conjugacy classes of `G`.
# The vector `conjugacy_classes(G)` (stored or not stored in `G`) can be
# independent of the vector of conjugacy classes of the corresponding
# group on the GAP side.
# Therefore we store a vector `conjugacy_classes(tbl)` via an attribute,
# which is compatible with the `GAP.Globals.ConjugacyClasses` value of
# the underlying GAP character table.
@attr function conjugacy_classes(tbl::GAPGroupCharacterTable)
    @req isdefined(tbl, :group) "character table stores no group"

    # If `GAPTable(tbl)` does not yet store conjugacy classes
    # then compute them.
    if characteristic(tbl) == 0
      ccl = GAP.Globals.ConjugacyClasses(GAPTable(tbl))::GapObj
      ccl = [conjugacy_class(group(tbl),
               preimage(isomorphism_to_GAP_group(tbl),
                        GAPWrap.Representative(x)))
             for x in ccl]
    else
      ordtbl = get_attribute(tbl, :ordinary_table)
      ccl = conjugacy_classes(ordtbl)
      known, fus = known_class_fusion(tbl, ordtbl)
      @assert known "the class fusion is not stored"
      ccl = ccl[fus]
    end
    return ccl
end


"""
    character_table(G::GAPGroup, p::Int = 0)

Return the ordinary (if `p == 0`) or `p`-modular character table of the
finite group `G`.
If the `p`-modular character table of `G` cannot be computed by GAP
then `nothing` is returned.

# Examples
```jldoctest
julia> Oscar.with_unicode() do
         show(character_table(symmetric_group(3)))
       end;
Sym( [ 1 .. 3 ] )

 2  1  1  .
 3  1  .  1
           
   1a 2a 3a
2P 1a 1a 3a
3P 1a 2a 1a
           
χ₁  1 -1  1
χ₂  2  . -1
χ₃  1  1  1

julia> Oscar.with_unicode() do
         show(character_table(symmetric_group(3), 2))
       end;
Sym( [ 1 .. 3 ] ) mod 2

 2  1  .
 3  1  1
        
   1a 3a
2P 1a 3a
3P 1a 1a
        
χ₁  1  1
χ₂  2 -1
```
"""
function character_table(G::Union{GAPGroup, GrpAbFinGen}, p::Int = 0)
    tbls = get_attribute!(() -> Dict{Int,Any}(), G, :character_tables)
    return get!(tbls, p) do
      if p == 0
        iso = isomorphism_to_GAP_group(G)
        gaptbl = GAP.Globals.CharacterTable(codomain(iso))::GapObj
      else
        # Create the `p`-modular table if possible.
        @req is_prime(p) "p must be 0 or a prime integer"
        ordtbl = character_table(G, 0)
        iso = isomorphism_to_GAP_group(ordtbl)
        gaptbl = GAP.Globals.mod(GAPTable(ordtbl), GAP.Obj(p))::GapObj
        gaptbl === GAP.Globals.fail && return nothing
      end
      return GAPGroupCharacterTable(G, gaptbl, iso, p)
    end
end

# A character table with stored group object is stored in this group.
# Character tables from the table library do not store groups,
# they are cached in the dictionary `character_tables_by_id`,
# in order to achieve that fetching the same table twice yields the same
# object.
const character_tables_by_id = Dict{String, Union{GAPGroupCharacterTable, Nothing}}()

"""
    character_table(id::String, p::Int = 0)

Return the ordinary (if `p == 0`) or `p`-modular character table
for which `id` is an admissible name in GAP's library of character tables.
If no such table is available then `nothing` is returned.

# Examples
```jldoctest
julia> println(character_table("A5"))
character_table("A5")

julia> println(character_table("A5", 2))
character_table("A5mod2")

julia> println(character_table("J5"))
nothing
```
"""
function character_table(id::String, p::Int = 0)
    if p == 0
      modid = id
    else
      @req is_prime(p) "p must be 0 or a prime integer"
      modid = "$(id)mod$(p)"
    end

    return get!(character_tables_by_id, modid) do
      tbl = GAP.Globals.CharacterTable(GapObj(modid))::GapObj
      tbl === GAP.Globals.fail && return nothing
      return GAPGroupCharacterTable(tbl, p)
    end
end


"""
    character_table(series::Symbol, parameter::Any)

Return the ordinary character table of the group described by the series
`series` and the parameter `parameter`.

# Examples
```jldoctest
julia> println(character_table(:Symmetric, 5))
character_table("Sym(5)")

julia> println(character_table(:WeylB, 3))
character_table("W(B3)")
```

Currently the following series are supported.

| Series | Parameter |
| ------ | ---------------- |
| `:Cyclic` | pos. integer |
| `:Dihedral` | even pos. integer |
| `:Symmetric` | pos. integer |
| `:Alternating` | integer `> 1` |
| `:WeylB` | pos. integer |
| `:WeylD` | integer `> 1` |
| `:DoubleCoverSymmetric` | pos. integer |
| `:DoubleCoverAlternating` | pos. integer |
| `:GL2` | prime power |
| `:SL2odd` | odd prime power |
| `:SL2even` | even prime power |
| `:PSL2odd` | odd prime power `q` s. t. `(q-1)/2` is odd |
| `:PSL2even` | odd prime power `q` s. t. `(q-1)/2` is even |
| `:Suzuki` | odd power of 2 |
| `:GU3` | prime power |
| `:SU3` | prime power |
| `Symbol("P:Q")` | array `[p, q]` with prime `p` and `q` dividing `p-1` |
| `:ExtraspecialPlusOdd` | odd power of odd prime |
"""
function character_table(series::Symbol, parameter::Union{Int, Vector{Int}})
    args = GAP.Obj([string(series), parameter], recursive = true)
    tbl = GAP.Globals.CallFuncList(GAP.Globals.CharacterTable, args)::GapObj
    tbl === GAP.Globals.fail && return nothing
    tbl = GAPGroupCharacterTable(tbl, 0)
    set_attribute!(tbl, :series, (series, parameter))
    return tbl
end


# Assume that `isomorphism_to_GAP_group(tbl)` maps elements via `.X` access.
function preimages(iso::MapFromFunc{T, GapObj}, H::GapObj) where T <: GAPGroup
    return _as_subgroup(domain(iso), H)
end

# Use the isomorphism.
function preimages(iso::MapFromFunc{GrpAbFinGen, GapObj}, H::GapObj)
  return sub(domain(iso), [preimage(iso, x) for x in GAP.Globals.GeneratorsOfGroup(H)])
end


# For character tables with stored group, we take the hash value of the group.
# For character tables without stored group, we take the table identifier.
function Base.hash(tbl::GAPGroupCharacterTable, h::UInt)
  if isdefined(tbl, :group)
    return Base.hash(group(tbl), h)
  else
    return Base.hash(identifier(tbl), h)
  end
end


##############################################################################
#
# admissible names of library character tables

"""
    is_duplicate_table(tbl::GAPGroupCharacterTable)

Return whether `tbl` is a table from the character table library
that was constructed from another library character table by permuting rows
and columns.

One application of this function is to restrict the search with
[`all_character_table_names`](@ref) to only one library character table for each
class of permutation equivalent tables.
"""
@gapattribute is_duplicate_table(tbl::GAPGroupCharacterTable) = GAP.Globals.IsDuplicateTable(GAPTable(tbl))::Bool

"""
    all_character_table_names(L...; ordered_by = nothing)

Return an array of strings that contains all those names of character tables
in the character table library that satisfy the conditions in the array `L`.

# Examples
```
julia> spor_names = all_character_table_names(is_sporadic_simple => true,
         is_duplicate_table => false);

julia> println(spor_names[1:5])
["B", "Co1", "Co2", "Co3", "F3+"]

julia> spor_names = all_character_table_names(is_sporadic_simple,
         !is_duplicate_table; ordered_by = order);

julia> println(spor_names[1:5])
["M11", "M12", "J1", "M22", "J2"]

julia> length(all_character_table_names(number_conjugacy_classes => 1))
1
```
"""
function all_character_table_names(L...; ordered_by = nothing)
    gapargs = translate_group_library_args(L; filter_attrs = _ctbl_filter_attrs)

    if ordered_by isa Function
      K = GAP.call_gap_func(GAP.Globals.AllCharacterTableNames, gapargs...;
            OrderedBy = find_index_function(ordered_by, _ctbl_filter_attrs)[2])::GapObj
    else
      K = GAP.Globals.AllCharacterTableNames(gapargs...)::GapObj
    end
    return Vector{String}(K)
end


##############################################################################
#
# `print` and `show` character tables

# Utility:
# Create strings in length-lexicographical ordering w.r.t. the
# alphabet 'alphabet'.
# (If `alphabet` is `"ABCDEFGHIJKLMNOPQRSTUVWXYZ"` then the strings
# have the form `"A", "B", ..., "Z", "AA", ...`.)
mutable struct WordsIterator
    alphabet::String
end

Base.iterate(wi::WordsIterator) = length(wi.alphabet) == 0 ? nothing : (string(wi.alphabet[1]), 2)

function Base.iterate(wi::WordsIterator, state::Int)
    name = ""
    n = state
    ll = length(wi.alphabet)
    while 0 < n
      n, r = divrem(n-1, ll)
      name = wi.alphabet[r+1] * name
    end
    return (name, state+1)
end


@doc raw"""
    as_sum_of_roots(val::nf_elem, root::String)

Return a string representing the element `val` of a cyclotomic field
as a sum of multiples of powers of the primitive root which is printed as
`root`.
"""
function as_sum_of_roots(val::nf_elem, root::String)
    F = parent(val)
    flag, N = Hecke.is_cyclotomic_type(F)
    @req flag "$val is not an element of a cyclotomic field"

    # `string` yields an expression of the right structure,
    # but not in terms of `root`
    # and without curly brackets for subscripts and superscripts.
    str = string(val)
    str = replace(str, "*" => "")
    str = replace(str, string(F.S) => "$(root)_{$N}")
    str = replace(str, r"\^([0-9]*)" => s"^{\1}")
    return str
end


@doc raw"""
    matrix_of_strings(tbl::GAPGroupCharacterTable; alphabet::String = "", root::String = "\\zeta")

Return `(mat, legend)` where `mat` is a matrix of strings that describe
the values of the irreducible characters of `tbl`,
and `legend` is a list of triples `(val, name, disp)` such that
`name` occurs in `mat`,
the character value `val` is represented by `name` in `mat`,
and `disp` is a string that describes `val`.

Integral character values are just turned into strings, except that `0` gets
replaced by `"."`.

If `alphabet` is empty then irrational values are written as sums of powers
of roots of unity,
where $\exp(2 \pi i/n)$ is shown as `root` followed by the index $n$.
Since all entries of the matrix are self-explanatory, `legend` is empty
in this case.

If `alphabet` is nonempty then irrational values are written as words in terms
of `alphabet`, and these words are the `name` entries in `legend`;
the corresponding `disp` entries are sums of powers of roots of unity
in terms of `root`.
"""
function matrix_of_strings(tbl::GAPGroupCharacterTable; alphabet::String = "", root::String = "\\zeta")
  n = nrows(tbl)
  m = Array{String}(undef, n, n)
  legend = []
  if alphabet != ""
    iter = WordsIterator(alphabet)
    state = 1
  end

  # Run column-wise through the matrix of irreducibles.
  # The same is done in GAP, thus the names are compatible with the ones
  # created by GAP,
  # except that relative names in the case of complex conjugation and
  # quadratic irrationalities are currently not handled here.
  for j in 1:n
    for i in 1:n
      val = tbl[i,j]
      if iszero(val)
        m[i, j] = "."
      elseif val.c == 1
        m[i, j] = string(val.data)
      elseif alphabet != ""
        # write irrationalities using symbolic names
        pos = findnext(x -> x[1] == val, legend, 1)
        if pos != nothing
          m[i,j] = legend[pos][2]
        else
          pos = findnext(x -> x[1] == -val, legend, 1)
          if pos != nothing
            # The negative of a known name is shown relative to that name.
            m[i,j] = "-" * legend[pos][2]
          else
            name, state = iterate(iter, state)
            disp = as_sum_of_roots(val.data, root)
            push!(legend, (val, name, disp))
            valbar = conj(val)
            if valbar != val
              # The complex conjugate of a known name is shown relative to
              # that name.
              disp = as_sum_of_roots(valbar.data, root)
              push!(legend, (valbar, "\\overline{$name}", disp))
            else
              info = AbelianClosure.quadratic_irrationality_info(val)
              if info != nothing
                # `val` generates a quadratic field extension.
                # Show the unique Galois conjugate of `A` different from
                # `A` as `A*`, and show both `A` and `A*` in the footer.
#TODO: For elements in quadratic fields, show also an expression in terms of
#      square roots in the footer.
                valstar = 2*info[1] - val
                disp = as_sum_of_roots(valstar.data, root)
                push!(legend, (valstar, name*"*", disp))
              end
            end
            m[i,j] = name
          end
        end
      else
        # write irrationalities in terms of `root`.
        m[i,j] = as_sum_of_roots(val.data, root)
      end
    end
  end
  return (m, legend)
end

# Produce LaTeX output if `"text/html"` is prescribed,
# via the `:TeX` attribute of the io context.
function Base.show(io::IO, ::MIME"text/latex", tbl::GAPGroupCharacterTable)
  print(io, "\$")
  show(IOContext(io, :TeX => true), tbl)
  print(io, "\$")
end

# Produce a screen format without LaTeX markup but with unicode characters
# and sub-/superscripts if LaTeX output is not requested..
function Base.show(io::IO, tbl::GAPGroupCharacterTable)
    n = nrows(tbl)
    gaptbl = GAPTable(tbl)
    size = order(ZZRingElem, tbl)
    primes = [x[1] for x in collect(factor(size))]
    sort!(primes)

    # Decide how to deal with irrationalities.
    alphabet = get(io, :alphabet, "")
    if alphabet == ""
      with_legend = get(io, :with_legend, false)
      if with_legend == true
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      end
    end

    # Create the strings of the values of the irreducibles.
    mat, legend = matrix_of_strings(tbl, alphabet = alphabet)

    # Compute the factored centralizer orders.
    cents = orders_centralizers(tbl)
    fcents = [collect(factor(x)) for x in cents]
    d = Dict([p => fill(".", n) for p in primes]...)
    for i in 1:n
      for pair in fcents[i]
        d[pair[1]][i] = string(pair[2])
      end
    end
    cents_strings = [d[p] for p in primes]

    # Compute display format for power maps.
    names = Vector{String}(GAP.Globals.ClassNames(gaptbl)::GapObj)
    pmaps = Vector{Any}(GAP.Globals.ComputedPowerMaps(gaptbl)::GapObj)
    power_maps_primes = String[]
    power_maps_strings = Vector{String}[]
    for i in 2:length(pmaps)
      map = pmaps[i]
      if map != nothing
        push!(power_maps_primes, string(i)*"P")
        push!(power_maps_strings, names[map])
      end
    end

    # Compute the indicator values if applicable.
    ind = get(io, :indicator, Int[])
    if ind == true
      ind = [2]
    end
    indicators = [[string(indicator(x, n)) for x in tbl] for n in ind]
    for i in 1:length(ind)
      if ind[i] == 2
        indicators[i] = replace( x -> x == "1" ? "+" :
                                    ( x == "0" ? "o" :
                                    ( x == "-1" ? "-" : "?" ) ), indicators[i])
      end
    end
    emptycor = ["" for i in 1:length(ind)]

    # Fetch the Orthogonal Discriminants if applicable.
    # (This is possible only if the OD database is available.
    OD = get(io, :OD, false)::Bool
    if OD && hasproperty(GAP.Globals, :OrthogonalDiscriminants)
      ODs = [replace(x -> isnothing(x) ? "" : string(x),
                     Vector{Any}(GAP.Globals.OrthogonalDiscriminants(gaptbl)::GapObj))]
      ODlabel = ["OD"]
      push!(emptycor, "")
    else
      ODs = []
      ODlabel = []
    end

    # Compute the degrees of the character fields if applicable.
    field_degrees = get(io, :character_field, false)::Bool
    if field_degrees
      p = characteristic(tbl)
      if p == 0
        field_degrees = [[string(GAP.Globals.Dimension(GAP.Globals.Field(GAP.Globals.Rationals, x.values))) for x in tbl]]
      else
        field_degrees = [[string(collect(factor(GAP.Globals.SizeOfFieldOfDefinition(x.values, p)))[1][2]) for x in tbl]]
      end
      field_label = ["d"]
      push!(emptycor, "")
    else
      field_degrees = []
      field_label = []
    end

    emptycol = ["" for i in 1:n]

    if isdefined(tbl, :group)
      headerstring = string(group(tbl))
      if characteristic(tbl) != 0
        headerstring = "$headerstring mod $(characteristic(tbl))"
      end
    else
      headerstring = identifier(tbl)
    end

    # Create the IO context.
    ioc = IOContext(io,
      # header (an array of strings):
      # name of the table and separating empty line
      :header => [headerstring, ""],

      # column labels:
      # centralizer orders (factored),
      # separating empty line,
      # class names,
      # separating empty line,
      # p-th power maps for known p-th power maps,
      # separating empty line,
      :labels_col => permutedims(hcat(
        cents_strings..., emptycol, names, power_maps_strings..., emptycol)),

      # row labels:
      # character names and perhaps indicators.
      :labels_row => hcat(["\\chi_{$i}" for i in 1:n], field_degrees..., ODs..., indicators...),

      # corner:
      # primes in the centralizer rows,
      # separating empty line,
      # separating empty line,
      # primes in the power map rows,
      # separating line perhaps containing indicator labels,
      :corner => permutedims(hcat(
                      [vcat(emptycor, [string(x)]) for x in primes]...,
                      vcat(emptycor, [""]),
                      vcat(emptycor, [""]),
                      [vcat(emptycor, [x]) for x in power_maps_primes]...,
                      vcat([""], field_label, ODlabel, [string(x) for x in ind]))),

      # footer (an array of strings)
      :footer => length(legend) == 0 ? [] :
                 vcat([""], [triple[2]*" = "*triple[3] for triple in legend]),
    )

    # print the table
    labelled_matrix_formatted(ioc, mat)
end

# print: abbreviated form
function Base.print(io::IO, tbl::GAPGroupCharacterTable)
    gaptbl = GAPTable(tbl)
    if isdefined(tbl, :group)
      id = string(group(tbl))
      if characteristic(tbl) != 0
        id = "$(id)mod$(characteristic(tbl))"
      end
    else
      id = "\"" * identifier(tbl) * "\""
    end
    print(io, "character_table($id)")
end


##############################################################################
#
length(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(GAPTable(tbl))::Int
nrows(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(GAPTable(tbl))::Int
ncols(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(GAPTable(tbl))::Int
number_conjugacy_classes(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(GAPTable(tbl))::Int

@doc raw"""
    order(::Type{T} = ZZRingElem, tbl::GAPGroupCharacterTable) where T <: IntegerUnion

Return the order of the group for which `tbl` is the character table.

# Examples
```jldoctest
julia> order(character_table(symmetric_group(4)))
24
```
"""
order(tbl::GAPGroupCharacterTable) = order(ZZRingElem, tbl)

function order(::Type{T}, tbl::GAPGroupCharacterTable) where T <: IntegerUnion
  return T(GAPWrap.Size(GAPTable(tbl)))
end

@doc raw"""
    orders_class_representatives(tbl::GAPGroupCharacterTable)

Return the array of the orders of conjugacy class representatives for `tbl`,
ordered according to the columns of `tbl`.

# Examples
```jldoctest
julia> println(orders_class_representatives(character_table("A5")))
[1, 2, 3, 5, 5]
```
"""
@gapattribute orders_class_representatives(tbl::GAPGroupCharacterTable) = Vector{Int}(GAP.Globals.OrdersClassRepresentatives(GAPTable(tbl))::GapObj)

@doc raw"""
    orders_centralizers(tbl::GAPGroupCharacterTable)

Return the array of the orders of centralizers of conjugacy class
representatives for `tbl` in the group of `tbl`,
ordered according to the columns of `tbl`.

# Examples
```jldoctest
julia> println(orders_centralizers(character_table("A5")))
ZZRingElem[60, 4, 3, 5, 5]
```
"""
@gapattribute orders_centralizers(tbl::GAPGroupCharacterTable) = Vector{ZZRingElem}(GAP.Globals.SizesCentralizers(GAPTable(tbl))::GAP.Obj)

@doc raw"""
    class_lengths(tbl::GAPGroupCharacterTable)

# Examples
```jldoctest
julia> println(class_lengths(character_table("A5")))
ZZRingElem[1, 15, 20, 12, 12]
```
"""
@gapattribute class_lengths(tbl::GAPGroupCharacterTable) = Vector{ZZRingElem}(GAP.Globals.SizesConjugacyClasses(GAPTable(tbl))::GapObj)

@doc raw"""
    maxes(tbl::GAPGroupCharacterTable)

Return either nothing (if the value is not known) or an array of identifiers
of the ordinary character tables of all maximal subgroups of `tbl`.
There is no default method to compute this value from `tbl`.

If the `maxes` value of `tbl` is stored then it lists exactly one
representative for each conjugacy class of maximal subgroups of the group of
`tbl`, and the character tables of these maximal subgroups are available
in the character table library, and compatible class fusions to `tbl` are
stored on these tables.

# Examples
```jldoctest
julia> println(maxes(character_table("M11")))
["A6.2_3", "L2(11)", "3^2:Q8.2", "A5.2", "2.S4"]

julia> maxes(character_table("M")) === nothing  # not (yet) known
true
```
"""
function maxes(tbl::GAPGroupCharacterTable)
  if GAP.Globals.HasMaxes(GAPTable(tbl))::Bool
    return Vector{String}(GAP.Globals.Maxes(GAPTable(tbl))::GapObj)
  end
  return nothing
end

@doc raw"""
    identifier(tbl::GAPGroupCharacterTable)

Return a string that identifies `tbl`.
It is used mainly for library tables.

# Examples
```jldoctest
julia> identifier(character_table("A5"))
"A5"
```
"""
@gapattribute identifier(tbl::GAPGroupCharacterTable) = string(GAP.Globals.Identifier(GAPTable(tbl))::GapObj)

@doc raw"""
    class_positions_of_pcore(tbl::GAPGroupCharacterTable, p::IntegerUnion)

Return the array of integers ``i`` such that the ``i``-th conjugacy class
of `tbl` is contained in the `p`-core of the group of `tbl`,
see [`pcore(G::GAPGroup, p::IntegerUnion)`](@ref).

# Examples
```jldoctest
julia> println(class_positions_of_pcore(character_table("2.A5"), 2))
[1, 2]
```
"""
class_positions_of_pcore(tbl::GAPGroupCharacterTable, p::IntegerUnion) = Vector{Int}(GAP.Globals.ClassPositionsOfPCore(GAPTable(tbl), GAP.Obj(p))::GapObj)

@doc raw"""
    pcore(tbl::GAPGroupCharacterTable, p::IntegerUnion)

Return the `p`-core of the group of `tbl`,
see [`pcore(G::GAPGroup, p::IntegerUnion)`](@ref),
but computed character-theoretically (see [`class_positions_of_pcore`](@ref)).

# Examples
```jldoctest
julia> order(pcore(character_table(symmetric_group(4)), 2)[1])
4
```
"""
function pcore(tbl::GAPGroupCharacterTable, p::IntegerUnion)
    @req isdefined(tbl, :group) "character table stores no group"
    t = GAPTable(tbl)
    pcorepos = GAP.Globals.ClassPositionsOfPCore(t, GAP.Obj(p))::GapObj
    P = GAP.Globals.NormalSubgroupClasses(t, pcorepos)::GapObj
    return preimages(isomorphism_to_GAP_group(tbl), P)
end

function class_positions_of_kernel(fus::Vector{Int})
  return filter(i -> fus[i] == fus[1], 1:length(fus))
end

function Base.getindex(tbl::GAPGroupCharacterTable, i::Int)
    irr = GAP.Globals.Irr(GAPTable(tbl))::GapObj
    return group_class_function(tbl, irr[i])
end
#TODO: cache the irreducibles in the table

# in order to make `tbl[end]` work
Base.lastindex(tbl::GAPGroupCharacterTable) = length(tbl)

# in order to make `findfirst` and `findall` work
function Base.keys(tbl::GAPGroupCharacterTable)
    return keys(1:length(tbl))
end

function Base.getindex(tbl::GAPGroupCharacterTable, i::Int, j::Int)
    irr = GAP.Globals.Irr(GAPTable(tbl))::GapObj
    val = irr[i, j]
    return QQAbElem(val)
end
#TODO: cache the values once they are known?

Base.iterate(tbl::GAPGroupCharacterTable, state = 1) = state > nrows(tbl) ? nothing : (tbl[state], state+1)

"""
    mod(tbl::GAPGroupCharacterTable, p::Int)

Return the `p`-modular character table of `tbl`,
or `nothing` if this table cannot be computed.

An exception is thrown if `tbl` is not an ordinary character table.
"""
function Base.mod(tbl::GAPGroupCharacterTable, p::Int)
    @req is_prime(p) "p must be a prime integer"
    characteristic(tbl) == 0 || error("tbl mod p only for ordinary table tbl")

    modtbls = get_attribute!(() -> Dict{Int,Any}(), tbl, :brauer_tables)
    if ! haskey(modtbls, p)
      modtblgap = mod(GAPTable(tbl), p)::GapObj
      if modtblgap === GAP.Globals.fail
        modtbls[p] = nothing
      elseif isdefined(tbl, :group)
        modtbls[p] = GAPGroupCharacterTable(group(tbl), modtblgap,
                       isomorphism_to_GAP_group(tbl), p)
      else
        modtbls[p] = GAPGroupCharacterTable(modtblgap, p)
      end
    end

    set_attribute!(modtbls[p], :ordinary_table, tbl)
    return modtbls[p]
end

"""
    decomposition_matrix(modtbl::GAPGroupCharacterTable)

Return the decomposition matrix (of type `ZZMatrix`) of the Brauer character
table `modtbl`.
The rows and columns are indexed by the irreducible characters of the ordinary
character table of `modtbl` and the irreducible characters of `modtbl`,
respectively,

# Examples
```jldoctest
julia> t = character_table("A5"); t2 = mod(t, 2);

julia> decomposition_matrix(t2)
[1   0   0   0]
[1   0   1   0]
[1   1   0   0]
[0   0   0   1]
[1   1   1   0]
```
"""
function decomposition_matrix(modtbl::GAPGroupCharacterTable)
    @req is_prime(characteristic(modtbl)) "characteristic of tbl must be a prime integer"
    return matrix(ZZ, GAP.Globals.DecompositionMatrix(GAPTable(modtbl))::GapObj)
end

@doc raw"""
    class_multiplication_coefficient(::Type{T} = ZZRingElem, tbl::GAPGroupCharacterTable, i::Int, j::Int, k::Int) where T <: IntegerUnion

Return the class multiplication coefficient of the classes `i`, `j`, and `k`
of the group ``G`` with ordinary character table `tbl`,
as an instance of `T`.

The class multiplication coefficient ``c_{i,j,k}`` of the classes
``i, j, k`` equals the  number of pairs ``(x, y)`` of elements ``x, y \in G``
such that ``x`` lies in class ``i``, ``y`` lies in class ``j``,
and their product ``xy`` is a fixed element of class ``k``.

In the center of the group algebra of ``G``, these numbers are found as
coefficients of the decomposition of the product of two class sums ``K_i``
and ``K_j`` into class sums:
```math
K_i K_j = \sum_k c_{ijk} K_k.
```
Given the character table of a finite group ``G``,
whose classes are ``C_1, \ldots, C_r`` with representatives ``g_i \in C_i``,
the class multiplication coefficient ``c_{ijk}`` can be computed
with the following formula:
```math
    c_{ijk} = |C_i| |C_j| / |G|
              \sum_{\chi \in Irr(G)} \chi(g_i) \chi(g_j) \chi(g_k^{-1})
              / \chi(1).
```
On the other hand the knowledge of the class multiplication coefficients
admits the computation of the irreducible characters of ``G``.

# Examples
```jldoctest
julia> class_multiplication_coefficient(character_table("A5"), 2, 3, 4)
5

julia> class_multiplication_coefficient(character_table("A5"), 2, 4, 4)
0
```
"""
function class_multiplication_coefficient(::Type{T}, tbl::GAPGroupCharacterTable, i::Int, j::Int, k::Int) where T <: IntegerUnion
  return T(GAP.Globals.ClassMultiplicationCoefficient(GAPTable(tbl), i, j, k)::GAP.Obj)
end

class_multiplication_coefficient(tbl::GAPGroupCharacterTable, i::Int, j::Int, k::Int) = class_multiplication_coefficient(ZZRingElem, tbl, i, j, k)

@doc raw"""
    possible_class_fusions(subtbl::GAPGroupCharacterTable, tbl::GAPGroupCharacterTable)

Return the array of possible class fusions from `subtbl` to `tbl`.
Each entry is an array of positive integers, where the value at position `i`
is the position of the conjugacy class in `tbl` that contains the `i`-th class
of `subtbl`.

# Examples
```jldoctest
julia> possible_class_fusions(character_table("A5"), character_table("A6"))
4-element Vector{Vector{Int64}}:
 [1, 2, 3, 6, 7]
 [1, 2, 3, 7, 6]
 [1, 2, 4, 6, 7]
 [1, 2, 4, 7, 6]
```
"""
function possible_class_fusions(subtbl::GAPGroupCharacterTable, tbl::GAPGroupCharacterTable)
  fus = GAP.Globals.PossibleClassFusions(GAPTable(subtbl), GAPTable(tbl))::GapObj
  return [Vector{Int}(x::GapObj) for x in fus]
end

#############################################################################
##
##  character parameters, class parameters
##
function _translate_parameter(para)
    if GAPWrap.IsChar(para)
      return Char(para)
    elseif GAPWrap.IsInt(para)
      return para
    elseif GAPWrap.IsCyc(para)
      # happens for the `P:Q` table, only roots of unity occur
      return [x for x in GAP.Globals.DescriptionOfRootOfUnity(para)::GapObj]
    elseif ! GAPWrap.IsList(para)
      # What can this parameter be?
      return GAP.gap_to_julia(para)
    elseif length(para) == 0
      return Int[]
    else
      return [_translate_parameter(x) for x in para]
    end
end

function _translate_parameter_list(paras)
    if all(x -> GAPWrap.IsList(x) && length(x) == 2 && x[1] == 1, paras)
      # If all parameters are lists of length 2 with first entry `1` then
      # take the second entry.
      paras = [x[2] for x in paras]
      return [_translate_parameter(x) for x in paras]
    else
      # Create tuples `(t, v)` where `t` is the parameter type
      # and `v` is the value for this type.
      return [(x[1], x[2]) for x in [_translate_parameter(x) for x in paras]]
    end
end

@doc raw"""
    character_parameters(tbl::GAPGroupCharacterTable)

Return a vector of character parameters for the rows of `tbl`
if such parameters are stored, and `nothing` otherwise.

# Examples
```jldoctest
julia> character_parameters(character_table("S5"))
7-element Vector{Vector{Int64}}:
 [5]
 [1, 1, 1, 1, 1]
 [3, 1, 1]
 [4, 1]
 [2, 1, 1, 1]
 [3, 2]
 [2, 2, 1]

julia> character_parameters(character_table("M11"))
```
"""
function character_parameters(tbl::GAPGroupCharacterTable)
    return get_attribute!(tbl, :character_parameters) do
      GAPt = GAPTable(tbl)
      GAP.Globals.HasCharacterParameters(GAPt)::Bool || return nothing
      paras = Vector{GAP.Obj}(GAP.Globals.CharacterParameters(GAPt)::GapObj)
      return _translate_parameter_list(paras)
    end
end

@doc raw"""
    class_parameters(tbl::GAPGroupCharacterTable)

Return a vector of class parameters for the columns of `tbl`
if such parameters are stored, and `nothing` otherwise.

# Examples
```jldoctest
julia> class_parameters(character_table("S5"))
7-element Vector{Vector{Int64}}:
 [1, 1, 1, 1, 1]
 [2, 2, 1]
 [3, 1, 1]
 [5]
 [2, 1, 1, 1]
 [4, 1]
 [3, 2]

julia> class_parameters(character_table("M11"))
```
"""
function class_parameters(tbl::GAPGroupCharacterTable)
    return get_attribute!(tbl, :class_parameters) do
      GAPt = GAPTable(tbl)
      GAP.Globals.HasClassParameters(GAPt)::Bool || return nothing
      paras = Vector{GAP.Obj}(GAP.Globals.ClassParameters(GAPt)::GapObj)
      return _translate_parameter_list(paras)
    end
end

@doc raw"""
    names_of_fusion_sources(tbl::GAPGroupCharacterTable)

Return the array of strings that are identifiers of those character tables
which store a class fusion to `tbl`.
"""
function names_of_fusion_sources(tbl::GAPGroupCharacterTable)
    return [string(name) for name in GAP.Globals.NamesOfFusionSources(GAPTable(tbl))::GapObj]
end

@doc raw"""
    known_class_fusion(tbl1::GAPGroupCharacterTable, tbl2::GAPGroupCharacterTable)

Return `(flag, fus)` where `flag == true` if a class fusion to `tbl2` is stored
on `tbl1`, and `flag == false` otherwise.

In the former case,
`fus` is the vector of integers, of length `ncols(tbl1)`,
such that the $i$-th conjugacy class of `tbl1` corresponds to the `fus`[$i$]-th
conjugacy class of `tbl2`, in the following sense.

If the group of `tbl1` is a *subgroup* of the group of `tbl2` then
the $i$-th conjugacy class of `tbl1` is contained in the `fus`[$i$]-th
conjugacy class of `tbl2`.
If the group of `tbl2` is a *factor group* of the group of `tbl1` then
the image of the $i$-th conjugacy class `tbl1` under the relevant epimorphism
is the `fus`[$i$]-th conjugacy class of `tbl2`.
"""
function known_class_fusion(subtbl::GAPGroupCharacterTable, tbl::GAPGroupCharacterTable)
    map = GAP.Globals.GetFusionMap(GAPTable(subtbl), GAPTable(tbl))::GapObj
    if map === GAP.Globals.fail
      return (false, Int[])
    else
      return (true, Vector{Int}(map))
    end
end


#############################################################################
##
##  class functions (and characters)
##
abstract type GroupClassFunction end

struct GAPGroupClassFunction <: GroupClassFunction
    table::GAPGroupCharacterTable
    values::GapObj
end

function Base.show(io::IO, chi::GAPGroupClassFunction)
    print(io, "group_class_function($(chi.table), $(values(chi)))")
end

function values(chi::GAPGroupClassFunction)
    gapvalues = GAP.Globals.ValuesOfClassFunction(chi.values)::GapObj
    return [QQAbElem(x) for x in gapvalues]
end

function group_class_function(tbl::GAPGroupCharacterTable, values::GapObj)
    @req GAPWrap.IsClassFunction(values) "values must be a class function"
    return GAPGroupClassFunction(tbl, values)
end

function group_class_function(tbl::GAPGroupCharacterTable, values::Vector{<:QQAbElem})
    gapvalues = GapObj([GAP.Obj(x) for x in values])
    return GAPGroupClassFunction(tbl, GAP.Globals.ClassFunction(GAPTable(tbl), gapvalues)::GapObj)
end

function group_class_function(G::GAPGroup, values::Vector{<:QQAbElem})
    return group_class_function(character_table(G), values)
end

@doc raw"""
    trivial_character(tbl::GAPGroupCharacterTable)

Return the character of `tbl` that has the value `QQAbElem(1)` in each position.
"""
function trivial_character(tbl::GAPGroupCharacterTable)
    val = QQAbElem(1)
    return group_class_function(tbl, [val for i in 1:ncols(tbl)])
end

@doc raw"""
    trivial_character(G::GAPGroup)

Return the character of (the ordinary character table of) `G`
that has the value `QQAbElem(1)` in each position.
"""
function trivial_character(G::GAPGroup)
    val = QQAbElem(1)
    return group_class_function(G, [val for i in 1:Int(number_conjugacy_classes(G))])
end

@doc raw"""
    natural_character(G::PermGroup)

Return the permutation character of degree `degree(G)`
that maps each element of `G` to the number of its fixed points.
"""
function natural_character(G::PermGroup)
    tbl = character_table(G)
    ccl = conjugacy_classes(tbl)
    FF = abelian_closure(QQ)[1]
    n = degree(G)
    vals = [FF(n - number_moved_points(representative(x))) for x in ccl]
    return group_class_function(G, vals)
end

@doc raw"""
    natural_character(G::Union{MatrixGroup{QQFieldElem}, MatrixGroup{nf_elem}})

Return the character that maps each element of `G` to its trace.
We assume that the entries of the elements of `G` are either of type `QQFieldElem`
or contained in a cyclotomic field.
"""
function natural_character(G::Union{MatrixGroup{ZZRingElem}, MatrixGroup{QQFieldElem}, MatrixGroup{nf_elem}})
    tbl = character_table(G)
    ccl = conjugacy_classes(tbl)
    FF = abelian_closure(QQ)[1]
    vals = [FF(tr(representative(x))) for x in ccl]
    return group_class_function(G, vals)
end

@doc raw"""
    induce(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable[, fusion::Vector{Int}])

Return the class function of `tbl` that is induced from `chi`,
which is a class function of a subgroup of the group of `tbl`.
The default for the class fusion `fus` is given either by the fusion of the
conjugacy classes of the two character tables (if groups are stored in the
tables) or by the class fusion given by `known_class_fusion` for the two
tables.

# Examples
```jldoctest
julia> s = character_table("A5");  t = character_table("A6");

julia> maps = possible_class_fusions(s, t);  length(maps)
4

julia> chi = trivial_character(s);

julia> ind = [induce(chi, t, x) for x in maps];

julia> length(Set(ind))
2
```
"""
function induce(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable)
  # If the class fusion is stored then use it.
  known, fus = known_class_fusion(chi.table, tbl)
  known && return induce(chi, tbl, fus)

  subtbl = chi.table
  if !(isdefined(tbl, :group) && isdefined(subtbl, :group))
    # If there is no stored group then let GAP try to find the result.
    fus = GAP.Globals.FusionConjugacyClasses(GAPTable(subtbl), GAPTable(tbl))
    @req fus !== GAP.Globals.fail "class fusion is not uniquely determinaed"
    ind = GAP.Globals.InducedClassFunctionsByFusionMap(GAPTable(chi.table),
          GAPTable(tbl), GapObj([chi.values]), GapObj(fus))::GapObj
    return GAPGroupClassFunction(tbl, ind[1])
  else
    # Dispatch on the types of the stored groups.
    return _induce(chi, tbl, group(subtbl), group(tbl))
  end
end

function induce(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable, fusion::Vector{Int})
  ind = GAP.Globals.InducedClassFunctionsByFusionMap(GAPTable(chi.table),
          GAPTable(tbl), GapObj([chi.values]), GapObj(fusion))::GapObj
  return GAPGroupClassFunction(tbl, ind[1])
end

# If `GAPGroup` groups are stored then we let GAP try to find the result.
function _induce(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable, G_chi::GAPGroup, G_tbl::GAPGroup)
  ind = GAP.Globals.InducedClassFunction(chi.values, GAPTable(tbl))::GapObj
  return GAPGroupClassFunction(tbl, ind)
end

# If `GrpAbFinGen` groups are stored then we have to work with explicit
# embeddings.
function _induce(chi::GAPGroupClassFunction,
              tbl::GAPGroupCharacterTable,
              G_chi::GrpAbFinGen, G_tbl::GrpAbFinGen)
  H, emb = is_subgroup(G_chi, G_tbl)
  Hreps = [emb(representative(C)) for C in conjugacy_classes(chi.table)]
  Greps = [representative(C) for C in conjugacy_classes(tbl)]
  fus = [findfirst(x -> x == y, Greps) for y in Hreps]
  return induce(chi, tbl, fus)
end


@doc raw"""
    induced_cyclic(tbl::GAPGroupCharacterTable)

Return the array of permutation characters of `tbl` that are induced from
cyclic subgroups.
"""
function induced_cyclic(tbl::GAPGroupCharacterTable)
    return [GAPGroupClassFunction(tbl, chi) for chi in GAP.Globals.InducedCyclic(GAPTable(tbl))::GapObj]
end

"""
    restrict(chi::GAPGroupClassFunction, subtbl::GAPGroupCharacterTable[, fusion::Vector{Int}])

Return the class function of `subtbl` that is the restriction of `chi`,
which is a class function of a supergroup of the group of `subtbl`.
The default for the class fusion `fus` is given either by the fusion of the
conjugacy classes of the two character tables (if groups are stored in the
tables) or by the class fusion given by `known_class_fusion` for the two
tables.

# Examples
```jldoctest
julia> s = character_table("A5");  t = character_table("A6");

julia> maps = possible_class_fusions(s, t);  length(maps)
4

julia> chi = t[2];  rest = [restrict(chi, s, x) for x in maps];

julia> length(Set(rest))
2
```
"""
function restrict(chi::GAPGroupClassFunction, subtbl::GAPGroupCharacterTable)
  # If the class fusion is stored then use it.
  known, fus = known_class_fusion(subtbl, chi.table)
  known && return restrict(chi, subtbl, fus)

  tbl = chi.table
  if !(isdefined(tbl, :group) && isdefined(subtbl, :group))
    # If there is no stored group then let GAP try to find the result.
    fus = GAP.Globals.FusionConjugacyClasses(GAPTable(subtbl), GAPTable(tbl))::GapObj
    @req fus !== GAP.Globals.fail "class fusion is not uniquely determinaed"
    rest = GAP.Globals.ELMS_LIST(chi.values, GapObj(fus))::GapObj
    rest = GAP.Globals.ClassFunction(GAPTable(subtbl), rest)::GapObj
    return GAPGroupClassFunction(subtbl, rest)
  else
    # Dispatch on the types of the stored groups.
    return _restrict(chi, subtbl, group(tbl), group(subtbl))
  end
end

function restrict(chi::GAPGroupClassFunction, subtbl::GAPGroupCharacterTable, fusion::Vector{Int})
  rest = GAP.Globals.ELMS_LIST(chi.values, GapObj(fusion))::GapObj
  rest = GAP.Globals.ClassFunction(GAPTable(subtbl), rest)::GapObj
  return GAPGroupClassFunction(subtbl, rest)
end

# If `GAPGroup` groups are stored then we let GAP try to find the result.
function _restrict(chi::GAPGroupClassFunction, subtbl::GAPGroupCharacterTable, G_chi::GAPGroup, G_tbl::GAPGroup)
  tbl = chi.table
  fus = GAP.Globals.FusionConjugacyClasses(GAPTable(subtbl), GAPTable(tbl))::GapObj
  rest = GAP.Globals.ELMS_LIST(chi.values, GapObj(fus))::GapObj
  rest = GAP.Globals.ClassFunction(GAPTable(subtbl), rest)::GapObj
  return GAPGroupClassFunction(subtbl, rest)
end

# If `GrpAbFinGen` groups are stored then we have to work with explicit
# embeddings.
function _restrict(chi::GAPGroupClassFunction,
              subtbl::GAPGroupCharacterTable,
              G_chi::GrpAbFinGen, G_subtbl::GrpAbFinGen)
  H, emb = is_subgroup(G_subtbl, G_chi)
  Hreps = [emb(representative(C)) for C in conjugacy_classes(subtbl)]
  Greps = [representative(C) for C in conjugacy_classes(chi.table)]
  fus = [findfirst(x -> x == y, Greps) for y in Hreps]
  return restrict(chi, subtbl, fus)
end


Base.length(chi::GAPGroupClassFunction) = length(chi.values)

Base.iterate(chi::GAPGroupClassFunction, state = 1) = state > length(chi.values) ? nothing : (chi[state], state+1)

@doc raw"""
    degree(::Type{T} = QQFieldElem, chi::GAPGroupClassFunction)
           where T <: Union{IntegerUnion, ZZRingElem, QQFieldElem, QQAbElem}

Return `chi[1]`, as an instance of `T`.
"""
Nemo.degree(chi::GAPGroupClassFunction) = Nemo.degree(QQFieldElem, chi)::QQFieldElem

Nemo.degree(::Type{QQFieldElem}, chi::GAPGroupClassFunction) = Nemo.coeff(values(chi)[1].data, 0)::QQFieldElem

Nemo.degree(::Type{ZZRingElem}, chi::GAPGroupClassFunction) = ZZ(Nemo.coeff(values(chi)[1].data, 0))::ZZRingElem

Nemo.degree(::Type{QQAbElem}, chi::GAPGroupClassFunction) = values(chi)[1]::QQAbElem{nf_elem}

Nemo.degree(::Type{T}, chi::GAPGroupClassFunction) where T <: IntegerUnion = T(Nemo.degree(ZZRingElem, chi))::T

# access character values
function Base.getindex(chi::GAPGroupClassFunction, i::Int)
  vals = GAP.Globals.ValuesOfClassFunction(chi.values)::GapObj
  return QQAbElem(vals[i])
end

# arithmetic with class functions
function Base.:(==)(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    @req chi.table === psi.table "character tables must be identical"
#T check_parent?
    return chi.values == psi.values
end

# Currently we cannot implement a `hash` method based on the values,
# since `hash(::QQAbElem)` is based on `objectid`.
function Base.hash(chi::GAPGroupClassFunction, h::UInt)
  return Base.hash(chi.table, h)
end

function Base.:+(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    @req chi.table === psi.table "character tables must be identical"
    return GAPGroupClassFunction(chi.table, chi.values + psi.values)
end

Base.:-(chi::GAPGroupClassFunction) = GAPGroupClassFunction(chi.table, - chi.values)

function Base.:-(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    @req chi.table === psi.table "character tables must be identical"
    return GAPGroupClassFunction(chi.table, chi.values - psi.values)
end

function Base.:*(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    @req chi.table === psi.table "character tables must be identical"
    return GAPGroupClassFunction(chi.table, chi.values * psi.values)
end

function Base.zero(chi::GAPGroupClassFunction)
    val = QQAbElem(0)
    return group_class_function(chi.table, [val for i in 1:length(chi)])
end

Base.one(chi::GAPGroupClassFunction) = trivial_character(chi.table)

@doc raw"""
    scalar_product(::Type{T} = QQFieldElem, chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
                   where T <: Union{IntegerUnion, ZZRingElem, QQFieldElem, QQAbElem}

Return $\sum_{g \in G}$ `chi`($g$) `conj(psi)`($g$) / $|G|$,
where $G$ is the group of both `chi` and `psi`.
The result is an instance of `T`.
"""
scalar_product(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction) = scalar_product(QQFieldElem, chi, psi)

function scalar_product(::Type{T}, chi::GAPGroupClassFunction, psi::GAPGroupClassFunction) where T <: Union{Integer, ZZRingElem, QQFieldElem, QQAbElem}
    @req chi.table === psi.table "character tables must be identical"
    return T(GAPWrap.ScalarProduct(chi.table.GAPTable, chi.values, psi.values))::T
end


@doc raw"""
    coordinates(::Type{T} = QQFieldElem, chi::GAPGroupClassFunction)
                   where T <: Union{IntegerUnion, ZZRingElem, QQFieldElem, QQAbElem}

Return the vector $[a_1, a_2, \ldots, a_n]$ of scalar products
(see [`scalar_product`](@ref)) of `chi` with the irreducible characters
$[t[1], t[2], \ldots, t[n]]$ of the character table $t$ of `chi`,
that is, `chi` is equal to $\sum_{i==1}^n a_i t[i]$.
The result is an instance of `Vector{T}`.

# Examples
```jldoctest
julia> g = symmetric_group(4)
Sym( [ 1 .. 4 ] )

julia> chi = natural_character(g);

julia> coordinates(Int, chi)
5-element Vector{Int64}:
 0
 0
 0
 1
 1

julia> t = chi.table;  t3 = mod(t, 3);  chi3 = restrict(chi, t3);

julia> coordinates(Int, chi3)
4-element Vector{Int64}:
 0
 1
 0
 1
```
"""
coordinates(chi::GAPGroupClassFunction) = coordinates(QQFieldElem, chi)

function coordinates(::Type{T}, chi::GAPGroupClassFunction) where T <: Union{Integer, ZZRingElem, QQFieldElem, QQAbElem}
    t = chi.table
    GAPt = t.GAPTable
    if characteristic(t) == 0
      # use scalar products for an ordinary character
      c = GAPWrap.MatScalarProducts(GAPt, GAPWrap.Irr(GAPt), GapObj([chi.values]))
    else
      # decompose a Brauer character
      c = GAPWrap.Decomposition(GAPWrap.Irr(GAPt), GapObj([chi.values]), GapObj("nonnegative"))
    end
    return Vector{T}(c[1])::Vector{T}
end

function Base.:*(n::IntegerUnion, chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(chi.table, n * chi.values)
end

function Base.:^(chi::GAPGroupClassFunction, n::IntegerUnion)
    return GAPGroupClassFunction(chi.table, chi.values ^ n)
end

function Base.:^(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable)
    return induce(chi, tbl)
end

function Base.:^(chi::GAPGroupClassFunction, g::Union{GAPGroupElem, GrpAbFinGenElem})
    tbl = chi.table
    @req isdefined(tbl, :group) "character table stores no group"
    G = group(tbl)
    ccl = conjugacy_classes(tbl)
    reps = [representative(c) for c in ccl]
    pi = [findfirst(x -> x^g in c, reps) for c in ccl]
    return group_class_function(tbl, values(chi)[pi])
end

@doc raw"""
    conj(chi::GAPGroupClassFunction)

Return the class function whose values are the complex conjugates of
the values of `chi`.
"""
function conj(chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(chi.table, GAP.Globals.GaloisCyc(chi.values, -1)::GAP.Obj)
end

@doc raw"""
    (sigma::QQAbAutomorphism)(chi::GAPGroupClassFunction)

Return the class function whose values are the images of the values of `chi`
under `sigma`.
"""
function (sigma::QQAbAutomorphism)(chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(chi.table, GAP.Globals.GaloisCyc(chi.values, sigma.exp)::GAP.Obj)
end

Base.:^(chi::GAPGroupClassFunction, sigma::QQAbAutomorphism) = sigma(chi)

@doc raw"""
    is_rational(chi::GAPGroupClassFunction)

Return `true` if all values of `chi` are rational, i.e., in `QQ`,
and `false` otherwise.

```jldoctest
julia> all(is_rational, character_table(symmetric_group(4)))
true

julia> all(is_rational, character_table(alternating_group(4)))
false
```
"""
function is_rational(chi::GAPGroupClassFunction)
    return all(is_rational, [x.data for x in values(chi)])
end

@doc raw"""
    is_irreducible(chi::GAPGroupClassFunction)

Return `true` if `chi` is an irreducible character, and `false` otherwise.

A character is irreducible if it cannot be written as the sum of two
characters.
For ordinary characters this can be checked using the scalar product of
class functions (see [`scalar_product`](@ref).
For Brauer characters there is no generic method for checking irreducibility.
"""
function is_irreducible(chi::GAPGroupClassFunction)
    return GAPWrap.IsIrreducibleCharacter(chi.values)
end

@doc raw"""
    is_faithful(chi::GAPGroupClassFunction)

Return `true` if the value of `chi` at the identity element does not occur
as value of `chi` at any other element, and `false` otherwise.

If `chi` is an ordinary character then `true` is returned if and only if
the representations affording `chi` have trivial kernel.

# Examples
```jldoctest
julia> show(map(is_faithful, character_table(symmetric_group(3))))
Bool[0, 1, 0]
```
"""
function is_faithful(chi::GAPGroupClassFunction)
    return length(class_positions_of_kernel(chi)) == 1
end

# Apply a class function to a group element.
function(chi::GAPGroupClassFunction)(g::Union{GAPGroupElem, GrpAbFinGenElem})
    tbl = chi.table

    # Identify the conjugacy class of `g`.
    ccl = conjugacy_classes(tbl)
    for i in 1:length(ccl)
      if g in ccl[i]
        return chi[i]
      end
    end
    error("$g is not an element in the underlying group")
end

@doc raw"""
    class_positions_of_kernel(chi::GAPGroupClassFunction)

Return the array of those integers `i` such that `chi[i] == chi[1]` holds.

# Examples
```jldoctest
julia> println(class_positions_of_kernel(character_table("2.A5")[2]))
[1, 2]
```
"""
function class_positions_of_kernel(chi::GAPGroupClassFunction)
    deg = chi[1]
    return filter(i -> chi[i] == deg, 1:length(chi))
end

function class_positions_of_kernel(list::Vector{T}) where T
    length(list) == 0 && return T[]
    deg = list[1]
    return filter(i -> list[i] == deg, 1:length(list))
end

@doc raw"""
    kernel(chi::GAPGroupClassFunction)

Return `C, f` where `C` is the kernel of `chi`
(i.e. the largest normal subgroup of the underlying group `G` of `chi`
such that `chi` maps each element of `C` to `chi[1]`)
and `f` is the embedding morphism of `C` into `G`.

# Examples
```jldoctest
julia> t = character_table(symmetric_group(4));

julia> chi = t[3];  chi[1]
2

julia> C, f = kernel(chi);  order(C)
4
```
"""
function kernel(chi::GAPGroupClassFunction)
    tbl = chi.table
    @req isdefined(tbl, :group) "character table stores no group"
    GAP_K = GAP.Globals.KernelOfCharacter(GAPTable(tbl), chi.values)::GapObj
    return preimages(isomorphism_to_GAP_group(tbl), GAP_K)
end

@doc raw"""
    class_positions_of_center(chi::GAPGroupClassFunction)

Return the array of those integers `i` such that `chi[i]` is `chi[1]` times
a root of unity.

# Examples
```jldoctest
julia> println(class_positions_of_center(character_table("2.A5")[2]))
[1, 2]
```
"""
function class_positions_of_center(chi::GAPGroupClassFunction)
    return Vector{Int}(GAP.Globals.ClassPositionsOfCentre(chi.values))
end

@doc raw"""
    center(chi::GAPGroupClassFunction)

Return `C, f` where `C` is the center of `chi`
(i.e. the largest normal subgroup of the underlying group `G` of `chi`
such that `chi` maps each element of `C` to `chi[1]` times a root of unity)
and `f` is the embedding morphism of `C` into `G`.

# Examples
```jldoctest
julia> t = character_table(symmetric_group(4));

julia> chi = t[3];  chi[1]
2

julia> C, f = center(chi);  order(C)
4
```
"""
function center(chi::GAPGroupClassFunction)
    tbl = chi.table
    @req isdefined(tbl, :group) "character table stores no group"
    G = group(tbl)
    C = GAP.Globals.CentreOfCharacter(GAPTable(tbl), chi.values)::GapObj
    return preimages(isomorphism_to_GAP_group(tbl), C)
end

@doc raw"""
    det(chi::GAPGroupClassFunction)

Return the determinant character of the character `chi`.
This is defined to be the character obtained by taking the determinant
of representing matrices of any representation affording `chi`.

# Examples
```jldoctest
julia> t = character_table(symmetric_group(4));

julia> all(chi -> det(chi) == exterior_power(chi, Int(degree(chi))), t)
true
```
"""
function det(chi::GAPGroupClassFunction)
    values = GAP.Globals.DeterminantOfCharacter(chi.values)
    return GAPGroupClassFunction(chi.table, values)
end

# Note that defining the determinantal order for arbitrary characters
# is not possible in GAP.
@doc raw"""
    order(::Type{T} = ZZRingElem, chi::GAPGroupClassFunction)
          where T <: IntegerUnion

Return the determinantal order of the character `chi`.
This is defined to be the multiplicative order of `det(chi)`.

# Examples
```jldoctest
julia> println([order(chi) for chi in character_table(symmetric_group(4))])
ZZRingElem[2, 1, 2, 2, 1]
```
"""
order(chi::GAPGroupClassFunction) = order(ZZRingElem, chi)::ZZRingElem

function order(::Type{T}, chi::GAPGroupClassFunction) where T <: IntegerUnion
    return T(GAP.Globals.Order(det(chi).values))::T
end

@doc raw"""
    indicator(chi::GAPGroupClassFunction, n::Int = 2)

Return the `n`-th Frobenius-Schur indicator of `chi`, that is,
the value $(∑_{g ∈ G} chi(g^n))/|G|$, where $G$ is the group of `chi`.

If `chi` is irreducible then `indicator(chi)` is
`0` if `chi` is not real-valued,
`1` if `chi` is afforded by a real representation of $G$, and
`-1` if `chi` is real-valued but not afforded by a real representation of $G$.
"""
function indicator(chi::GAPGroupClassFunction, n::Int = 2)
    ind = GAP.Globals.Indicator(GAPTable(chi.table), GapObj([chi.values]), n)::GapObj
    return ind[1]::Int
end

@doc raw"""
    character_field(chi::GAPGroupClassFunction)

Return the pair `(F, phi)` where `F` is a number field that is generated
by the character values of `chi`, and `phi` is the embedding of `F` into
`abelian_closure(QQ)`.
"""
function character_field(chi::GAPGroupClassFunction)
    values = chi.values  # a list of GAP cyclotomics
    gapfield = GAP.Globals.Field(values)::GapObj
    N = GAPWrap.Conductor(gapfield)
    FF, = abelian_closure(QQ)
    if GAPWrap.IsCyclotomicField(gapfield)
      # In this case, the want to return a field that knows to be cyclotomic
      # (and the embedding is easy).
      F, z = AbelianClosure.cyclotomic_field(FF, N)
      f = x::nf_elem -> QQAbElem(x, N)
      finv = function(x::QQAbElem)
        g = gcd(x.c, N)
        K, = AbelianClosure.cyclotomic_field(FF, g)
        x = Hecke.force_coerce_cyclo(K, x.data)
        x = Hecke.force_coerce_cyclo(F, x)
        return x
      end
    else
      # In the general case, we have to work for the embedding.
      gapgens = GAP.Globals.GeneratorsOfField(gapfield)::GapObj
      @assert length(gapgens) == 1
      gappol = GAP.Globals.MinimalPolynomial(GAP.Globals.Rationals, gapgens[1])::GapObj
      gapcoeffs = GAP.Globals.CoefficientsOfUnivariatePolynomial(gappol)::GapObj
      v = Vector{QQFieldElem}(gapcoeffs)
      R, = polynomial_ring(QQ, "x")
      f = R(v)
      F, z = number_field(f, "z"; cached = true, check = false)
      K, zz = AbelianClosure.cyclotomic_field(FF, N)

      nfelm = QQAbElem(gapgens[1]).data

      # Compute the expression of powers of `z` as sums of roots of unity (once).
      powers = [coefficients(Hecke.force_coerce_cyclo(K, nfelm^i)) for i in 0:length(v)-2]
      c = transpose(matrix(QQ, powers))

      f = function(x::nf_elem)
        return QQAbElem(evaluate(R(x), nfelm), N)
      end

      finv = function(x::QQAbElem)
        # Write `x` w.r.t. the N-th cyclotomic field ...
        g = gcd(x.c, N)
        Kg, = AbelianClosure.cyclotomic_field(FF, g)
        x = Hecke.force_coerce_cyclo(Kg, x.data)
        x = Hecke.force_coerce_cyclo(K, x)

        # ... and then w.r.t. `F`
        a = coefficients(x)
        b = transpose(solve(c, matrix(QQ,length(a),1,a)))
        b = [b[i] for i in 1:length(b)]
        return F(b)
      end
    end

    return F, MapFromFunc(f, finv, F, FF)
end

@doc raw"""
    schur_index(chi::GAPGroupClassFunction)

Return the minimal integer `m` such that the character `m * chi`
is afforded by a representation over the character field of `chi`,
or throw an exception if the currently used character theoretic criteria
do not suffice for computing `m`.
"""
function schur_index(chi::GAPGroupClassFunction, recurse::Bool = true)
    deg = numerator(degree(chi))
    deg == 1 && return 1
    indicator(chi) == -1 && return 2

    # The character field contains an `m`-th root of unity.
    values = chi.values
    if conj(chi) == chi
      bound = ZZRingElem(2)
    else
      # Compute the conductor of the largest cyclotomic field
      # that is contained in the character field of `chi`.
      gapfield = GAP.Globals.Field(values)::GapObj
      N = GAPWrap.Conductor(gapfield)
      for n in reverse(sort(divisors(N)))
        if GAP.Globals.E(n)::GAP.Obj in gapfield
          if Base.isodd(n)
            bound = ZZRingElem(2*n)
          else
            bound = ZZRingElem(n)
          end
          break
        end
      end
    end

    # `m` divides `deg`
    bound = gcd(bound, deg)
    bound == 1 && return 1

    # `m` divides the multiplicity of `chi` in any rational character
    # with trivial Schur index.
    # - Consider permutation characters induced from cyclic subgroups.
    tbl = chi.table
    for psi in induced_cyclic(tbl)
      bound = gcd(bound, scalar_product(ZZRingElem, chi, psi))
      bound == 1 && return 1
    end
    # - Consider characters induced from other known subgroups.
    for name in names_of_fusion_sources(tbl)
      s = character_table(name)
      if s !== nothing
        known, fus = known_class_fusion(s, tbl)
        @assert known "the class fusion is not stored"
        if length(class_positions_of_kernel(fus)) == 1
          psi = trivial_character(s)^(tbl)
          bound = gcd(bound, scalar_product(ZZRingElem, chi, psi))
          bound == 1 && return 1
        end
      end
    end

    if recurse
      # Consider tensor products of rational characters with Schur index 1.
      cand = filter(psi -> degree(character_field(psi)[1]) == 1 &&
                           schur_index(psi, false) == 1, collect(tbl))
      for i in 1:length(cand)
        for j in 1:i
          bound = gcd(bound, scalar_product(ZZRingElem, chi, cand[i] * cand[j]))
          bound == 1 && return 1
        end
      end
    end

    if isdefined(tbl, :group) && hasproperty(GAP.Globals, :SchurIndexByCharacter)
      g = group(tbl)
      return GAP.Globals.SchurIndexByCharacter(GAP.Globals.Rationals, codomain(isomorphism_to_GAP_group(tbl)), values)
    end

    # For the moment, we do not have more character theoretic criteria.
    error("cannot determine the Schur index with the currently used criteria")
end

@doc raw"""
    symmetrizations(characters::Vector{GAPGroupClassFunction}, n::Int)

Return the vector of symmetrizations of `characters` with the ordinary
irreducible characters of the symmetric group of degree `n`.

The symmetrization $\chi^{[\lambda]}$ of the character $\chi$
with the character $\lambda$ of the symmetric group $S_n$ of degree `n`
is defined by

```math
\chi^{[\lambda]}(g) =
(\sum_{\rho\in S_n} \lambda(\rho) \prod_{k=1}^n \chi(g^k)^{a_k(\rho)} ) / n!,
```

where $a_k(\rho)$ is the number of cycles of length $k$ in $\rho$.

Note that the returned list may contain zero class functions,
and duplicates are not deleted.

For special kinds of symmetrizations, see [`symmetric_parts`](@ref),
[`anti_symmetric_parts`](@ref), [`orthogonal_components`](@ref),
[`symplectic_components`](@ref), [`exterior_power`](@ref),
[`symmetric_power`](@ref).
"""
function symmetrizations(characters::Vector{GAPGroupClassFunction}, n::Int)
    length(characters) == 0 && return eltype(typeof(characters))[]
    tbl = characters[1].table
    return [group_class_function(tbl, chi)
            for chi in GAP.Globals.Symmetrizations(GAPTable(tbl),
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

@doc raw"""
    symmetric_parts(characters::Vector{GAPGroupClassFunction}, n::Int)

Return the vector of symmetrizations of `characters`
with the trivial character of the symmetric group of degree `n`,
see [`symmetrizations`](@ref).
"""
function symmetric_parts(characters::Vector{GAPGroupClassFunction}, n::Int)
    length(characters) == 0 && return eltype(typeof(characters))[]
    tbl = characters[1].table
    return [group_class_function(tbl, chi)
            for chi in GAP.Globals.SymmetricParts(GAPTable(tbl),
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

@doc raw"""
    anti_symmetric_parts(characters::Vector{GAPGroupClassFunction}, n::Int)

Return the vector of symmetrizations of `characters`
with the sign character of the symmetric group of degree `n`,
see [`symmetrizations`](@ref).
"""
function anti_symmetric_parts(characters::Vector{GAPGroupClassFunction}, n::Int)
    length(characters) == 0 && return eltype(typeof(characters))[]
    tbl = characters[1].table
    return [group_class_function(tbl, chi)
            for chi in GAP.Globals.AntiSymmetricParts(GAPTable(tbl),
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

@doc raw"""
    exterior_power(chi::GAPGroupClassFunction, n::Int)

Return the class function of the `n`-th exterior power of the module that is
afforded by `chi`.

This exterior power is the symmetrization of `chi` with the sign character
of the symmetric group of degree `n`,
see also [`symmetrizations`](@ref) and [`anti_symmetric_parts`](@ref).
"""
function exterior_power(chi::GAPGroupClassFunction, n::Int)
#T when GAP's `ExteriorPower` method becomes available then use it
    tbl = chi.table
    return group_class_function(tbl,
      GAP.Globals.AntiSymmetricParts(GAPTable(tbl), GAP.Obj([chi.values]), n)[1])
end

@doc raw"""
    symmetric_power(chi::GAPGroupClassFunction, n::Int)

Return the class function of the `n`-th symmetric power of the module that is
afforded by `chi`.

This symmetric power is the symmetrization of `chi` with the trivial character
of the symmetric group of degree `n`,
see also [`symmetrizations`](@ref) and [`symmetric_parts`](@ref).
"""
function symmetric_power(chi::GAPGroupClassFunction, n::Int)
#T when GAP's `SymmetricPower` method becomes available then use it
    tbl = chi.table
    return group_class_function(tbl,
      GAP.Globals.SymmetricParts(GAPTable(tbl), GAP.Obj([chi.values]), n)[1])
end

@doc raw"""
    orthogonal_components(characters::Vector{GAPGroupClassFunction}, n::Int)

Return the vector of the so-called Murnaghan components of the $m$-th
tensor powers of the entries of `characters`, for $m$ up to `n`,
where `n` must be at least 2 and at most 6
and where we assume that the entries of `characters` are irreducible
characters with Frobenius-Schur indicator +1, see [`indicator`](@ref).
"""
function orthogonal_components(characters::Vector{GAPGroupClassFunction}, n::Int)
    @req (2 <= n && n <= 6) "the second ergument must be in 2:6"
    length(characters) == 0 && return eltype(typeof(characters))[]
    tbl = characters[1].table
    return [group_class_function(tbl, chi)
            for chi in GAP.Globals.OrthogonalComponents(GAPTable(tbl),
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

@doc raw"""
    symplectic_components(characters::Vector{GAPGroupClassFunction}, n::Int)

Return the vector of the Murnaghan components of the $m$-th tensor powers
of the entries of `characters`, for $m$ up to `n`,
where `n` must be at least 2 and at most 6
and where we assume that the entries of `characters` are irreducible
characters with Frobenius-Schur indicator -1, see [`indicator`](@ref).
"""
function symplectic_components(characters::Vector{GAPGroupClassFunction}, n::Int)
    @req (2 <= n && n <= 6) "the second ergument must be in 2:6"
    length(characters) == 0 && return eltype(typeof(characters))[]
    tbl = characters[1].table
    return [group_class_function(tbl, chi)
            for chi in GAP.Globals.SymplecticComponents(GAPTable(tbl),
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

function character_table_complex_reflection_group(m::Int, p::Int, n::Int)
    @req p == 1 "the case G(m,p,n) with p != 1 is not (yet) supported"
    tbl = GAP.Globals.CharacterTableWreathSymmetric(
            GAP.Globals.CharacterTable(GapObj("Cyclic"), m)::GapObj, n):GapObj
    tbl = GAPGroupCharacterTable(tbl, 0)
    set_attribute!(tbl, :type, (m, p, n))

    return tbl
end
