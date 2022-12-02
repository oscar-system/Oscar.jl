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

export
    all_character_table_names,
    anti_symmetric_parts,
    atlas_irrationality,
    character_field,
    character_parameters,
    character_table,
    class_lengths,
    class_multiplication_coefficient,
    class_parameters,
    class_positions_of_kernel,
    class_positions_of_pcore,
    decomposition_matrix,
    exterior_power,
    identifier,
    indicator,
    induced_class_function,
    induced_cyclic,
    is_duplicate_table,
    known_class_fusion,
    maxes,
    names_of_fusion_sources,
    natural_character,
    orders_centralizers,
    orders_class_representatives,
    orthogonal_components,
    possible_class_fusions,
    scalar_product,
    schur_index,
    symmetric_parts,
    symmetric_power,
    symmetrizations,
    symplectic_components,
    trivial_character


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

An object of type `GAPGroup` can (but need not) be stored
in the field `GAPGroup`.

The value of the field `characteristic` determines whether the table
is an ordinary one (value `0`) or a `p`-modular one (value `p`).

Objects of type `GAPGroupCharacterTable` support [`get_attribute`](@ref),
for example in order to store the already computed `p`-modular tables
in an ordinary table, and to store the corresponding ordinary table
in a `p`-modular table.
"""
@attributes mutable struct GAPGroupCharacterTable <: GroupCharacterTable
    GAPTable::GapObj  # the character table object
    characteristic::Int
    GAPGroup::GAPGroup    # the underlying group, if any

    function GAPGroupCharacterTable(G::GAPGroup, tab::GapObj, char::Int)
      return new(tab, char, G)
    end

    function GAPGroupCharacterTable(tab::GapObj, char::Int)
      #GAPGroup is left undefined
      return new(tab, char)
    end
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
function character_table(G::GAPGroup, p::Int = 0)
    tbls = get_attribute!(() -> Dict{Int,Any}(), G, :character_tables)
    return get!(tbls, p) do
      gaptbl = GAP.Globals.CharacterTable(G.X)::GapObj
      if p != 0
        # Create the `p`-modular table if possible.
        is_prime(p) || error("p must be 0 or a prime integer")
        gaptbl = GAP.Globals.mod(gaptbl, GAP.Obj(p))::GapObj
        gaptbl === GAP.Globals.fail && return nothing
      end
      return GAPGroupCharacterTable(G, gaptbl, p)
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
    hasproperty(GAP.Globals, :CTblLib) || error("no character table library available")

    if p == 0
      modid = id
    else
      is_prime(p) || error("p must be 0 or a prime integer")
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
    hasproperty(GAP.Globals, :CTblLib) || error("no character table library available")
    args = GAP.Obj([string(series), parameter], recursive = true)
    tbl = GAP.Globals.CallFuncList(GAP.Globals.CharacterTable, args)::GapObj
    tbl === GAP.Globals.fail && return nothing
    tbl = GAPGroupCharacterTable(tbl, 0)
    set_attribute!(tbl, :series, (series, parameter))
    return tbl
end


# For character tables with stored group, we take the hash value of the group.
# For character tables without stored group, we take the table identifier.
function Base.hash(tbl::GAPGroupCharacterTable, h::UInt)
  if isdefined(tbl, :GAPGroup)
    return Base.hash(tbl.GAPGroup, h)
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
@gapattribute is_duplicate_table(tbl::GAPGroupCharacterTable) = GAP.Globals.IsDuplicateTable(G.GAPTable)::Bool

"""
    all_character_table_names(L...; ordered_by = nothing)

Return an array of strings that contains all those names of character tables
in the character table library that satisfy the conditions in the array `L`.

# Examples
```jldoctest
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


@doc Markdown.doc"""
    as_sum_of_roots(val::nf_elem, root::String)

Return a string representing the element `val` of a cyclotomic field
as a sum of multiples of powers of the primitive root which is printed as
`root`.
"""
function as_sum_of_roots(val::nf_elem, root::String)
    F = parent(val)
    flag, N = Hecke.is_cyclotomic_type(F)
    flag || error("$val is not an element of a cyclotomic field")

    # `string` yields an expression of the right structure,
    # but not in terms of `root`
    # and without curly brackets for subscripts and superscripts.
    str = string(val)
    str = replace(str, "*" => "")
    str = replace(str, string(F.S) => "$(root)_{$N}")
    str = replace(str, r"\^([0-9]*)" => s"^{\1}")
    return str
end


@doc Markdown.doc"""
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
              info = Oscar.AbelianClosure.quadratic_irrationality_info(val)
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
    gaptbl = tbl.GAPTable
    size = fmpz(GAPWrap.Size(gaptbl))
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
    cents = Vector{fmpz}(GAP.Globals.SizesCentralizers(gaptbl)::GapObj)
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
      p = tbl.characteristic
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

    if isdefined(tbl, :GAPGroup)
      headerstring = string(tbl.GAPGroup)
      if tbl.characteristic != 0
        headerstring = "$headerstring mod $(tbl.characteristic)"
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
    gaptbl = tbl.GAPTable
    if isdefined(tbl, :GAPGroup)
      id = string(tbl.GAPGroup)
      if tbl.characteristic != 0
        id = "$(id)mod$(tbl.characteristic)"
      end
    else
      id = "\"" * identifier(tbl) * "\""
    end
    print(io, "character_table($id)")
end


##############################################################################
#
length(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(tbl.GAPTable)::Int
Oscar.nrows(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(tbl.GAPTable)::Int
Oscar.ncols(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(tbl.GAPTable)::Int
number_conjugacy_classes(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(tbl.GAPTable)::Int

@doc Markdown.doc"""
    order(::Type{T} = fmpz, tbl::GAPGroupCharacterTable) where T <: IntegerUnion

Return the order of the group for which `tbl` is the character table.

# Examples
```jldoctest
julia> order(character_table(symmetric_group(4)))
24

```
"""
order(tbl::GAPGroupCharacterTable) = order(fmpz, tbl)

function order(::Type{T}, tbl::GAPGroupCharacterTable) where T <: IntegerUnion
  return T(GAPWrap.Size(tbl.GAPTable))
end

@doc Markdown.doc"""
    orders_class_representatives(tbl::GAPGroupCharacterTable)

Return the array of the orders of conjugacy class representatives for `tbl`,
ordered according to the columns of `tbl`.

# Examples
```jldoctest
julia> println(orders_class_representatives(character_table("A5")))
[1, 2, 3, 5, 5]

```
"""
@gapattribute orders_class_representatives(tbl::GAPGroupCharacterTable) = Vector{Int}(GAP.Globals.OrdersClassRepresentatives(tbl.GAPTable)::GapObj)

@doc Markdown.doc"""
    orders_centralizers(tbl::GAPGroupCharacterTable)

Return the array of the orders of centralizers of conjugacy class
representatives for `tbl` in the group of `tbl`,
ordered according to the columns of `tbl`.

# Examples
```jldoctest
julia> println(orders_centralizers(character_table("A5")))
fmpz[60, 4, 3, 5, 5]

```
"""
@gapattribute orders_centralizers(tbl::GAPGroupCharacterTable) = Vector{fmpz}(GAP.Globals.SizesCentralizers(tbl.GAPTable)::GAP.Obj)

@doc Markdown.doc"""
    class_lengths(tbl::GAPGroupCharacterTable)

# Examples
```jldoctest
julia> println(class_lengths(character_table("A5")))
fmpz[1, 15, 20, 12, 12]

```
"""
@gapattribute class_lengths(tbl::GAPGroupCharacterTable) = Vector{fmpz}(GAP.Globals.SizesConjugacyClasses(tbl.GAPTable)::GapObj)

@doc Markdown.doc"""
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

julia> maxes(character_table("M")) == nothing  # not (yet) known
true

```
"""
function maxes(tbl::GAPGroupCharacterTable)
  if GAP.Globals.HasMaxes(tbl.GAPTable)::Bool
    return Vector{String}(GAP.Globals.Maxes(tbl.GAPTable)::GapObj)
  end
  return nothing
end

@doc Markdown.doc"""
    identifier(tbl::GAPGroupCharacterTable)

Return a string that identifies `tbl`.
It is used mainly for library tables.

# Examples
```jldoctest
julia> identifier(character_table("A5"))
"A5"

```
"""
@gapattribute identifier(tbl::GAPGroupCharacterTable) = string(GAP.Globals.Identifier(tbl.GAPTable)::GapObj)

@doc Markdown.doc"""
    class_positions_of_pcore(tbl::GAPGroupCharacterTable, p::IntegerUnion)

Return the array of integers ``i`` such that the ``i``-th conjugacy class
of `tbl` is contained in the `p`-core of the group of `tbl`,
see [`pcore`](@ref).

# Examples
```jldoctest
julia> println(Oscar.class_positions_of_pcore(character_table("2.A5"), 2))
[1, 2]

```
"""
class_positions_of_pcore(tbl::GAPGroupCharacterTable, p::IntegerUnion) = Vector{Int}(GAP.Globals.ClassPositionsOfPCore(tbl.GAPTable, GAP.Obj(p))::GapObj)

function class_positions_of_kernel(fus::Vector{Int})
  return filter(i -> fus[i] == fus[1], 1:length(fus))
end

function Base.getindex(tbl::GAPGroupCharacterTable, i::Int)
    irr = GAP.Globals.Irr(tbl.GAPTable)::GapObj
    return group_class_function(tbl, irr[i])
end
#TODO: cache the irreducibles in the table

function Base.getindex(tbl::GAPGroupCharacterTable, i::Int, j::Int)
    irr = GAP.Globals.Irr(tbl.GAPTable)::GapObj
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
    is_prime(p) || error("p must be a prime integer")
    tbl.characteristic == 0 || error("tbl mod p only for ordinary table tbl")

    modtbls = get_attribute!(() -> Dict{Int,Any}(), tbl, :brauer_tables)
    if ! haskey(modtbls, p)
      modtblgap = mod(tbl.GAPTable, p)::GapObj
      if modtblgap === GAP.Globals.fail
        modtbls[p] = nothing
      elseif isdefined(tbl, :GAPGroup)
        modtbls[p] = GAPGroupCharacterTable(tbl.GAPGroup, modtblgap, p)
      else
        modtbls[p] = GAPGroupCharacterTable(modtblgap, p)
      end
    end

    set_attribute!(modtbls[p], :ordinary_table, tbl)
    return modtbls[p]
end

"""
    decomposition_matrix(modtbl::GAPGroupCharacterTable)

Return the decomposition matrix (of type `fmpz_mat`) of the Brauer character
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
    is_prime(modtbl.characteristic) || error("characteristic of tbl must be a prime integer")
    return matrix(ZZ, GAP.Globals.DecompositionMatrix(modtbl.GAPTable)::GapObj)
end

@doc Markdown.doc"""
    class_multiplication_coefficient(::Type{T} = fmpz, tbl::GAPGroupCharacterTable, i::Int, j::Int, k::Int) where T <: IntegerUnion

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
  return T(GAP.Globals.ClassMultiplicationCoefficient(tbl.GAPTable, i, j, k)::GAP.Obj)
end

class_multiplication_coefficient(tbl::GAPGroupCharacterTable, i::Int, j::Int, k::Int) = class_multiplication_coefficient(fmpz, tbl, i, j, k)

@doc Markdown.doc"""
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
  fus = GAP.Globals.PossibleClassFusions(subtbl.GAPTable, tbl.GAPTable)::GapObj
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

@doc Markdown.doc"""
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
      GAPt = tbl.GAPTable
      GAP.Globals.HasCharacterParameters(GAPt)::Bool || return nothing
      paras = Vector{GAP.Obj}(GAP.Globals.CharacterParameters(GAPt)::GapObj)
      return _translate_parameter_list(paras)
    end
end

@doc Markdown.doc"""
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
      GAPt = tbl.GAPTable
      GAP.Globals.HasClassParameters(GAPt)::Bool || return nothing
      paras = Vector{GAP.Obj}(GAP.Globals.ClassParameters(GAPt)::GapObj)
      return _translate_parameter_list(paras)
    end
end

@doc Markdown.doc"""
    names_of_fusion_sources(tbl::GAPGroupCharacterTable)

Return the array of strings that are identifiers of those character tables
which store a class fusion to `tbl`.
"""
function names_of_fusion_sources(tbl::GAPGroupCharacterTable)
    return [string(name) for name in GAP.Globals.NamesOfFusionSources(tbl.GAPTable)::GapObj]
end

@doc Markdown.doc"""
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
    map = GAP.Globals.GetFusionMap(subtbl.GAPTable, tbl.GAPTable)::GapObj
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
    GAPWrap.IsClassFunction(values) || error("values must be a class function")
    return GAPGroupClassFunction(tbl, values)
end

function group_class_function(tbl::GAPGroupCharacterTable, values::Vector{<:QQAbElem})
    gapvalues = GapObj([GAP.Obj(x) for x in values])
    return GAPGroupClassFunction(tbl, GAP.Globals.ClassFunction(tbl.GAPTable, gapvalues)::GapObj)
end

function group_class_function(G::GAPGroup, values::Vector{<:QQAbElem})
    return group_class_function(character_table(G), values)
end

@doc Markdown.doc"""
    trivial_character(tbl::GAPGroupCharacterTable)

Return the character of `tbl` that has the value `QQAbElem(1)` in each position.
"""
function trivial_character(tbl::GAPGroupCharacterTable)
    val = QQAbElem(1)
    return group_class_function(tbl, [val for i in 1:ncols(tbl)])
end

@doc Markdown.doc"""
    trivial_character(G::GAPGroup)

Return the character of (the ordinary character table of) `G`
that has the value `QQAbElem(1)` in each position.
"""
function trivial_character(G::GAPGroup)
    val = QQAbElem(1)
    return group_class_function(G, [val for i in 1:Int(number_conjugacy_classes(G))])
end

@doc Markdown.doc"""
    natural_character(G::PermGroup)

Return the permutation character of degree `degree(G)`
that maps each element of `G` to the number of its fixed points.
"""
function natural_character(G::PermGroup)
    ccl = conjugacy_classes(G)
    FF = abelian_closure(QQ)[1]
    n = degree(G)
    vals = [FF(n - number_moved_points(representative(x))) for x in ccl]
    return group_class_function(G, vals)
end

@doc Markdown.doc"""
    natural_character(G::Union{MatrixGroup{fmpq}, MatrixGroup{nf_elem}})

Return the character that maps each element of `G` to its trace.
We assume that the entries of the elements of `G` are either of type `fmpq`
or contained in a cyclotomic field.
"""
function natural_character(G::Union{MatrixGroup{fmpq}, MatrixGroup{nf_elem}})
    ccl = conjugacy_classes(G)
    FF = abelian_closure(QQ)[1]
    vals = [FF(tr(representative(x))) for x in ccl]
    return group_class_function(G, vals)
end

@doc Markdown.doc"""
    induced_class_function(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable[, fusion::Vector{Int}])

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

julia> ind = [induced_class_function(chi, t, x) for x in maps];  length(ind)
4
julia> length(Set(ind))
2

```
"""
function induced_class_function(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable)
  ind = GAP.Globals.InducedClassFunction(chi.values, tbl.GAPTable)::GapObj
  return GAPGroupClassFunction(tbl, ind)
end

function induced_class_function(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable, fusion::Vector{Int})
  ind = GAP.Globals.InducedClassFunctionsByFusionMap(chi.table.GAPTable,
          tbl.GAPTable, GapObj([chi.values]), GapObj(fusion))::GapObj
  return GAPGroupClassFunction(tbl, ind[1])
end

@doc Markdown.doc"""
    induced_cyclic(tbl::GAPGroupCharacterTable)

Return the array of permutation characters of `tbl` that are induced from
cyclic subgroups.
"""
function induced_cyclic(tbl::GAPGroupCharacterTable)
    return [GAPGroupClassFunction(tbl, chi) for chi in GAP.Globals.InducedCyclic(tbl.GAPTable)::GapObj]
end

Base.length(chi::GAPGroupClassFunction) = length(chi.values)

Base.iterate(chi::GAPGroupClassFunction, state = 1) = state > length(chi.values) ? nothing : (chi[state], state+1)

@doc Markdown.doc"""
    degree(::Type{T} = fmpq, chi::GAPGroupClassFunction)
           where T <: Union{IntegerUnion, fmpz, mpq, QQAbElem}

Return `chi[1]`, as an instance of `T`.
"""
Nemo.degree(chi::GAPGroupClassFunction) = Nemo.degree(fmpq, chi)::fmpq

Nemo.degree(::Type{fmpq}, chi::GAPGroupClassFunction) = Nemo.coeff(values(chi)[1].data, 0)::fmpq

Nemo.degree(::Type{fmpz}, chi::GAPGroupClassFunction) = ZZ(Nemo.coeff(values(chi)[1].data, 0))::fmpz

Nemo.degree(::Type{QQAbElem}, chi::GAPGroupClassFunction) = values(chi)[1]::QQAbElem{nf_elem}

Nemo.degree(::Type{T}, chi::GAPGroupClassFunction) where T <: IntegerUnion = T(Nemo.degree(fmpz, chi))::T

# access character values
function Base.getindex(chi::GAPGroupClassFunction, i::Int)
  vals = GAP.Globals.ValuesOfClassFunction(chi.values)::GapObj
  return QQAbElem(vals[i])
end

# arithmetics with class functions
function Base.:(==)(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    chi.table === psi.table || error("character tables must be identical")
#T check_parent?
    return chi.values == psi.values
end

# Currently we cannot implement a `hash` method based on the values,
# since `hash(::QQAbElem)` is based on `objectid`.
function Base.hash(chi::GAPGroupClassFunction, h::UInt)
  return Base.hash(chi.table, h)
end

function Base.:+(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    chi.table === psi.table || error("character tables must be identical")
    return GAPGroupClassFunction(chi.table, chi.values + psi.values)
end

Base.:-(chi::GAPGroupClassFunction) = GAPGroupClassFunction(chi.table, - chi.values)

function Base.:-(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    chi.table === psi.table || error("character tables must be identical")
    return GAPGroupClassFunction(chi.table, chi.values - psi.values)
end

function Base.:*(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    chi.table === psi.table || error("character tables must be identical")
    return GAPGroupClassFunction(chi.table, chi.values * psi.values)
end

function Base.zero(chi::GAPGroupClassFunction)
    val = QQAbElem(0)
    return group_class_function(chi.table, [val for i in 1:length(chi)])
end

Base.one(chi::GAPGroupClassFunction) = trivial_character(chi.table)

@doc Markdown.doc"""
    scalar_product(::Type{T} = fmpq, chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
                   where T <: Union{IntegerUnion, fmpz, fmpq, QQAbElem}

Return $\sum_{g \in G}$ `chi`($g$) `conj(psi)`($g$) / $|G|$,
where $G$ is the group of both `chi` and `psi`.
"""
scalar_product(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction) = scalar_product(fmpq, chi, psi)

function scalar_product(::Type{T}, chi::GAPGroupClassFunction, psi::GAPGroupClassFunction) where T <: Union{Integer, fmpz, fmpq, QQAbElem}
    chi.table === psi.table || error("character tables must be identical")
    return T(GAP.Globals.ScalarProduct(chi.values, psi.values)::GAP.Obj)::T
end

function Base.:*(n::IntegerUnion, chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(chi.table, n * chi.values)
end

function Base.:^(chi::GAPGroupClassFunction, n::IntegerUnion)
    return GAPGroupClassFunction(chi.table, chi.values ^ n)
end

function Base.:^(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable)
    return induced_class_function(chi, tbl)
end

function Base.:^(chi::GAPGroupClassFunction, g::GAPGroupElem)
    tbl = chi.table
    isdefined(tbl, :GAPGroup) || error("character table stores no group")
    G = tbl.GAPGroup
    ccl = conjugacy_classes(G)
    if tbl.characteristic != 0
      known, fus = known_class_fusion(tbl, get_attribute(tbl, :ordinary_table))
      @assert known "the class fusion is not stored"
      ccl = ccl[fus]
    end
    reps = [representative(c) for c in ccl]
    pi = [findfirst(x -> x^g in c, reps) for c in ccl]
    return group_class_function(tbl, values(chi)[pi])
end

@doc Markdown.doc"""
    conj(chi::GAPGroupClassFunction)

Return the class function whose values are the complex conjugates of
the values of `chi`.
"""
function conj(chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(chi.table, GAP.Globals.GaloisCyc(chi.values, -1)::GAP.Obj)
end

@doc Markdown.doc"""
    (sigma::QQAbAutomorphism)(chi::GAPGroupClassFunction)

Return the class function whose values are the images of the values of `chi`
under `sigma`.
"""
function (sigma::QQAbAutomorphism)(chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(chi.table, GAP.Globals.GaloisCyc(chi.values, sigma.exp)::GAP.Obj)
end

Base.:^(chi::Oscar.GAPGroupClassFunction, sigma::QQAbAutomorphism) = sigma(chi)

@doc Markdown.doc"""
    is_irreducible(chi::GAPGroupClassFunction)

Return `true` if `chi` is an irreducible character, and `false` otherwise.

A character is irreducible if it cannot be written as the sum of two
characters.
For ordinary characters this can be checked using the scalar product of
class functions (see [`scalar_product`](@ref).
For Brauer characters there is no generic method for checking irreducibility.
"""
function is_irreducible(chi::GAPGroupClassFunction)
    return GAPWrap.IsIrreducibleCharacter(chi.table, chi.values)
end

# Apply a class function to a group element.
function(chi::GAPGroupClassFunction)(g::GAPGroupElem)
    tbl = chi.table
    if tbl.characteristic != 0
      known, fus = known_class_fusion(tbl, get_attribute(tbl, :ordinary_table))
      @assert known "the class fusion is not stored"
    else
      fus = 1:length(tbl)
    end

    # Identify the conjugacy class of `g`.
    ccl = conjugacy_classes(tbl.GAPGroup)
    for i in 1:length(fus)
      if g in ccl[fus[i]]
        return chi[i]
      end
    end
    error("$g is not an element in the underlying group")
end

@doc Markdown.doc"""
    class_positions_of_kernel(chi::GAPGroupClassFunction)

Return the array of those integers `i` such that `chi[i] == chi[1]` holds.
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

@doc Markdown.doc"""
    indicator(chi::GAPGroupClassFunction, n::Int = 2)

Return the `n`-th Frobenius-Schur indicator of `chi`, that is,
the value $(∑_{g ∈ G} chi(g^n))/|G|$, where $G$ is the group of `chi`.

If `chi` is irreducible then `indicator(chi)` is
`0` if `chi` is not real-valued,
`1` if `chi` is afforded by a real representation of $G$, and
`-1` if `chi` is real-valued but not afforded by a real representation of $G$.
"""
function indicator(chi::GAPGroupClassFunction, n::Int = 2)
    ind = GAP.Globals.Indicator(chi.table.GAPTable, GapObj([chi.values]), n)::GapObj
    return ind[1]::Int
end

@doc Markdown.doc"""
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
      F, z = Oscar.AbelianClosure.cyclotomic_field(FF, N)
      f = x::nf_elem -> QQAbElem(x, N)
      finv = function(x::QQAbElem)
        g = gcd(x.c, N)
        K, = Oscar.AbelianClosure.cyclotomic_field(FF, g)
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
      v = Vector{fmpq}(gapcoeffs)
      R, = PolynomialRing(QQ, "x")
      f = R(v)
      F, z = NumberField(f, "z"; cached = true, check = false)
      K, zz = Oscar.AbelianClosure.cyclotomic_field(FF, N)

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
        Kg, = Oscar.AbelianClosure.cyclotomic_field(FF, g)
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

@doc Markdown.doc"""
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
      bound = fmpz(2)
    else
      # Compute the conductor of the largest cyclotomic field
      # that is contained in the character field of `chi`.
      gapfield = GAP.Globals.Field(values)::GapObj
      N = GAPWrap.Conductor(gapfield)
      for n in reverse(sort(divisors(N)))
        if GAP.Globals.E(n)::GAP.Obj in gapfield
          if isodd(n)
            bound = fmpz(2*n)
          else
            bound = fmpz(n)
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
      bound = gcd(bound, scalar_product(fmpz, chi, psi))
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
          bound = gcd(bound, scalar_product(fmpz, chi, psi))
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
          bound = gcd(bound, scalar_product(fmpz, chi, cand[i] * cand[j]))
          bound == 1 && return 1
        end
      end
    end

    if isdefined(tbl, :GAPGroup) && hasproperty(GAP.Globals, :SchurIndexByCharacter)
      g = tbl.GAPGroup
      return GAP.Globals.SchurIndexByCharacter(GAP.Globals.Rationals, g.X, values)
    end

    # For the moment, we do not have more character theoretic criteria.
    error("cannot determine the Schur index with the currently used criteria")
end

@doc Markdown.doc"""
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
            for chi in GAP.Globals.Symmetrizations(tbl.GAPTable,
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

@doc Markdown.doc"""
    symmetric_parts(characters::Vector{GAPGroupClassFunction}, n::Int)

Return the vector of symmetrizations of `characters`
with the trivial character of the symmetric group of degree `n`,
see [`symmetrizations`](@ref).
"""
function symmetric_parts(characters::Vector{GAPGroupClassFunction}, n::Int)
    length(characters) == 0 && return eltype(typeof(characters))[]
    tbl = characters[1].table
    return [group_class_function(tbl, chi)
            for chi in GAP.Globals.SymmetricParts(tbl.GAPTable,
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

@doc Markdown.doc"""
    anti_symmetric_parts(characters::Vector{GAPGroupClassFunction}, n::Int)

Return the vector of symmetrizations of `characters`
with the sign character of the symmetric group of degree `n`,
see [`symmetrizations`](@ref).
"""
function anti_symmetric_parts(characters::Vector{GAPGroupClassFunction}, n::Int)
    length(characters) == 0 && return eltype(typeof(characters))[]
    tbl = characters[1].table
    return [group_class_function(tbl, chi)
            for chi in GAP.Globals.AntiSymmetricParts(tbl.GAPTable,
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

@doc Markdown.doc"""
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
      GAP.Globals.AntiSymmetricParts(tbl.GAPTable, GAP.Obj([chi.values]), n)[1])
end

@doc Markdown.doc"""
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
      GAP.Globals.SymmetricParts(tbl.GAPTable, GAP.Obj([chi.values]), n)[1])
end

@doc Markdown.doc"""
    orthogonal_components(characters::Vector{GAPGroupClassFunction}, n::Int)

Return the vector of the so-called Murnaghan components of the $m$-th
tensor powers of the entries of `characters`, for $m$ up to `n`,
where `n` must be at least 2 and at most 6
and where we assume that the entries of `characters` are irreducible
characters with Frobenius-Schur indicator +1, see [`indicator`](@ref).
"""
function orthogonal_components(characters::Vector{GAPGroupClassFunction}, n::Int)
    (2 <= n && n <= 6) || error("the second ergument must be in 2:6")
    length(characters) == 0 && return eltype(typeof(characters))[]
    tbl = characters[1].table
    return [group_class_function(tbl, chi)
            for chi in GAP.Globals.OrthogonalComponents(tbl.GAPTable,
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

@doc Markdown.doc"""
    symplectic_components(characters::Vector{GAPGroupClassFunction}, n::Int)

Return the vector of the Murnaghan components of the $m$-th tensor powers
of the entries of `characters`, for $m$ up to `n`,
where `n` must be at least 2 and at most 6
and where we assume that the entries of `characters` are irreducible
characters with Frobenius-Schur indicator -1, see [`indicator`](@ref).
"""
function symplectic_components(characters::Vector{GAPGroupClassFunction}, n::Int)
    (2 <= n && n <= 6) || error("the second ergument must be in 2:6")
    length(characters) == 0 && return eltype(typeof(characters))[]
    tbl = characters[1].table
    return [group_class_function(tbl, chi)
            for chi in GAP.Globals.SymplecticComponents(tbl.GAPTable,
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

function character_table_complex_reflection_group(m::Int, p::Int, n::Int)
    p == 1 || error("the case G(m,p,n) with p != 1 is not (yet) supported")
    tbl = GAP.Globals.CharacterTableWreathSymmetric(
            GAP.Globals.CharacterTable(GapObj("Cyclic"), m)::GapObj, n):GapObj
    tbl = GAPGroupCharacterTable(tbl, 0)
    set_attribute!(tbl, :type, (m, p, n))

    return tbl
end
