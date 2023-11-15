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
@doc raw"""
    atlas_irrationality([F::AnticNumberField, ]description::String)

Return the value encoded by `description`.
If `F` is given and is a cyclotomic field that contains the value then
the result is in `F`,
if `F` is not given then the result has type `QQAbElem`.

`description` is assumed to have the format defined in
[CCNPW85](@cite), Chapter 6, Section 10.

# Examples
```jldoctest
julia> Oscar.with_unicode() do
         show(atlas_irrationality("r5"))
       end;
-2*ζ(5)^3 - 2*ζ(5)^2 - 1

julia> atlas_irrationality(cyclotomic_field(5)[1], "r5")
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
    return F(GAPWrap.AtlasIrrationality(GapObj(description)))
end

function atlas_irrationality(description::String)
    F = abelian_closure(QQ)[1]
    return F(GAPWrap.AtlasIrrationality(GapObj(description)))
end


@doc raw"""
    atlas_description(val::QQAbElem)

Return a string in the format defined in
[CCNPW85](@cite), Chapter 6, Section 10,
describing `val`.
Applying [`atlas_irrationality`](@ref) to the result yields `val`.

# Examples
```jldoctest
julia> K, z = abelian_closure(QQ);

julia> val = z(5) + z(5)^4;

julia> str = Oscar.atlas_description(val)
"b5"

julia> val == atlas_irrationality(str)
true
```
"""
function atlas_description(val::QQAbElem)
    iso = Oscar.iso_oscar_gap(parent(val))
    val = iso(val)
    return string(GAP.Globals.CTblLib.StringOfAtlasIrrationality(val))
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
    characteristic::T where T <: IntegerUnion
    group::Union{GAPGroup, GrpAbFinGen}    # the underlying group, if any
    isomorphism::Map  # isomorphism from `group` to a group in GAP

    function GAPGroupCharacterTable(G::Union{GAPGroup, GrpAbFinGen}, tab::GapObj, iso::Map, char::T) where T <: IntegerUnion
      return new(tab, char, G, iso)
    end

    function GAPGroupCharacterTable(tab::GapObj, char::T) where T <: IntegerUnion
      # group and isomorphism are left undefined
      return new(tab, char)
    end
end

# access to field values via functions
GAPTable(tbl::GAPGroupCharacterTable) = tbl.GAPTable

@doc raw"""
    characteristic(::Type{T} = Int, tbl::GAPGroupCharacterTable) where T <: IntegerUnion

Return `T(0)` if `tbl` is an ordinary character table,
and `T(p)` if `tbl` is a `p`-modular character table.

# Examples
```jldoctest
julia> tbl = character_table("A5");

julia> characteristic(tbl)
0

julia> characteristic(tbl % 2)
2
```
"""
characteristic(tbl::GAPGroupCharacterTable) = characteristic(Int, tbl)

function characteristic(::Type{T}, tbl::GAPGroupCharacterTable) where T <: IntegerUnion
    return T(tbl.characteristic)
end


function group(tbl::GAPGroupCharacterTable)
  @req isdefined(tbl, :group) "character table stores no group"
  return tbl.group
end

# backwards compatibility:
# return `group(tbl)` when one accesses `tbl.GAPGroup`
function Base.getproperty(tbl::GAPGroupCharacterTable, sym::Symbol)
   sym === :GAPGroup && return group(tbl)
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
    return MapFromFunc(G, G.X, f, finv)
end

function isomorphism_to_GAP_group(G::GrpAbFinGen)
    @req isfinite(G) "the group is not finite"
    iso = isomorphism(PcGroup, G)
    C = codomain(iso)
    @assert C isa GAPGroup
    f = function(x) return iso(x).X; end
    finv = function(x::GAP.Obj) return preimage(iso, group_element(C, x)); end
    return MapFromFunc(G, C.X, f, finv)
end

function isomorphism_to_GAP_group(tbl::GAPGroupCharacterTable)
    @req isdefined(tbl, :group) "character table stores no group"
    @assert isdefined(tbl, :isomorphism) "character table with group must store also an isomorphism"
    return tbl.isomorphism
end

@doc raw"""
    ordinary_table(tbl::GAPGroupCharacterTable)

Return the ordinary character table of `tbl`,
provided that `tbl` is a Brauer character table.

# Examples
```jldoctest
julia> tbl = character_table("A5");

julia> ordinary_table(tbl % 2) === tbl
true
```
"""
function ordinary_table(tbl::GAPGroupCharacterTable)
    @req characteristic(tbl) != 0 "the character table must be a Brauer table"
    return get_attribute(tbl, :ordinary_table)
end


@doc raw"""
    conjugacy_classes(tbl::GAPGroupCharacterTable)

Return the vector of conjugacy classes of `group(tbl)`,
ordered such that they correspond to the columns of `tbl`
and to the `GAPWrap.ConjugacyClasses` value of the underlying
GAP character table.
Note that the vector `conjugacy_classes(group(tbl))` can be independent
of the vector of conjugacy classes stored in the group of the underlying
GAP character table.

An error is thrown if `tbl` does  not store a group.

# Examples
```jldoctest
julia> g = symmetric_group(4);  tbl = character_table(g);

julia> [length(c) for c in conjugacy_classes(tbl)] == class_lengths(tbl)
true
```
"""
@attr function conjugacy_classes(tbl::GAPGroupCharacterTable)
    G = group(tbl)

    # If `GAPTable(tbl)` does not yet store conjugacy classes
    # then compute them.
    if characteristic(tbl) == 0
      ccl = GAPWrap.ConjugacyClasses(GAPTable(tbl))::GapObj
      ccl = [conjugacy_class(G,
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
    character_table(G::GAPGroup, p::T = 0) where T <: IntegerUnion

Return the ordinary (if `p == 0`) or `p`-modular character table of the
finite group `G`.
If the `p`-modular character table of `G` cannot be computed by GAP
then `nothing` is returned.

# Examples
```jldoctest
julia> Oscar.with_unicode() do
         show(stdout, MIME("text/plain"), character_table(symmetric_group(3)))
       end;
Character table of permutation group of degree 3 and order 6

 2  1  1  .
 3  1  .  1

   1a 2a 3a
2P 1a 1a 3a
3P 1a 2a 1a

χ₁  1 -1  1
χ₂  2  . -1
χ₃  1  1  1

julia> Oscar.with_unicode() do
         show(stdout, MIME("text/plain"), character_table(symmetric_group(3), 2))
       end;
2-modular Brauer table of permutation group of degree 3 and order 6

 2  1  .
 3  1  1

   1a 3a
2P 1a 3a
3P 1a 1a

χ₁  1  1
χ₂  2 -1
```
"""
function character_table(G::Union{GAPGroup, GrpAbFinGen}, p::T = 0) where T <: IntegerUnion
    tbls = get_attribute!(() -> Dict{Int,Any}(), G, :character_tables)
    return get!(tbls, p) do
      p != 0 && return mod(character_table(G, 0), p)
      iso = isomorphism_to_GAP_group(G)
      gaptbl = GAPWrap.CharacterTable(codomain(iso))
      return GAPGroupCharacterTable(G, gaptbl, iso, p)
    end
end

# A character table with stored group object is stored in this group.
# Character tables from the table library do not store groups;
# the ordinary ones are cached in the dictionary `character_tables_by_id`,
# in order to achieve that fetching the same table twice yields the same
# object, and the modular ones are stored in their ordinary tables.
const character_tables_by_id = Dict{String, Union{GAPGroupCharacterTable, Nothing}}()

"""
    character_table(id::String, p::Int = 0)

Return the ordinary (if `p == 0`) or `p`-modular character table
for which `id` is an admissible name in GAP's library of character tables.
If no such table is available then `nothing` is returned.

# Examples
```jldoctest
julia> println(character_table("A5"))
character table of A5

julia> println(character_table("A5", 2))
2-modular Brauer table of A5

julia> println(character_table("J5"))
nothing
```
"""
function character_table(id::String, p::Int = 0)
    p != 0 && return mod(character_table(id, 0), p)
    return get!(character_tables_by_id, id) do
      tbl = GAPWrap.CharacterTable(GapObj(id))
      tbl === GAP.Globals.fail && return nothing
      return GAPGroupCharacterTable(tbl, 0)
    end
end


"""
    character_table(series::Symbol, parameter::Any)

Return the ordinary character table of the group described by the series
`series` and the parameter `parameter`.

# Examples
```jldoctest
julia> println(character_table(:Symmetric, 5))
character table of Sym(5)

julia> println(character_table(:WeylB, 3))
character table of W(B3)
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
| `Symbol("P:Q")` | vector `[p, q]` with prime `p` and `q` dividing `p-1` |
| `:ExtraspecialPlusOdd` | odd power of odd prime |
"""
function character_table(series::Symbol, parameter::Union{Int, Vector{Int}})
    paras = GAP.Obj(parameter, recursive = true)
    tbl = GAPWrap.CharacterTable(GapObj(series), paras)
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
  return sub(domain(iso),
             GrpAbFinGenElem[preimage(iso, x) for x in GAPWrap.GeneratorsOfGroup(H)])
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

Return whether `tbl` is an ordinary table from the character table library
that was constructed from another library character table by permuting rows
and columns.

One application of this function is to restrict the search with
[`all_character_table_names`](@ref) to only one library character table for each
class of permutation equivalent tables.


# Examples
```jldoctest
julia> is_duplicate_table(character_table("A5"))
false

julia> is_duplicate_table(character_table("A6M2"))
true
```
"""
@gapattribute is_duplicate_table(tbl::GAPGroupCharacterTable) = GAP.Globals.IsDuplicateTable(GAPTable(tbl))::Bool

"""
    all_character_table_names(L...; ordered_by = nothing)

Return an vector of strings that contains all those names of character tables
in the character table library that satisfy the conditions in the vector `L`.

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
      K = GAP.Globals.AllCharacterTableNames(gapargs...;
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
        if pos !== nothing
          m[i,j] = legend[pos][2]
        else
          pos = findnext(x -> x[1] == -val, legend, 1)
          if pos !== nothing
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
              if info !== nothing
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

# supercompact and one-line printing
function Base.show(io::IO, tbl::GAPGroupCharacterTable)
  if get(io, :supercompact, false)
    # no nested printing
    if characteristic(tbl) == 0
      print(io, "character table of a group")
    else
      print(io, "$(characteristic(tbl))-modular Brauer table of a group")
    end
  else
    # nested printing allowed, preferably supercompact
    if isdefined(tbl, :group)
      if characteristic(tbl) == 0
        print(io, "character table of ")
      else
        print(io, "$(characteristic(tbl))-modular Brauer table of ")
      end
      io = pretty(io)
      print(IOContext(io, :supercompact => true), Lowercase(), group(tbl))
    elseif characteristic(tbl) == 0
      print(io, "character table of ", identifier(tbl))
    else
      print(io, "$(characteristic(tbl))-modular Brauer table of ",
        identifier(ordinary_table(tbl)))
    end
  end
end

# Produce LaTeX output if `"text/latex"` is prescribed,
# via the `:TeX` attribute of the io context.
function Base.show(io::IO, ::MIME"text/latex", tbl::GAPGroupCharacterTable)
  show(IOContext(io, :TeX => true), MIME("text/plain"), tbl)
end

@doc raw"""
    Base.show(io::IO, ::MIME"text/plain", tbl::GAPGroupCharacterTable)

Display the irreducible characters of `tbl` and context information
as a two-dimensional array.

- First a *header* is shown.
  If `tbl` stores a group then the header describes this group,
  otherwise it is equal to the
  [`identifier(tbl::GAPGroupCharacterTable)`](@ref) value of `tbl`.

- Then the irreducible characters of `tbl` are shown in column portions
  that fit on the screen,
  together with *column labels* above each portion
  and *row labels* on the left of each portion.

  The column labels consist of
  the factored centralizer orders (see [`orders_centralizers`](@ref),
  one row for each prime divisor of the group order),
  followed by one row showing the class names (see [`class_names`](@ref)),
  followed by the power maps (one row for each stored power map).

  The row labels are `X_1`, `X_2`, ...
  (or χ with subscripts `1`, `2`, ... if unicode output is allowed).
  If `io` is an `IOContext` with key `:indicator` set to `true` then
  a second column of row labels shows the 2nd Frobenius-Schur indicator
  of the irreducibles (see [`indicator`](@ref));
  analogously, setting the key `:character_field` to `true` yields a
  column showing the degrees of the character fields
  (see [`character_field`](@ref)),
  and setting the key `:OD` to `true` yields a column showing the known
  orthogonal discriminants of those irreducibles that have indicator `+`
  and even degree.

- Depending on the way how irrational character values are shown,
  a *footer* may be shown in the end.
  By default, irrationalities are shown as sums of roots of unity,
  where `z_n` (or ζ with subscript `n` if unicode output is allowed)
  denotes the primitive `n`-th root $\exp(2 \pi i/n)$.
  If `io` is an `IOContext` with key `:with_legend` set to `true`
  then irrationalities are abbreviated as `A`, `B`, ...,
  and these names together with the corresponding expression as
  sums of roots of unity appear in the footer.

Output in ``\LaTeX`` syntax can be created by calling `show`
with second argument `MIME("text/latex")`.

# Examples
```jldoctest
julia> tbl = character_table(:Cyclic, 3);

julia> Oscar.with_unicode() do
         show(stdout, MIME("text/plain"), tbl)
       end;
C3

 3  1       1       1
                     
   1a      3a      3b
3P 1a      1a      1a
                     
χ₁  1       1       1
χ₂  1      ζ₃ -ζ₃ - 1
χ₃  1 -ζ₃ - 1      ζ₃

julia> Oscar.with_unicode() do
         show(IOContext(stdout, :with_legend => true), MIME("text/plain"), tbl)
       end;
C3

 3  1  1  1
           
   1a 3a 3b
3P 1a 1a 1a
           
χ₁  1  1  1
χ₂  1  A  A̅
χ₃  1  A̅  A

A = ζ₃
A̅ = -ζ₃ - 1

julia> Oscar.with_unicode() do
         show(IOContext(stdout, :indicator => true), MIME("text/plain"), tbl)
       end;
C3

    3  1       1       1
                        
      1a      3a      3b
   3P 1a      1a      1a
    2                   
χ₁  +  1       1       1
χ₂  o  1      ζ₃ -ζ₃ - 1
χ₃  o  1 -ζ₃ - 1      ζ₃

julia> Oscar.with_unicode() do
         show(IOContext(stdout, :character_field => true), MIME("text/plain"), tbl)
       end;
C3

    3  1       1       1
                        
      1a      3a      3b
   2P 1a      3b      3a
   3P 1a      1a      1a
    d                   
χ₁  1  1       1       1
χ₂  2  1      ζ₃ -ζ₃ - 1
χ₃  2  1 -ζ₃ - 1      ζ₃

julia> Oscar.with_unicode() do
         show(IOContext(stdout, :with_legend => true), MIME("text/latex"), tbl)
       end;
C3

$\begin{array}{rrrr}
3 & 1 & 1 & 1 \\
 &  &  &  \\
 & 1a & 3a & 3b \\
2P & 1a & 3b & 3a \\
3P & 1a & 1a & 1a \\
 &  &  &  \\
\chi_{1} & 1 & 1 & 1 \\
\chi_{2} & 1 & A & \overline{A} \\
\chi_{3} & 1 & \overline{A} & A \\
\end{array}

\begin{array}{l}
A = \zeta_{3} \\
\overline{A} = -\zeta_{3} - 1 \\
\end{array}
$
```
"""
function Base.show(io::IO, ::MIME"text/plain", tbl::GAPGroupCharacterTable)
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
    names = class_names(tbl)
    pmaps = Vector{Any}(GAPWrap.ComputedPowerMaps(gaptbl))
    power_maps_primes = String[]
    power_maps_strings = Vector{String}[]
    for i in 2:length(pmaps)
      map = pmaps[i]
      if map !== nothing
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
    # (This is possible only if the OD database is available.)
    OD = get(io, :OD, false)::Bool
    if OD
      ODs = [orthogonal_discriminants(tbl)]
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
        field_degrees = [[string(degree(character_field(x)[1])) for x in tbl]]
      else
        field_degrees = [[string(collect(factor(order_field_of_definition(x)))[1][2]) for x in tbl]]
      end
      field_label = ["d"]
      push!(emptycor, "")
    else
      field_degrees = []
      field_label = []
    end

    emptycol = ["" for i in 1:n]

    if isdefined(tbl, :group)
      headerstring = lowercasefirst(string(group(tbl)))
      if characteristic(tbl) != 0
        headerstring = "$(characteristic(tbl))-modular Brauer table of $(headerstring)"
      else
        headerstring = "Character table of $(headerstring)"
      end
    else
      headerstring = identifier(tbl)
    end

    # Create the IO context.
    ioc = IOContext(io,
      # header (a vector of strings):
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

      # footer (a vector of strings)
      :footer => length(legend) == 0 ? [] :
                 vcat([""], [triple[2]*" = "*triple[3] for triple in legend]),
    )

    # print the table
    labelled_matrix_formatted(ioc, mat)
end


##############################################################################
#
length(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(GAPTable(tbl))::Int
nrows(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(GAPTable(tbl))::Int
ncols(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(GAPTable(tbl))::Int
number_conjugacy_classes(tbl::GAPGroupCharacterTable) = GAPWrap.NrConjugacyClasses(GAPTable(tbl))::Int

@doc raw"""
    order(::Type{T} = ZZRingElem, tbl::GAPGroupCharacterTable) where T <: IntegerUnion

Return the order of the group for which `tbl` is the character table,
as an instance of `T`.

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

Return the vector of the orders of conjugacy class representatives for `tbl`,
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

Return the vector of the orders of centralizers of conjugacy class
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

Return either nothing (if the value is not known) or a vector of identifiers
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
  if GAPWrap.HasMaxes(GAPTable(tbl))
    return Vector{String}(GAPWrap.Maxes(GAPTable(tbl)))
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
    class_positions_of_normal_subgroups(tbl::GAPGroupCharacterTable)

Return a vector whose entries describe all normal subgroups of the group
of `tbl`.
Each entry is a vector of those class positions such that the union of
these classes forms a normal subgroup.

# Examples
```jldoctest
julia> t = character_table("2.A5");

julia> class_positions_of_normal_subgroups(t)
3-element Vector{Vector{Int64}}:
 [1]
 [1, 2]
 [1, 2, 3, 4, 5, 6, 7, 8, 9]
```
"""
function class_positions_of_normal_subgroups(tbl::GAPGroupCharacterTable)
    return Vector{Vector{Int}}(GAPWrap.ClassPositionsOfNormalSubgroups(GAPTable(tbl)))
end


@doc raw"""
    class_positions_of_center(tbl::GAPGroupCharacterTable)

Return the vector of integers ``i`` such that the ``i``-th conjugacy class
of `tbl` is contained in the center of the group $G$ of `tbl`,
i.e., the ``i``-th class has length `1`.

# Examples
```jldoctest
julia> tbl = character_table(dihedral_group(8));

julia> println(class_positions_of_center(tbl))
[1, 4]
```
"""
function class_positions_of_center(tbl::GAPGroupCharacterTable)
    return Vector{Int}(GAPWrap.ClassPositionsOfCentre(GAPTable(tbl)))
end


@doc raw"""
    class_positions_of_derived_subgroup(tbl::GAPGroupCharacterTable)

Return the vector of integers ``i`` such that the ``i``-th conjugacy class
of `tbl` is contained in the derived subgroup of the group $G$ of `tbl`.

# Examples
```jldoctest
julia> tbl = character_table(dihedral_group(8));

julia> println(class_positions_of_derived_subgroup(tbl))
[1, 4]
```
"""
function class_positions_of_derived_subgroup(tbl::GAPGroupCharacterTable)
    return Vector{Int}(GAPWrap.ClassPositionsOfDerivedSubgroup(GAPTable(tbl)))
end


@doc raw"""
    class_positions_of_solvable_residuum(tbl::GAPGroupCharacterTable)

Return the vector of integers ``i`` such that the ``i``-th conjugacy class
of `tbl` is contained in the solvable residuum of the group $G$ of `tbl`,
i.e., the smallest normal subgroup $N$ of $G$ such that the factor group
$G/N$ is solvable.
This normal subgroup is equal to the last term of the derived series of $G$,
see [`derived_series`](@ref).

# Examples
```jldoctest
julia> tbl = character_table(symmetric_group(4));

julia> println(class_positions_of_solvable_residuum(tbl))
[1]
```
"""
function class_positions_of_solvable_residuum(tbl::GAPGroupCharacterTable)
    return Vector{Int}(GAPWrap.ClassPositionsOfSolvableResiduum(GAPTable(tbl)))
end


@doc raw"""
    class_positions_of_pcore(tbl::GAPGroupCharacterTable, p::IntegerUnion)

Return the vector of integers ``i`` such that the ``i``-th conjugacy class
of `tbl` is contained in the `p`-core of the group of `tbl`,
see [`pcore(G::GAPGroup, p::IntegerUnion)`](@ref).

# Examples
```jldoctest
julia> println(class_positions_of_pcore(character_table("2.A5"), 2))
[1, 2]
```
"""
class_positions_of_pcore(tbl::GAPGroupCharacterTable, p::IntegerUnion) = Vector{Int}(GAPWrap.ClassPositionsOfPCore(GAPTable(tbl), GAP.Obj(p)))

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
    iso = isomorphism_to_GAP_group(tbl)
    t = GAPTable(tbl)
    pcorepos = GAPWrap.ClassPositionsOfPCore(t, GAP.Obj(p))
    P = GAPWrap.NormalSubgroupClasses(t, pcorepos)
    return preimages(iso, P)
end

function class_positions_of_kernel(fus::Vector{Int})
  return filter(i -> fus[i] == fus[1], 1:length(fus))
end

function Base.getindex(tbl::GAPGroupCharacterTable, i::Int)
    irr = GAPWrap.Irr(GAPTable(tbl))
    return class_function(tbl, irr[i])
end
#TODO: cache the irreducibles in the table

# in order to make `tbl[end]` work
Base.lastindex(tbl::GAPGroupCharacterTable) = length(tbl)

# in order to make `findfirst` and `findall` work
function Base.keys(tbl::GAPGroupCharacterTable)
    return keys(1:length(tbl))
end

function Base.getindex(tbl::GAPGroupCharacterTable, i::Int, j::Int)
    irr = GAPWrap.Irr(GAPTable(tbl))
    val = irr[i, j]
    return QQAbElem(val)
end
#TODO: cache the values once they are known?

Base.iterate(tbl::GAPGroupCharacterTable, state = 1) = state > nrows(tbl) ? nothing : (tbl[state], state+1)

"""
    mod(tbl::GAPGroupCharacterTable, p::T) where T <: IntegerUnion
    rem(tbl::GAPGroupCharacterTable, p::T) where T <: IntegerUnion

Return the `p`-modular character table of `tbl`,
or `nothing` if this table cannot be computed.

The syntax `tbl % p` is also supported.

An exception is thrown if `tbl` is not an ordinary character table.

# Examples
```jldoctest
julia> show(character_table("A5") % 2)
2-modular Brauer table of A5
```
"""
function Base.mod(tbl::GAPGroupCharacterTable, p::T) where T <: IntegerUnion
    @req is_prime(p) "p must be a prime integer"
    characteristic(tbl) == 0 || error("tbl mod p only for ordinary table tbl")

    modtbls = get_attribute!(() -> Dict{Int,Any}(), tbl, :brauer_tables)
    if ! haskey(modtbls, p)
      modtblgap = mod(GAPTable(tbl), GAP.Obj(p))::GapObj
      if modtblgap === GAP.Globals.fail
        modtbls[p] = nothing
      elseif isdefined(tbl, :group)
        modtbls[p] = GAPGroupCharacterTable(group(tbl), modtblgap,
                       isomorphism_to_GAP_group(tbl), p)
        set_attribute!(modtbls[p], :ordinary_table, tbl)
      else
        modtbls[p] = GAPGroupCharacterTable(modtblgap, p)
        set_attribute!(modtbls[p], :ordinary_table, tbl)
      end
    end

    return modtbls[p]
end

rem(tbl::GAPGroupCharacterTable, p::T) where T <: IntegerUnion = mod(tbl, p)

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
    return matrix(ZZ, GAPWrap.DecompositionMatrix(GAPTable(modtbl)))
end


@doc raw"""
    quo(tbl::GAPGroupCharacterTable, nclasses::Vector{Int})

Return the pair `(fact, proj)` where `fact` is the character table of the
factor of `tbl` modulo the normal subgroup that is the normal closure of
the conjugacy classes whose positions are listed in `nclasses`,
and `proj` is the class fusion from `tbl` to `fact`.

# Examples
```jldoctest
julia> t = character_table("2.A5");

julia> n = class_positions_of_center(t);  println(n)
[1, 2]

julia> fact, proj = quo(t, n);

julia> order(fact)
60

julia> println(proj)
[1, 1, 2, 3, 3, 4, 4, 5, 5]
```
"""
function quo(tbl::GAPGroupCharacterTable, nclasses::Vector{Int})
  @req characteristic(tbl) == 0 "supported only for ordinary character tables"
  gap_fact = GAPWrap.CharacterTableFactorGroup(GAPTable(tbl), GapObj(nclasses))
  fact = GAPGroupCharacterTable(gap_fact, 0)
  flag, fus = known_class_fusion(tbl, fact)
  @assert flag
  return fact, fus
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
  return T(GAPWrap.ClassMultiplicationCoefficient(GAPTable(tbl), i, j, k))
end

class_multiplication_coefficient(tbl::GAPGroupCharacterTable, i::Int, j::Int, k::Int) = class_multiplication_coefficient(ZZRingElem, tbl, i, j, k)

@doc raw"""
    approximate_class_fusion(subtbl::GAPGroupCharacterTable,
                             tbl::GAPGroupCharacterTable)

Compute for each class of `subtbl` all those classes in `tbl` to which it can
fuse under an embedding of the group of `subtbl` into the group of `tbl`,
according to element orders and centralizer orders in the two tables.

If no embedding is possible then return an empty vector.
Otherwise return a vector of length equal to the number of classes of
`subtbl`, such that the entry at position $i$ either an integer (if there is
a unique possible image class) or the vector of the positions of possible
image classes.

# Examples
```jldoctest
julia> subtbl = character_table("A5"); tbl = character_table("A6");

julia> println(approximate_class_fusion(subtbl, tbl))
Union{Int64, Vector{Int64}}[1, 2, [3, 4], [6, 7], [6, 7]]

```
"""
function approximate_class_fusion(subtbl::GAPGroupCharacterTable,
                                  tbl::GAPGroupCharacterTable)
  fus = GAPWrap.InitFusion(GAPTable(subtbl), GAPTable(tbl))
  res = Union{Int, Vector{Int}}[]
  fus == GAP.Globals.fail && return res
  for i in 1:length(fus)
    if fus[i] isa Int
      push!(res, fus[i])
    else
      push!(res, Vector{Int}(fus[i]))
    end
  end
  return res
end


@doc raw"""
    possible_class_fusions(subtbl::GAPGroupCharacterTable,
                           tbl::GAPGroupCharacterTable;
                           decompose::Bool = true,
                           fusionmap::Vector = [])

Return the vector of possible class fusions from `subtbl` to `tbl`.
Each entry is an vector of positive integers, where the value at position `i`
is the position of the conjugacy class in `tbl` that contains the `i`-th class
of `subtbl`.

If `decompose` is set to `true` then the strategy is changed:
The decomposability of restricted characters of `tbl` as integral linear
combinations of characters of `subtbl` (with perhaps negative coefficients)
is not checked;
this does not change the result,
but in certain situations it is faster to omit this step.

If `fusionmap` is set to a vector of integers and integer vectors then
only those maps are returned that are compatible with the prescribed value.

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
function possible_class_fusions(subtbl::GAPGroupCharacterTable,
                                tbl::GAPGroupCharacterTable;
                                decompose::Bool = true,
                                fusionmap::Vector = [])
  cond = Dict{Symbol, Any}(:decompose => decompose)
  if length(fusionmap) != 0
    cond[:fusionmap] = GapObj(fusionmap, recursive = true)
  end
  fus = GAPWrap.PossibleClassFusions(GAPTable(subtbl), GAPTable(tbl),
            GapObj(cond))
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
      return [x for x in GAPWrap.DescriptionOfRootOfUnity(para)]
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
      GAPWrap.HasCharacterParameters(GAPt) || return nothing
      paras = Vector{GAP.Obj}(GAPWrap.CharacterParameters(GAPt))
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
      GAPWrap.HasClassParameters(GAPt) || return nothing
      paras = Vector{GAP.Obj}(GAPWrap.ClassParameters(GAPt))
      return _translate_parameter_list(paras)
    end
end


@doc raw"""
    class_names(tbl::GAPGroupCharacterTable)

Return a vector of strings corresponding to the columns of `tbl`.
The $i$-th entry consists of the element order for the $i$-th column,
followed by at least one distinguishing letter.
For example, the classes of elements of order two have class names
`"2a", "2b", and so on.

# Examples
```jldoctest
julia> println(class_names(character_table("S5")))
["1a", "2a", "3a", "5a", "2b", "4a", "6a"]
```
"""
function class_names(tbl::GAPGroupCharacterTable)
    return get_attribute!(tbl, :class_names) do
      return Vector{String}(GAPWrap.ClassNames(GAPTable(tbl)))
    end
end


@doc raw"""
    names_of_fusion_sources(tbl::GAPGroupCharacterTable)

Return the vector of strings that are identifiers of those character tables
which store a class fusion to `tbl`,
which must be an ordinary character table.

# Examples
```jldoctest
julia> tbl = character_table("A5");

julia> println(maxes(tbl))
["a4", "D10", "S3"]

julia> all(x -> x in names_of_fusion_sources(tbl), maxes(tbl))
true
```
"""
function names_of_fusion_sources(tbl::GAPGroupCharacterTable)
    return [string(name) for name in GAPWrap.NamesOfFusionSources(GAPTable(tbl))]
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
    map = GAPWrap.GetFusionMap(GAPTable(subtbl), GAPTable(tbl))
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

parent(chi::GAPGroupClassFunction) = chi.table

function Base.show(io::IO, chi::GAPGroupClassFunction)
    print(io, "class_function($(parent(chi)), $(values(chi)))")
end

function values(chi::GAPGroupClassFunction)
    gapvalues = GAPWrap.ValuesOfClassFunction(chi.values)
    return [QQAbElem(x) for x in gapvalues]
end

group(chi::GAPGroupClassFunction) = group(parent(chi))

characteristic(chi::GAPGroupClassFunction) = characteristic(parent(chi))

function class_function(tbl::GAPGroupCharacterTable, values::GapObj)
    @req GAPWrap.IsClassFunction(values) "values must be a class function"
    return GAPGroupClassFunction(tbl, values)
end

function class_function(tbl::GAPGroupCharacterTable, values::Vector{<:QQAbElem})
    gapvalues = GapObj([GAP.Obj(x) for x in values])
    return GAPGroupClassFunction(tbl, GAPWrap.ClassFunction(GAPTable(tbl), gapvalues))
end

function class_function(G::GAPGroup, values::GapObj)
    @req GAPWrap.IsClassFunction(values) "values must be a class function"
    return GAPGroupClassFunction(character_table(G), values)
end

function class_function(G::GAPGroup, values::Vector{<:QQAbElem})
    return class_function(character_table(G), values)
end

@doc raw"""
    trivial_character(tbl::GAPGroupCharacterTable)

Return the character of `tbl` that has the value `QQAbElem(1)` in each position.

# Examples
```jldoctest
julia> t = character_table(symmetric_group(4));

julia> all(x -> x == 1, trivial_character(t))
true
```
"""
function trivial_character(tbl::GAPGroupCharacterTable)
    val = QQAbElem(1)
    return class_function(tbl, [val for i in 1:ncols(tbl)])
end

@doc raw"""
    trivial_character(G::GAPGroup)

Return the character of (the ordinary character table of) `G`
that has the value `QQAbElem(1)` in each position.

# Examples
```jldoctest
julia> g = symmetric_group(4);

julia> all(x -> x == 1, trivial_character(g))
true
```
"""
function trivial_character(G::GAPGroup)
    val = QQAbElem(1)
    return class_function(G, [val for i in 1:Int(number_conjugacy_classes(G))])
end

@doc raw"""
    natural_character(G::PermGroup)

Return the permutation character of degree `degree(G)`
that maps each element of `G` to the number of its fixed points.

# Examples
```jldoctest
julia> g = symmetric_group(4);

julia> degree(natural_character(g))
4

julia> degree(natural_character(stabilizer(g, 4)[1]))
4
```
"""
function natural_character(G::PermGroup)
    tbl = character_table(G)
    ccl = conjugacy_classes(tbl)
    FF = abelian_closure(QQ)[1]
    n = degree(G)
    vals = [FF(n - number_moved_points(representative(x))) for x in ccl]
    return class_function(G, vals)
end

@doc raw"""
    natural_character(G::Union{MatrixGroup{QQFieldElem}, MatrixGroup{nf_elem}})

Return the character that maps each element of `G` to its trace.
We assume that the entries of the elements of `G` are either of type `QQFieldElem`
or contained in a cyclotomic field.

# Examples
```jldoctest
julia> g = matrix_group(matrix(ZZ, [0 1; 1 0]));

julia> println(values(natural_character(g)))
QQAbElem{nf_elem}[2, 0]
```
"""
function natural_character(G::Union{MatrixGroup{ZZRingElem}, MatrixGroup{QQFieldElem}, MatrixGroup{nf_elem}})
    tbl = character_table(G)
    ccl = conjugacy_classes(tbl)
    FF = abelian_closure(QQ)[1]
    vals = [FF(tr(representative(x))) for x in ccl]
    return class_function(G, vals)
end

@doc raw"""
    natural_character(G::MatrixGroup{FinFieldElem})

Return the character that maps each $p$-regular element of `G`,
where $p$ is the characteristic of the base field of `G`,
to its Brauer character value.

# Examples
```jldoctest
julia> g = general_linear_group(2, 2);

julia> println(values(natural_character(g)))
QQAbElem{nf_elem}[2, -1]
```
"""
function natural_character(G::MatrixGroup{T, MT}) where T <: FinFieldElem where MT
    p = characteristic(base_ring(G))
    tbl = character_table(G, p)
    ccl = conjugacy_classes(tbl)
    vals = [GAP.Globals.BrauerCharacterValue(representative(x).X) for x in ccl]
    vals = GAPWrap.ClassFunction(GAPTable(tbl), GapObj(vals))

    return class_function(G, vals)
end

@doc raw"""
    natural_character(rho::GAPGroupHomomorphism)

Return the character of `domain(rho)` that is afforded by the representation
`rho`, where `codomain(rho)` must be a permutation group or a matrix group.
In the latter case, an ordinary character is returned if the characteristic
of the base field is zero, and a $p$-modular Brauer character is returned
if the characteristic is $p > 0$.

# Examples
```jldoctest
julia> g = symmetric_group(3);  h = general_linear_group(2, 2);

julia> mp = hom(g, h, [g([2,1]), g([1, 3, 2])], gens(h));

julia> println(values(natural_character(mp)))
QQAbElem{nf_elem}[2, -1]
```
"""
function natural_character(rho::GAPGroupHomomorphism)
    G = domain(rho)
    M = codomain(rho)
    tbl = character_table(G)
    FF = abelian_closure(QQ)[1]

    if M isa PermGroup
      modtbl = tbl
      ccl = conjugacy_classes(tbl)
      n = degree(M)
      vals = [FF(n - number_moved_points(rho(representative(x)))) for x in ccl]
    elseif M isa MatrixGroup
      p = characteristic(base_ring(M))
      if p == 0
        # ordinary character
        modtbl = tbl
        ccl = conjugacy_classes(modtbl)
        vals = [FF(tr(rho(representative(x)))) for x in ccl]
      else
        # Brauer character
        modtbl = mod(tbl, p)
        ccl = conjugacy_classes(modtbl)  # p-regular classes
        vals = [GAP.Globals.BrauerCharacterValue(rho(representative(x)).X) for x in ccl]
        vals = GAPWrap.ClassFunction(GAPTable(modtbl), GapObj(vals))
      end
    else
      throw(ArgumentError("codomain must be a PermGroup or MatrixGroup"))
    end

    return class_function(modtbl, vals)
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
  known, fus = known_class_fusion(parent(chi), tbl)
  known && return induce(chi, tbl, fus)

  subtbl = parent(chi)
  if !(isdefined(tbl, :group) && isdefined(subtbl, :group))
    # If there is no stored group then let GAP try to find the result.
    fus = GAPWrap.FusionConjugacyClasses(GAPTable(subtbl), GAPTable(tbl))
    @req fus !== GAP.Globals.fail "class fusion is not uniquely determinaed"
    ind = GAPWrap.InducedClassFunctionsByFusionMap(GAPTable(subtbl),
          GAPTable(tbl), GapObj([chi.values]), GapObj(fus))
    return GAPGroupClassFunction(tbl, ind[1])
  else
    # Dispatch on the types of the stored groups.
    return _induce(chi, tbl, group(subtbl), group(tbl))
  end
end

function induce(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable, fusion::Vector{Int})
  ind = GAPWrap.InducedClassFunctionsByFusionMap(GAPTable(parent(chi)),
          GAPTable(tbl), GapObj([chi.values]), GapObj(fusion))
  return GAPGroupClassFunction(tbl, ind[1])
end

# If `GAPGroup` groups are stored then we let GAP try to find the result.
function _induce(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable, G_chi::GAPGroup, G_tbl::GAPGroup)
  ind = GAPWrap.InducedClassFunction(chi.values, GAPTable(tbl))
  return GAPGroupClassFunction(tbl, ind)
end

# If `GrpAbFinGen` groups are stored then we have to work with explicit
# embeddings.
function _induce(chi::GAPGroupClassFunction,
              tbl::GAPGroupCharacterTable,
              G_chi::GrpAbFinGen, G_tbl::GrpAbFinGen)
  H, emb = is_subgroup(G_chi, G_tbl)
  Hreps = [emb(representative(C)) for C in conjugacy_classes(parent(chi))]
  Greps = [representative(C) for C in conjugacy_classes(tbl)]
  fus = [findfirst(x -> x == y, Greps) for y in Hreps]
  return induce(chi, tbl, fus)
end


@doc raw"""
    induced_cyclic(tbl::GAPGroupCharacterTable)

Return the vector of permutation characters of `tbl` that are induced from
cyclic subgroups.
"""
function induced_cyclic(tbl::GAPGroupCharacterTable)
    return [GAPGroupClassFunction(tbl, chi) for chi in GAPWrap.InducedCyclic(GAPTable(tbl))]
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
  tbl = parent(chi)
  known, fus = known_class_fusion(subtbl, tbl)
  known && return restrict(chi, subtbl, fus)

  if !(isdefined(tbl, :group) && isdefined(subtbl, :group))
    # If there is no stored group then let GAP try to find the result.
    fus = GAPWrap.FusionConjugacyClasses(GAPTable(subtbl), GAPTable(tbl))
    @req fus !== GAP.Globals.fail "class fusion is not uniquely determinaed"
    rest = GAPWrap.ELMS_LIST(chi.values, GapObj(fus))
    rest = GAPWrap.ClassFunction(GAPTable(subtbl), rest)
    return GAPGroupClassFunction(subtbl, rest)
  else
    # Dispatch on the types of the stored groups.
    return _restrict(chi, subtbl, group(tbl), group(subtbl))
  end
end

function restrict(chi::GAPGroupClassFunction, subtbl::GAPGroupCharacterTable, fusion::Vector{Int})
  rest = GAPWrap.ELMS_LIST(chi.values, GapObj(fusion))
  rest = GAPWrap.ClassFunction(GAPTable(subtbl), rest)
  return GAPGroupClassFunction(subtbl, rest)
end

# If `GAPGroup` groups are stored then we let GAP try to find the result.
function _restrict(chi::GAPGroupClassFunction, subtbl::GAPGroupCharacterTable, G_chi::GAPGroup, G_tbl::GAPGroup)
  tbl = parent(chi)
  fus = GAPWrap.FusionConjugacyClasses(GAPTable(subtbl), GAPTable(tbl))
  rest = GAPWrap.ELMS_LIST(chi.values, GapObj(fus))
  rest = GAPWrap.ClassFunction(GAPTable(subtbl), rest)
  return GAPGroupClassFunction(subtbl, rest)
end

# If `GrpAbFinGen` groups are stored then we have to work with explicit
# embeddings.
function _restrict(chi::GAPGroupClassFunction,
              subtbl::GAPGroupCharacterTable,
              G_chi::GrpAbFinGen, G_subtbl::GrpAbFinGen)
  H, emb = is_subgroup(G_subtbl, G_chi)
  Hreps = [emb(representative(C)) for C in conjugacy_classes(subtbl)]
  Greps = [representative(C) for C in conjugacy_classes(parent(chi))]
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
  vals = GAPWrap.ValuesOfClassFunction(chi.values)
  return QQAbElem(vals[i])
end

# arithmetic with class functions
function Base.:(==)(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    @req parent(chi) === parent(psi) "character tables must be identical"
    return chi.values == psi.values
end

# Currently we cannot implement a `hash` method based on the values,
# since `hash(::QQAbElem)` is based on `objectid`.
function Base.hash(chi::GAPGroupClassFunction, h::UInt)
  return Base.hash(parent(chi), h)
end

function Base.:+(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    @req parent(chi) === parent(psi) "character tables must be identical"
    return GAPGroupClassFunction(parent(chi), chi.values + psi.values)
end

Base.:-(chi::GAPGroupClassFunction) = GAPGroupClassFunction(parent(chi), - chi.values)

function Base.:-(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    @req parent(chi) === parent(psi) "character tables must be identical"
    return GAPGroupClassFunction(parent(chi), chi.values - psi.values)
end

function Base.:*(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    @req parent(chi) === parent(psi) "character tables must be identical"
    return GAPGroupClassFunction(parent(chi), chi.values * psi.values)
end

function Base.zero(chi::GAPGroupClassFunction)
    val = QQAbElem(0)
    return class_function(parent(chi), [val for i in 1:length(chi)])
end

Base.one(chi::GAPGroupClassFunction) = trivial_character(parent(chi))

@doc raw"""
    scalar_product(::Type{T} = QQFieldElem, chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
                   where T <: Union{IntegerUnion, ZZRingElem, QQFieldElem, QQAbElem}

Return $\sum_{g \in G}$ `chi`($g$) `conj(psi)`($g$) / $|G|$,
where $G$ is the group of both `chi` and `psi`.
The result is an instance of `T`.

Note that we do not support `dot(chi, psi)` and its infix notation
because the documentation of `dot` states that the result is equal to the
sum of `dot` results of corresponding entries,
which does not hold for the scalar product of characters.
"""
scalar_product(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction) = scalar_product(QQFieldElem, chi, psi)

function scalar_product(::Type{T}, chi::GAPGroupClassFunction, psi::GAPGroupClassFunction) where T <: Union{Integer, ZZRingElem, QQFieldElem, QQAbElem}
    @req parent(chi) === parent(psi) "character tables must be identical"
    return T(GAPWrap.ScalarProduct(GAPTable(parent(chi)), chi.values, psi.values))::T
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
Permutation group of degree 4 and order 24

julia> chi = natural_character(g);

julia> coordinates(Int, chi)
5-element Vector{Int64}:
 0
 0
 0
 1
 1

julia> t = parent(chi);  t3 = mod(t, 3);  chi3 = restrict(chi, t3);

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
    t = parent(chi)
    GAPt = GAPTable(t)
    if characteristic(t) == 0
      # use scalar products for an ordinary character
      c = GAPWrap.MatScalarProducts(GAPt, GAPWrap.Irr(GAPt), GapObj([chi.values]))
    else
      # decompose a Brauer character
      c = GAPWrap.Decomposition(GAPWrap.Irr(GAPt), GapObj([chi.values]), GapObj("nonnegative"))
    end
    return Vector{T}(c[1])::Vector{T}
end


@doc raw"""
    multiplicities_eigenvalues(::Type{T} = Int, chi::GAPGroupClassFunction, i::Int) where T <: IntegerUnion

Let $M$ be a representing matrix of an element in the `i`-th conjugacy class
of the character table of `chi`,
in a representation affording the character `chi`,
and let $n$ be the order of the elements in this conjugacy class.

Return the vector $(m_1, m_2, \ldots, m_n)$ of integers of type `T`
such that $m_j$ is the multiplicity of $\zeta_n^j$ as an eigenvalue of $M$.

# Examples
```jldoctest
julia> t = character_table("A5");  chi = t[4];

julia> println(values(chi))
QQAbElem{nf_elem}[4, 0, 1, -1, -1]

julia> println(multiplicities_eigenvalues(chi, 5))
[1, 1, 1, 1, 0]
```
"""
function multiplicities_eigenvalues(::Type{T}, chi::GAPGroupClassFunction, i::Int) where T <: IntegerUnion
    return Vector{T}(GAPWrap.EigenvaluesChar(chi.values, i))
end

multiplicities_eigenvalues(chi::GAPGroupClassFunction, i::Int) = multiplicities_eigenvalues(Int, chi, i)


function Base.:*(n::IntegerUnion, chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(parent(chi), n * chi.values)
end

function Base.:^(chi::GAPGroupClassFunction, n::IntegerUnion)
    return GAPGroupClassFunction(parent(chi), chi.values ^ n)
end

function Base.:^(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable)
    return induce(chi, tbl)
end

function Base.:^(chi::GAPGroupClassFunction, g::Union{GAPGroupElem, GrpAbFinGenElem})
    tbl = parent(chi)
    ccl = conjugacy_classes(tbl)
    reps = [representative(c) for c in ccl]
    pi = [findfirst(x -> x^g in c, reps) for c in ccl]
    return class_function(tbl, values(chi)[pi])
end

@doc raw"""
    conj(chi::GAPGroupClassFunction)

Return the class function whose values are the complex conjugates of
the values of `chi`.

# Examples
```jldoctest
julia> tbl = character_table(alternating_group(4));

julia> println([findfirst(y -> y == conj(x), tbl) for x in tbl])
[1, 3, 2, 4]
```
"""
function conj(chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(parent(chi), GAPWrap.GaloisCyc(chi.values, -1))
end

@doc raw"""
    (sigma::QQAbAutomorphism)(chi::GAPGroupClassFunction)

Return the class function whose values are the images of the values of `chi`
under `sigma`.
"""
function (sigma::QQAbAutomorphism)(chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(parent(chi), GAPWrap.GaloisCyc(chi.values, sigma.exp))
end

Base.:^(chi::GAPGroupClassFunction, sigma::QQAbAutomorphism) = sigma(chi)

@doc raw"""
    is_rational(chi::GAPGroupClassFunction)

Return `true` if all values of `chi` are rational, i.e., in `QQ`,
and `false` otherwise.

# Examples
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

# Examples
```jldoctest
julia> g = symmetric_group(4);

julia> all(is_irreducible, character_table(g))
true

julia> is_irreducible(natural_character(g))
false
```
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
julia> println(map(is_faithful, character_table(symmetric_group(3))))
Bool[0, 1, 0]
```
"""
function is_faithful(chi::GAPGroupClassFunction)
    return length(class_positions_of_kernel(chi)) == 1
end

# Apply a class function to a group element.
function(chi::GAPGroupClassFunction)(g::Union{GAPGroupElem, GrpAbFinGenElem})
    tbl = parent(chi)

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

Return the vector of those integers `i` such that `chi[i] == chi[1]` holds.

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
    tbl = parent(chi)
    iso = isomorphism_to_GAP_group(tbl)
    GAP_K = GAPWrap.KernelOfCharacter(GAPTable(tbl), chi.values)
    return preimages(iso, GAP_K)
end

@doc raw"""
    class_positions_of_center(chi::GAPGroupClassFunction)

Return the vector of those integers `i` such that `chi[i]` is `chi[1]` times
a root of unity.

# Examples
```jldoctest
julia> println(class_positions_of_center(character_table("2.A5")[2]))
[1, 2]
```
"""
function class_positions_of_center(chi::GAPGroupClassFunction)
    return Vector{Int}(GAPWrap.ClassPositionsOfCentre(chi.values))
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
    tbl = parent(chi)
    iso = isomorphism_to_GAP_group(tbl)
    C = GAPWrap.CentreOfCharacter(GAPTable(tbl), chi.values)
    return preimages(iso, C)
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
    values = GAPWrap.DeterminantOfCharacter(chi.values)
    return GAPGroupClassFunction(parent(chi), values)
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
    return T(GAPWrap.Order(det(chi).values))::T
end

@doc raw"""
    indicator(chi::GAPGroupClassFunction, n::Int = 2)

Return the `n`-th Frobenius-Schur indicator of `chi`, that is,
the value $(∑_{g ∈ G} chi(g^n))/|G|$, where $G$ is the group of `chi`.

If `chi` is irreducible then `indicator(chi)` is
`0` if `chi` is not real-valued,
`1` if `chi` is afforded by a real representation of $G$, and
`-1` if `chi` is real-valued but not afforded by a real representation of $G$.

# Examples
```jldoctest
julia> tbl = character_table("U3(3)");

julia> println([indicator(chi) for chi in tbl])
[1, -1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0]
```
"""
function indicator(chi::GAPGroupClassFunction, n::Int = 2)
    tbl = parent(chi)
    if characteristic(tbl) == 0
      # The indicator can be computed for any character.
      ind = GAPWrap.Indicator(GAPTable(tbl), GapObj([chi.values]), n)
      return ind[1]::Int
    else
      # The indicator is defined only for `n = 2` and irreducible characters.
      @req n == 2 "defined for Brauer characters only for n = 2"
      chipos = findfirst(isequal(chi), tbl)
      @req chipos !== nothing "defined only for irreducible Brauer characters"
      ind = GAPWrap.Indicator(GAPTable(tbl), n)
      return ind[chipos]::Int
    end
end

@doc raw"""
    character_field(chi::GAPGroupClassFunction)

If `chi` is an ordinary character then
return the pair `(F, phi)` where `F` is a number field that is generated
by the character values of `chi`, and `phi` is the embedding of `F` into
`abelian_closure(QQ)`.

If `chi` is a Brauer character in characteristic `p` then
return the pair `(F, phi)` where `F` is the finite field that is generated
by the `p`-modular reductions of the values of `chi`,
and `phi` is the identity map on `F`.

# Examples
```jldoctest
julia> t = character_table("A5");

julia> character_field(t[2])[1]
Number field with defining polynomial x^2 + x - 1
  over rational field

julia> flds_2 = map(character_field, mod(t, 2));

julia> println([degree(x[1]) for x in flds_2])
[1, 2, 2, 1]
```
"""
function character_field(chi::GAPGroupClassFunction)
    p = characteristic(chi)
    if p != 0
      # Brauer character, construct a finite field
      q = order_field_of_definition(chi)
      flag, e, pp = is_prime_power_with_data(q)
      (flag && p == pp) || error("something is wrong with 'GAPWrap.SizeOfFieldOfDefinition'")
      F = Nemo._GF(p, e)
      return (F, identity_map(F))
    end

    values = chi.values  # a list of GAP cyclotomics
    gapfield = GAPWrap.Field(values)
    N = GAPWrap.Conductor(gapfield)
    FF, _ = abelian_closure(QQ)
    if GAPWrap.IsCyclotomicField(gapfield)
      # In this case, the want to return a field that knows to be cyclotomic
      # (and the embedding is easy).
      F, z = AbelianClosure.cyclotomic_field(FF, N)
      nfelm = FF(z)
    else
      # In the general case, we have to work for the embedding.
      gapgens = GAPWrap.GeneratorsOfField(gapfield)
      @assert length(gapgens) == 1
      gappol = GAPWrap.MinimalPolynomial(GAP.Globals.Rationals, gapgens[1])
      gapcoeffs = GAPWrap.CoefficientsOfUnivariatePolynomial(gappol)
      v = Vector{QQFieldElem}(gapcoeffs)
      R, = polynomial_ring(QQ, "x")
      f = R(v)
      F, _ = number_field(f, "z"; cached = true, check = false)
      nfelm = QQAbElem(gapgens[1])
    end

    return F, AbelianClosure._embedding(F, FF, nfelm)
end

@doc raw"""
    number_field(::QQField, chi::GAPGroupClassFunction; cached::Bool = false)

Return a pair `(F, a)` where `F` is a number field that describes the
character field of the ordinary character `chi`,
and `a` is a generator of `F`.
The syntax `QQ[chi]` is also supported.

No embedding of `F` into the abelian closure of `QQ` is computed,
use [`character_field`](@ref) if you are interested in the relations
between the character fields of several characters.
"""
function Hecke.number_field(::QQField, chi::GAPGroupClassFunction; cached::Bool = false)
  @req characteristic(chi) == 0 "character must be ordinary"
  return number_field(QQ, values(chi), cached = cached)
end

# support the syntax `QQ[chi]`
Base.getindex(::QQField, chi::Oscar.GAPGroupClassFunction) = number_field(QQ, chi)


@doc raw"""
    order_field_of_definition(::Type{T} = ZZRingElem, chi::GAPGroupClassFunction) where T <: IntegerUnion

Return $p^n$, as an instance of `T`, if `chi` is a $p$-modular Brauer character
such that the $p$-modular reductions of the values of `chi`
span the field with $p^n$ elements.

Note that one need not compute the [`character_field`](@ref) value of `chi`
in order to compute `order_field_of_definition(chi)`.

# Examples
```jldoctest
julia> tbl = character_table("A5", 2);

julia> println([order_field_of_definition(chi) for chi in tbl])
ZZRingElem[2, 4, 4, 2]
```
"""
order_field_of_definition(chi::GAPGroupClassFunction) =
    order_field_of_definition(ZZRingElem, chi)

function order_field_of_definition(::Type{T}, chi::GAPGroupClassFunction) where
  T <: IntegerUnion
    p = characteristic(chi)
    @req p != 0 "the character must be a Brauer character"
    return T(GAPWrap.SizeOfFieldOfDefinition(chi.values, p))
end


@doc raw"""
    schur_index(chi::GAPGroupClassFunction) -> Int

For an ordinary irreducible character `chi`,
return the minimal integer `m` such that the character `m * chi`
is afforded by a representation over the character field of `chi`,
or throw an exception if the currently used character theoretic criteria
do not suffice for computing `m`.

# Examples
```jldoctest
julia> t = character_table(quaternion_group(8));

julia> println(map(schur_index, t))
[1, 1, 1, 1, 2]
```
"""
function schur_index(chi::GAPGroupClassFunction, recurse::Bool = true)
    @req characteristic(chi) == 0 "defined only for ordinary characters"
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
      gapfield = GAPWrap.Field(values)
      N = GAPWrap.Conductor(gapfield)
      for n in reverse(sort(divisors(N)))
        if GAPWrap.E(n) in gapfield
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
    tbl = parent(chi)
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
      # The function is defined in the Wedderga package.
      g = group(tbl)
      return GAPWrap.SchurIndexByCharacter(GAP.Globals.Rationals, codomain(isomorphism_to_GAP_group(tbl)), values)
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
    tbl = parent(characters[1])
    return [class_function(tbl, chi)
            for chi in GAPWrap.Symmetrizations(GAPTable(tbl),
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
    tbl = parent(characters[1])
    return [class_function(tbl, chi)
            for chi in GAPWrap.SymmetricParts(GAPTable(tbl),
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
    tbl = parent(characters[1])
    return [class_function(tbl, chi)
            for chi in GAPWrap.AntiSymmetricParts(GAPTable(tbl),
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
    tbl = parent(chi)
    return class_function(tbl,
      GAPWrap.AntiSymmetricParts(GAPTable(tbl), GAP.Obj([chi.values]), n)[1])
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
    tbl = parent(chi)
    return class_function(tbl,
      GAPWrap.SymmetricParts(GAPTable(tbl), GAP.Obj([chi.values]), n)[1])
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
    tbl = parent(characters[1])
    return [class_function(tbl, chi)
            for chi in GAPWrap.OrthogonalComponents(GAPTable(tbl),
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
    tbl = parent(characters[1])
    return [class_function(tbl, chi)
            for chi in GAPWrap.SymplecticComponents(GAPTable(tbl),
                         GAP.GapObj([chi.values for chi in characters]), n)]
end

function character_table_complex_reflection_group(m::Int, p::Int, n::Int)
    @req p == 1 "the case G(m,p,n) with p != 1 is not (yet) supported"
    tbl = GAPWrap.CharacterTableWreathSymmetric(
            GAPWrap.CharacterTable(GapObj("Cyclic"), m), n)
    tbl = GAPGroupCharacterTable(tbl, 0)
    set_attribute!(tbl, :type, (m, p, n))

    return tbl
end
