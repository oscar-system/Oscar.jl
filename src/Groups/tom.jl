#############################################################################
##
##  tables of marks in Oscar
##
abstract type GroupTableOfMarks end

"""
    GAPGroupTableOfMarks <: GroupTableOfMarks

This is the type of tables of marks that can delegate
tasks to an underlying table of marks object in the GAP system
(field `GAPTable`).

A group can (but need not) be stored in the field `group`.
If it is available then also the field `isomorphism` is available,
its value is a bijective map from the `group` value to a group in GAP.

Objects of type `GAPGroupTableOfMarks` support [`get_attribute`](@ref).
"""
@attributes mutable struct GAPGroupTableOfMarks <: GroupTableOfMarks
  GAPTable::GapObj  # the GAP table of marks object
  group::Union{GAPGroup, FinGenAbGroup}    # the underlying group, if any
  isomorphism::Map  # isomorphism from `group` to a group in GAP

  function GAPGroupTableOfMarks(G::Union{GAPGroup, FinGenAbGroup}, tom::GapObj, iso::Map)
    @req GAPWrap.IsTableOfMarks(tom) "tom must be a GAP table of marks"
    return new(tom, G, iso)
  end

  function GAPGroupTableOfMarks(tom::GapObj)
    # group and isomorphism are left undefined
    @req GAPWrap.IsTableOfMarks(tom) "tom must be a GAP table of marks"
    return new(tom)
  end
end


# If `tblid` is an admissible name for a library character table
# then translate it to the GAP name of the corresponding
# library table of marks if there is one.
function _table_of_marks_name_in_GAP(tblid::String)
  info = GAPWrap.LibInfoCharacterTable(GapObj(tblid))
  info === GAP.Globals.fail && return false, ""

  tblid = string(info.firstName::GapObj)
  info = GAP.Globals.LIBLIST.TOM_TBL_INFO::GapObj
  pos = findfirst(isequal(GapObj(lowercase(tblid))), Vector{GapObj}(info[2]))
  pos === nothing && return false, ""

  return true, string(info[1][pos]::GapObj)
end

# If `tomid` is the GAP name for a library table of marks
# then translate it to the GAP name of the corresponding
# library character table.
# (There is always such a character table.)
function _character_table_name_in_GAP(tomid::String)
  info = GAP.Globals.LIBLIST.TOM_TBL_INFO
  pos = findfirst(isequal(GapObj(lowercase(tomid))), Vector{GapObj}(info[1]))
  pos === nothing && return false, ""
  info = GAPWrap.LibInfoCharacterTable(info[2][pos])
  return true, string(info.firstName)
end


# access to field values via functions
GapObj(tom::GAPGroupTableOfMarks) = tom.GAPTable

function group(tom::GAPGroupTableOfMarks)
  @req isdefined(tom, :group) "table of marks stores no group"
  return tom.group
end

# If the table of marks stores an Oscar group `G` then
# an isomorphism from `G` to a GAP group is stored in the field `isomorphism`.
# Different kinds of isomorphisms can be stored,
# depending on the type of the stored group,
# like for character tables.
#
function isomorphism_to_GAP_group(tom::GAPGroupTableOfMarks)
  @req isdefined(tom, :group) "table of marks stores no group"
  @assert isdefined(tom, :isomorphism) "table of marks with group must store also an isomorphism"
  return tom.isomorphism
end


"""
    table_of_marks(G::GAPGroup)

Return the table of marks of the finite group `G`.

# Examples
```jldoctest
julia> show(stdout, MIME("text/plain"), table_of_marks(symmetric_group(3)))
Table of marks of symmetric group of degree 3 and order 6

1: 6
2: 3 1
3: 2 . 2
4: 1 1 1 1
```
"""
@attr GAPGroupTableOfMarks function table_of_marks(G::Union{GAPGroup, FinGenAbGroup})
  iso = isomorphism_to_GAP_group(G)
  gaptom = GAPWrap.TableOfMarks(codomain(iso))
  return GAPGroupTableOfMarks(G, gaptom, iso)
end

# A table of marks constructed from a group is stored in this group.
# Tables of marks from the library are cached in the dictionary
# `tables_of_marks_by_id`,
# in order to achieve that fetching the same table twice yields the same
# object.
const tables_of_marks_by_id = Dict{String, Union{GAPGroupTableOfMarks, Nothing}}()

"""
    table_of_marks(id::String)

Return the table of marks
for which `id` is an admissible name in GAP's library of tables of marks.
If no such table is available then `nothing` is returned.

# Examples
```jldoctest
julia> println(table_of_marks("A5"))
table of marks of A5

julia> println(table_of_marks("J5"))
nothing
```
"""
function table_of_marks(id::String)
  flag, tomid =_table_of_marks_name_in_GAP(id)
  flag || return nothing
  return get!(tables_of_marks_by_id, tomid) do
    tom = GAPWrap.TableOfMarks(GapObj(tomid))
    @assert tom !== GAP.Globals.fail "no table of marks with name $tomid?"
    if GAP.Globals.HasUnderlyingGroup(tom)
      G = _oscar_group(GAPWrap.UnderlyingGroup(tom))
      iso = isomorphism_to_GAP_group(G)
      return GAPGroupTableOfMarks(G, tom, iso)
    else
      return GAPGroupTableOfMarks(tom)
    end
  end::GAPGroupTableOfMarks
end

"""
    all_table_of_marks_names(L...; ordered_by = nothing)

Return a vector of strings that contains an admissible name of each
table of marks in the library of tables of marks that satisfies the
conditions in the vector `L`.

The supported conditions in `L` are the same as for
[`all_character_table_names`](@ref),
and the returned vector contains the subset of those names returned by
`all_character_table_names` with the same input for which a table of marks
is available in the library.

# Examples
```
julia> spor_names = all_table_of_marks_names(is_sporadic_simple => true);

julia> println(spor_names[1:5])
["Co3", "HS", "He", "J1", "J2"]

julia> spor_names = all_table_of_marks_names(is_sporadic_simple;
         ordered_by = order);

julia> println(spor_names[1:5])
["M11", "M12", "J1", "M22", "J2"]

julia> length(all_table_of_marks_names(number_of_conjugacy_classes => 5))
4
```
"""
function all_table_of_marks_names(L...; ordered_by = nothing)
  names = all_character_table_names(L...; ordered_by = ordered_by)
  return filter(x -> _table_of_marks_name_in_GAP(x)[1], names)
end
#TODO: As soon as the character table library provides an *a priori*
#      filtering of names for which a table of marks is available,
#      use this instead of filtering afterwards.


"""
    is_table_of_marks_name(name::String)

Return `true` if `table_of_marks(name)` returns a table of marks,
and `false` otherwise

# Examples
```jldoctest
julia> is_table_of_marks_name("J1")
true

julia> is_table_of_marks_name("J4")
false
```
"""
is_table_of_marks_name(name::String) = _table_of_marks_name_in_GAP(name)[1]


# For tables of marks with stored group, we take the hash value of the group.
# For tables of marks without stored group, we take the table identifier.
function Base.hash(tom::GAPGroupTableOfMarks, h::UInt)
  if isdefined(tom, :group)
    return Base.hash(group(tom), h)
  else
    return Base.hash(identifier(tom), h)
  end
end


##############################################################################
#
# `print` and `show` tables of marks

function matrix_of_strings(tom::GAPGroupTableOfMarks)
  n = nrows(tom)
  m = Array{String}(undef, n, n)
  for j in 1:n
    # The matrix is lower triangular.
    for i in 1:(j-1)
      m[i, j] = ""
    end
    for i in j:n
      val = tom[i,j]
      m[i, j] = iszero(val) ?  "." : string(val)
    end
  end
  return m
end

function Base.show(io::IO, tom::GAPGroupTableOfMarks)
  if is_terse(io)
    print(io, "table of marks of a group")
  else
    id = identifier(tom)
    if id != ""
      print(io, "table of marks of ", identifier(tom))
    elseif isdefined(tom, :group)
      print(io, "table of marks of ")
      io = pretty(io)
      print(terse(io), Lowercase(), group(tom))
    else
      print(io, "table of marks of a group")
    end
  end
end

# Produce LaTeX output if `"text/latex"` is prescribed,
# via the `:TeX` attribute of the io context.
function Base.show(io::IO, ::MIME"text/latex", tom::GAPGroupTableOfMarks)
  show(IOContext(io, :TeX => true), MIME("text/plain"), tom)
end

@doc raw"""
    Base.show(io::IO, ::MIME"text/plain", tom::GAPGroupTableOfMarks)

Display the marks of `tom` and context information
as a two-dimensional array.

- First a *header* is shown.
  If `tom` stores a group then the header describes this group,
  otherwise it is equal to the
  [`identifier(tom::GAPGroupTableOfMarks)`](@ref) value of `tom`.

- Then the matrix of marks of `tom` is shown in column portions
  that fit on the screen,
  together with *row labels* (the number `i` for the `i`-th row)
  on the left of each portion.

Output in ``\LaTeX`` syntax can be created by calling `show`
with second argument `MIME("text/latex")`.

# Examples
```jldoctest
julia> tom = table_of_marks("A5");

julia> show(stdout, MIME("text/plain"), tom)
A5

1: 60                
2: 30 2              
3: 20 . 2            
4: 15 3 . 3          
5: 12 . . . 2        
6: 10 2 1 . . 1      
7:  6 2 . . 1 . 1    
8:  5 1 2 1 . . . 1  
9:  1 1 1 1 1 1 1 1 1
```
"""
function Base.show(io::IO, ::MIME"text/plain", tom::GAPGroupTableOfMarks)
  n = nrows(tom)
  mat = matrix_of_strings(tom)
  id = identifier(tom)
  if id != ""
    headerstring = id
  elseif isdefined(tom, :group)
    str_io = IOBuffer()
    print(pretty(str_io), Lowercase(), group(tom))
    headerstring = String(take!(str_io))
    headerstring = "Table of marks of $(headerstring)"
  else
    headerstring = "Table of marks"
  end

  ioc = IOContext(io,
    :header => [headerstring, ""],
    :labels_row => ["$i:" for i in 1:n],
  )

  # print the table
  labeled_matrix_formatted(ioc, mat)
end


##############################################################################
#
length(tom::GAPGroupTableOfMarks) = length(GAPWrap.MarksTom(GapObj(tom)))
number_of_rows(tom::GAPGroupTableOfMarks) = length(GAPWrap.MarksTom(GapObj(tom)))
number_of_columns(tom::GAPGroupTableOfMarks) = length(GAPWrap.MarksTom(GapObj(tom)))

matrix(tom::GAPGroupTableOfMarks) = matrix(ZZ, GAPWrap.MatTom(GapObj(tom)))

@doc raw"""
    order(::Type{T} = ZZRingElem, tom::GAPGroupTableOfMarks) where T <: IntegerUnion

Return the order of the group for which `tom` is the table of marks,
as an instance of `T`.

# Examples
```jldoctest
julia> order(table_of_marks(symmetric_group(4)))
24
```
"""
order(tom::GAPGroupTableOfMarks) = order(ZZRingElem, tom)

function order(::Type{T}, tom::GAPGroupTableOfMarks) where T <: IntegerUnion
  return T(GAPWrap.MarksTom(GapObj(tom))[1][1])
end

@doc raw"""
    orders_class_representatives(tom::GAPGroupTableOfMarks)

Return the vector of the orders of conjugacy class representatives for `tom`,
ordered according to the rows and columns of `tom`.

# Examples
```jldoctest
julia> println(orders_class_representatives(table_of_marks("A5")))
ZZRingElem[1, 2, 3, 4, 5, 6, 10, 12, 60]
```
"""
@gapattribute orders_class_representatives(tom::GAPGroupTableOfMarks) = Vector{ZZRingElem}(GAP.Globals.OrdersTom(GapObj(tom)))


@doc raw"""
    class_lengths(tom::GAPGroupTableOfMarks)

Return the vector of the lengths of the conjugacy classes of subgroups
for `tom`.

# Examples
```jldoctest
julia> println(class_lengths(table_of_marks("A5")))
ZZRingElem[1, 15, 10, 5, 6, 10, 6, 5, 1]
```
"""
@gapattribute class_lengths(tom::GAPGroupTableOfMarks) = Vector{ZZRingElem}(GAP.Globals.LengthsTom(GapObj(tom)))


@doc raw"""
    identifier(tom::GAPGroupTableOfMarks)

Return a string that identifies `tom` if `tom` belongs to the library
of tables of marks, and an empty string otherwise.

# Examples
```jldoctest
julia> identifier(table_of_marks("A5"))
"A5"

julia> identifier(table_of_marks(symmetric_group(3)))
""
```

!!! warning "Note:"
    The `identifier` of a table of marks from the library is equal to
    the `identifier` of the corresponding library character table.
    In a few cases, this value differs from the `GAP.Globals.Identifier`
    value of the underlying `GapObj` of  the table of marks.
"""
@attr String function identifier(tom::GAPGroupTableOfMarks)
  gaptom = GapObj(tom)
  GAPWrap.IsLibTomRep(gaptom) || return ""
  flag, id = _character_table_name_in_GAP(string(GAPWrap.Identifier(gaptom)))
  @assert flag "no character table for GAP table of marks?"
  return id
end


function Base.getindex(tom::GAPGroupTableOfMarks, i::Int)
  marks = GAPWrap.MarksTom(GapObj(tom))[i]
  subs = GAPWrap.SubsTom(GapObj(tom))[i]
  v = zeros(ZZRingElem, i)
  for j in 1:length(subs)
    v[subs[j]] = marks[j]
  end
  return marks_vector(tom, v)
end

# in order to make `tom[end]` work
Base.lastindex(tom::GAPGroupTableOfMarks) = length(tom)

# in order to make `findfirst` and `findall` work
function Base.keys(tom::GAPGroupTableOfMarks)
  return keys(1:length(tom))
end

function Base.getindex(tom::GAPGroupTableOfMarks, i::Int, j::Int)
  i < j && return ZZ(0)
  gaptom = GapObj(tom)
  subs = Vector{Int}(GAPWrap.SubsTom(gaptom)[i])
  pos = findfirst(isequal(j), subs)
  pos === nothing && return ZZ(0)
  marks = GAPWrap.MarksTom(gaptom)[i]
  return ZZ(marks[pos])
end

Base.iterate(tom::GAPGroupTableOfMarks, state = 1) = state > nrows(tom) ? nothing : (tom[state], state+1)


@doc raw"""
    representative(tom::GAPGroupTableOfMarks, i::Int)

Return a representative from the `i`-th class of subgroups of `tom`.

An exception is thrown if `tom` does not store a group,
or if `i` is larger than `length(tom)`.

# Examples
```jldoctest
julia> representative(table_of_marks("A5"), 2)
Permutation group of degree 5 and order 2 with 1 generator
  (2,3)(4,5)
```
"""
function representative(tom::GAPGroupTableOfMarks, i::Int)
  G = group(tom)
  @req 0 < i <= length(tom) "table of marks has only $(length(tom)) columns"
  return _oscar_subgroup(GAPWrap.RepresentativeTom(GapObj(tom), i), G)
end


#############################################################################
##
##  marks vectors
##
##  In order to describe Burnside rings over the integers,
##  we implement marks vectors as wrapped GAP vectors,
##  with parent the table of marks in question.
##  Note that marks vectors can be shorter than the number of columns
##  of the table of marks,
##  meaning that the values at larger positions are zero.
##
abstract type GroupMarksVector end

struct GAPGroupMarksVector <: GroupMarksVector
  table::GAPGroupTableOfMarks
  values::GapObj
end

GAP.@install GapObj(chi::GAPGroupMarksVector) = chi.values

parent(chi::GAPGroupMarksVector) = chi.table

values(chi::GAPGroupMarksVector) = Vector{ZZRingElem}(GapObj(chi))

group(chi::GAPGroupMarksVector) = group(parent(chi))

function Base.show(io::IO, chi::GAPGroupMarksVector)
  print(io, "marks_vector($(parent(chi)), $(values(chi)))")
end

function marks_vector(tom::GAPGroupTableOfMarks, values::GapObj)
  @req GAPWrap.IsHomogeneousList(values) "values must be a list of integers"
  GAPWrap.ShrinkRowVector(values)
  vector = Vector{ZZRingElem}(values)
  @req length(vector) <= ncols(tom) "the length of values must be at most $(ncols(tom))"
  return GAPGroupMarksVector(tom, values)
end

function marks_vector(tom::GAPGroupTableOfMarks, values::Vector{<:IntegerUnion})
  vector = GapObj([GAP.Obj(x) for x in values])
  GAPWrap.ShrinkRowVector(vector)
  return GAPGroupMarksVector(tom, vector)
end

Base.length(chi::GAPGroupMarksVector) = ncols(parent(chi))

Base.iterate(chi::GAPGroupMarksVector, state = 1) = state > length(chi) ? nothing : (chi[state], state+1)

# access values by position
function Base.getindex(chi::GAPGroupMarksVector, i::Int)
  v = GapObj(chi)
  return i <= length(v) ? ZZ(v[i]) : ZZ(0)
end

function Base.getindex(chi::GAPGroupMarksVector, l::Vector{Int})
  v = GapObj(chi)
  res = zeros(ZZRingElem, length(l))
  for i in 1:length(l)
    res[i] = l[i] <= length(v) ? ZZ(v[l[i]]) : ZZ(0)
  end
  return res
end

# arithmetic with marks vectors
function Base.:(==)(chi::GAPGroupMarksVector, psi::GAPGroupMarksVector)
  @req parent(chi) === parent(psi) "tables of marks must be identical"
  return GapObj(chi) == GapObj(psi)
end

# `hash` method based on the values
function Base.hash(chi::GAPGroupMarksVector, h::UInt)
  return Base.hash(values(chi), h)
end

function Base.:+(chi::GAPGroupMarksVector, psi::GAPGroupMarksVector)
  @req parent(chi) === parent(psi) "tables of marks must be identical"
  return marks_vector(parent(chi), GapObj(chi) + GapObj(psi))
end

Base.:-(chi::GAPGroupMarksVector) = marks_vector(parent(chi), -GapObj(chi))

function Base.:-(chi::GAPGroupMarksVector, psi::GAPGroupMarksVector)
  @req parent(chi) === parent(psi) "tables of marks must be identical"
  return marks_vector(parent(chi), GapObj(chi) - GapObj(psi))
end

function Base.:*(chi::GAPGroupMarksVector, psi::GAPGroupMarksVector)
  @req parent(chi) === parent(psi) "tables of marks must be identical"
  vchi = GapObj(chi)
  vpsi = GapObj(psi)
  n = min(length(vchi), length(vpsi))
  vector = GapObj([vchi[i] * vpsi[i] for i in 1:n])
  return marks_vector(parent(chi), vector)
end

function Base.zero(chi::GAPGroupMarksVector)
  return marks_vector(parent(chi), GapObj([]))
end

function Base.one(chi::GAPGroupMarksVector)
  val = ZZ(1)
  return marks_vector(parent(chi), [val for i in 1:length(chi)])
end

function Base.:*(n::IntegerUnion, chi::GAPGroupMarksVector)
  return marks_vector(parent(chi), n * values(chi))
end

function Base.:^(chi::GAPGroupMarksVector, n::IntegerUnion)
  vchi = GapObj(chi)
  return marks_vector(parent(chi), GapObj([vchi[i]^n for i in 1:length(vchi)]))
end

@doc raw"""
    coordinates(::Type{T} = ZZRingElem, chi::GAPGroupMarksVector)
                   where T <: Union{IntegerUnion, QQFieldElem}

Return the vector $[a_1, a_2, \ldots, a_n]$ such that `chi`
is equal to $\sum_{i=1}^n a_i t[i]$ where $t$ is `parent(chi)`.

The result is an instance of `Vector{T}`.
Note that the result can be shorter than `ncols(parent(chi))`.

# Examples
```jldoctest
julia> tom = table_of_marks(symmetric_group(4));

julia> chi = tom[3] * tom[6]
marks_vector(table of marks of Sym(4), ZZRingElem[72, 0, 4])

julia> println(coordinates(Int, chi))
[2, 0, 2]
```
"""
coordinates(chi::GAPGroupMarksVector) = coordinates(ZZRingElem, chi)

function coordinates(::Type{T}, chi::GAPGroupMarksVector) where T <: Union{Integer, ZZRingElem, QQFieldElem}
  c = GAPWrap.DecomposedFixedPointVector(GapObj(parent(chi)), GapObj(chi))
  return Vector{T}(c)::Vector{T}
end


#############################################################################
##
##  interface between tables of marks and character tables
##
##  A table of marks and the corresponding character table store each other
##  in attributes.
##

@doc raw"""
    character_table(tom::GAPGroupTableOfMarks)

Return the character table of the group of `tom`.
If `tom` belongs to the library of tables of marks then the corresponding
character table from the library of character tables is returned,
otherwise the character table of `group(tom)`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  tom = table_of_marks(g);

julia> character_table(tom) == character_table(g)
true
```
"""
@attr GAPGroupCharacterTable function character_table(tom::GAPGroupTableOfMarks)
  GAPWrap.IsLibTomRep(GapObj(tom)) && return character_table(identifier(tom))
  tbl = character_table(group(tom))
  set_attribute!(tbl, :table_of_marks, tom)
  return tbl
end

@doc raw"""
    table_of_marks(tbl::GAPGroupCharacterTable)

Return the table of marks of the group of `tbl`.
If `tbl` does not store a group and if the library of tables of marks
contains the table of marks corresponding to `tbl` then this is returned,
otherwise `nothing`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  tbl = character_table(g);

julia> table_of_marks(tbl) == table_of_marks(g)
true
```
"""
@attr Union{Nothing, GAPGroupTableOfMarks} function table_of_marks(tbl::GAPGroupCharacterTable)
  isdefined(tbl, :group) && return table_of_marks(group(tbl))
  tom = table_of_marks(identifier(tbl))
  tom === nothing && return nothing
  set_attribute!(tom, :character_table, tbl)
  return tom
end


@doc raw"""
    restrict(chi::GAPGroupMarksVector, tbl::GAPGroupCharacterTable)

Return the class function with parent `tbl` that is the restriction of `chi`.
For that, `parent(chi)` and `tbl` must belong to the same group ``G``.

If `chi` is the `i`-th row in the table of marks `parent(chi)`
then the result is the permutation character of the action of ``G``
on the right cosets of its subgroup `representative(parent(chi), i)`.

# Examples
```jldoctest
julia> tom = table_of_marks("A5");  tbl = character_table(tom);

julia> chi = tom[5]
marks_vector(table of marks of A5, ZZRingElem[12, 0, 0, 0, 2])

julia> println(values(restrict(chi, tbl)))
QQAbFieldElem{AbsSimpleNumFieldElem}[12, 0, 0, 2, 2]
```
"""
function restrict(chi::GAPGroupMarksVector, tbl::GAPGroupCharacterTable)
  tom = parent(chi)
  p = characteristic(tbl)
  ordtbl = p == 0 ? tbl : ordinary_table(tbl)
  fus = Vector{Int}(GAPWrap.FusionCharTableTom(GapObj(ordtbl), GapObj(tom)))
  if p != 0
    flag, modfus = known_class_fusion(tbl, ordtbl)
    @assert flag
    fus = fus[modfus]
  end
  return class_function(tbl, chi[fus])
end
