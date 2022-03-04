##  This is a first attempt to implement group characters in Oscar.
##  
##  The idea is that the available GAP objects (groups, character tables,
##  class functions) are used in a first step, and that access to character
##  values yields `QabElem` objects.
##  
##  Once we agree on the functionality and the integration into Oscar,
##  this setup can in a second step be replaced by one that uses
##  native Julia objects for representing class functions,
##  but character tables and groups still have some counterpart in GAP.
##  
##  In a third step, we replace the character table objects by native Julia
##  objects.

# character values are elements from QabField

export
    atlas_irrationality,
    character_field,
    character_table,
    decomposition_matrix,
    indicator,
    induced_cyclic,
    known_class_fusion,
    natural_character,
    scalar_product,
    schur_index,
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
if `F` is not given then the result has type `QabElem`.

`description` is assumed to have the format defined in
[CCNPW85](@cite), Chapter 6, Section 10.

```jldoctest
julia> atlas_irrationality("r5")
-2*ζ(5)^3 - 2*ζ(5)^2 - 1

julia> atlas_irrationality(CyclotomicField(5)[1], "r5")
-2*z_5^3 - 2*z_5^2 - 1

julia> atlas_irrationality("i")
ζ(4)

julia> atlas_irrationality("b7*3")
-ζ(7)^4 - ζ(7)^2 - ζ(7) - 1

julia> atlas_irrationality("3y'''24*13-2&5")
-5*ζ(24)^7 - 2*ζ(24)^5 + 2*ζ(24)^3 - 3*ζ(24)

```
"""
function atlas_irrationality(F::AnticNumberField, description::String)
    return F(GAP.Globals.AtlasIrrationality(GAP.GapObj(description)))
end

function atlas_irrationality(description::String)
    F = abelian_closure(QQ)[1]
    return F(GAP.Globals.AtlasIrrationality(GAP.GapObj(description)))
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
    GAPGroup::GAPGroup    # the underlying group, if any
    GAPTable::GAP.GapObj  # the character table object
    characteristic::Int

    function GAPGroupCharacterTable(G::GAPGroup, tab::GAP.GapObj, char::Int)
      ct = new()
      ct.GAPGroup = G
      ct.GAPTable = tab
      ct.characteristic = char
      return ct
    end

    function GAPGroupCharacterTable(tab::GAP.GapObj, char::Int)
      ct = new()
      #ct.GAPGroup is left undefined
      ct.GAPTable = tab
      ct.characteristic = char
      return ct
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
julia> character_table( symmetric_group(3) )
Sym( [ 1 .. 3 ] )

 2  1  1  .
 3  1  .  1
           
   1a 2a 3a
2P 1a 1a 3a
3P 1a 2a 1a
           
χ₁  1 -1  1
χ₂  2  . -1
χ₃  1  1  1


julia> character_table( symmetric_group(3), 2 )
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
      gaptbl = GAP.Globals.CharacterTable(G.X)
      if p != 0
        # Create the `p`-modular table if possible.
        isprime(p) || error("p must be 0 or a prime integer")
        gaptbl = GAP.Globals.mod(gaptbl, GAP.Obj(p))
        gaptbl == GAP.Globals.fail && return nothing
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
      isprime(p) || error("p must be 0 or a prime integer")
      modid = "$(id)mod$(p)"
    end

    return get!(character_tables_by_id, modid) do
      tbl = GAP.Globals.CharacterTable(GAP.GapObj(modid))
      tbl == GAP.Globals.fail && return nothing
      return GAPGroupCharacterTable(tbl, p)
    end
end

##############################################################################
#
# admissible names of library character tables

function all_character_table_names()
    K = GAP.Globals.CallFuncList(GAP.Globals.AllCharacterTableNames, GAP.GapObj([]))
    return Vector{String}(K)
end
#TODO:
# Support function/value pairs as arguments, similar to (but more general
# than) `all_small_groups` etc.
# This makes sense only if GAP's Browse package is available (and has been
# loaded at the time when the character table library got loaded),
# otherwise everything is too slow.
# Currently this cannot be assumed.


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
    flag, N = Hecke.iscyclotomic_type(F)
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
    cents = Vector{fmpz}(GAP.Globals.SizesCentralizers(gaptbl))
    fcents = [collect(factor(x)) for x in cents]
    d = Dict([p => fill(".", n) for p in primes]...)
    for i in 1:n
      for pair in fcents[i]
        d[pair[1]][i] = string(pair[2])
      end
    end
    cents_strings = [d[p] for p in primes]

    # Compute display format for power maps.
    names = Vector{String}(GAP.Globals.ClassNames(gaptbl))
    pmaps = Vector{Any}(GAP.Globals.ComputedPowerMaps(gaptbl))
    power_maps_primes = String[]
    power_maps_strings = Vector{String}[]
    for i in 2:length(pmaps)
      map = pmaps[i]
      if map != nothing
        push!(power_maps_primes, string(i)*"P")
        push!(power_maps_strings, names[map])
      end
    end

    empty = ["" for i in 1:n]

    if isdefined(tbl, :GAPGroup)
      headerstring = string(tbl.GAPGroup)
      if tbl.characteristic != 0
        headerstring = "$headerstring mod $(tbl.characteristic)"
      end
    else
      headerstring = String(GAP.Globals.Identifier(gaptbl))
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
        cents_strings..., empty, names, power_maps_strings..., empty)),

      # row labels:
      # character names (a column vector is sufficient)
      :labels_row => ["\\chi_{$i}" for i in 1:n],

      # corner (a column vector is sufficient):
      # primes in the centralizer rows,
      # separating empty line,
      # separating empty line,
      # primes in the power map rows,
      # separating empty line,
      :corner => vcat( string.(primes), ["", ""], power_maps_primes, [""] ),

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
      id = "\"" * String(GAP.Globals.Identifier(gaptbl)) * "\""
    end
    print(io, "character_table($id)")
end


##############################################################################
#
length(tbl::GAPGroupCharacterTable) = GAP.Globals.NrConjugacyClasses(tbl.GAPTable)
Oscar.nrows(tbl::GAPGroupCharacterTable) = GAP.Globals.NrConjugacyClasses(tbl.GAPTable)
Oscar.ncols(tbl::GAPGroupCharacterTable) = GAP.Globals.NrConjugacyClasses(tbl.GAPTable)

function Base.getindex(tbl::GAPGroupCharacterTable, i::Int)
    return group_class_function(tbl, GAP.Globals.Irr(tbl.GAPTable)[i])
end
#TODO: cache the irreducibles in the table

function Base.getindex(tbl::GAPGroupCharacterTable, i::Int, j::Int)
    val = GAP.Globals.Irr(tbl.GAPTable)[i, j]
    return QabElem(val)
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
    isprime(p) || error("p must be a prime integer")
    tbl.characteristic == 0 || error("tbl mod p only for ordinary table tbl")

    modtbls = get_attribute!(() -> Dict{Int,Any}(), tbl, :brauer_tables)
    if ! haskey(modtbls, p)
      modtblgap = mod(tbl.GAPTable, p)
      if modtblgap == GAP.Globals.fail
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
    isprime(modtbl.characteristic) || error("characteristic of tbl must be a prime integer")
    return matrix(ZZ, GAP.Globals.DecompositionMatrix(modtbl.GAPTable))
end

@doc Markdown.doc"""
    names_of_fusion_sources(tbl::GAPGroupCharacterTable)

Return the array of strings that are identifiers of those character tables
which store a class fusion to `tbl`.
"""
function names_of_fusion_sources(tbl::GAPGroupCharacterTable)
    return [string(name) for name in GAP.Globals.NamesOfFusionSources(tbl.GAPTable)]
end

@doc Markdown.doc"""
    known_class_fusion(tbl1::GAPGroupCharacterTable, tbl2::GAPGroupCharacterTable)

Return `(flag, fus)` where `flag == true` if a class fusion to `tbl2` is stored
on `tbl1`, and `flag == false` otherwise.

In the former case,
`fus` is the vector of integers, of length `number_conjugacy_classes(tbl1)`,
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
    map = GAP.Globals.GetFusionMap(subtbl.GAPTable, tbl.GAPTable)
    if map == GAP.Globals.fail
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
    values::GAP.GapObj
end

function Base.show(io::IO, chi::GAPGroupClassFunction)
    print(io, "group_class_function($(chi.table), $(values(chi)))")
end

function values(chi::GAPGroupClassFunction)
    gapvalues = GAP.Globals.ValuesOfClassFunction(chi.values)
    return [QabElem(x) for x in gapvalues]
end

function group_class_function(tbl::GAPGroupCharacterTable, values::GAP.GapObj)
    GAPWrap.IsClassFunction(values) || error("values must be a class function")
    return GAPGroupClassFunction(tbl, values)
end

function group_class_function(tbl::GAPGroupCharacterTable, values::Vector{QabElem})
    gapvalues = GAP.GapObj([GAP.Obj(x) for x in values])
    return GAPGroupClassFunction(tbl, GAP.Globals.ClassFunction(tbl.GAPTable, gapvalues))
end

function group_class_function(G::GAPGroup, values::Vector{QabElem})
    return group_class_function(character_table(G), values)
end

@doc Markdown.doc"""
    trivial_character(tbl::GAPGroupCharacterTable)

Return the character of `tbl` that has the value `QabElem(1)` in each position.
"""
function trivial_character(tbl::GAPGroupCharacterTable)
    val = QabElem(1)
    return group_class_function(tbl, [val for i in 1:ncols(tbl)])
end

@doc Markdown.doc"""
    trivial_character(G::GAPGroup)

Return the character of (the ordinary character table of) `G`
that has the value `QabElem(1)` in each position.
"""
function trivial_character(G::GAPGroup)
    val = QabElem(1)
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
    induced_cyclic(tbl::GAPGroupCharacterTable)

Return the array of permutation characters of `tbl` that are induced from
cyclic subgroups.
"""
function induced_cyclic(tbl::GAPGroupCharacterTable)
    return [GAPGroupClassFunction(tbl, chi) for chi in GAP.Globals.InducedCyclic(tbl.GAPTable)]
end

Base.length(chi::GAPGroupClassFunction) = length(chi.values)

Base.iterate(chi::GAPGroupClassFunction, state = 1) = state > length(chi.values) ? nothing : (chi[state], state+1)

@doc Markdown.doc"""
    degree(::Type{T} = fmpq, chi::GAPGroupClassFunction)
           where T <: Union{IntegerUnion, fmpz, mpq, QabElem}

Return `chi[1]`, as an instance of `T`.
"""
Nemo.degree(chi::GAPGroupClassFunction) = Nemo.degree(fmpq, chi)::fmpq

Nemo.degree(::Type{fmpq}, chi::GAPGroupClassFunction) = Nemo.coeff(values(chi)[1].data, 0)::fmpq

Nemo.degree(::Type{fmpz}, chi::GAPGroupClassFunction) = ZZ(Nemo.coeff(values(chi)[1].data, 0))::fmpz

Nemo.degree(::Type{QabElem}, chi::GAPGroupClassFunction) = values(chi)[1]::QabElem

Nemo.degree(::Type{T}, chi::GAPGroupClassFunction) where T <: IntegerUnion = T(Nemo.degree(fmpz, chi))::T

# access character values
Base.getindex(chi::GAPGroupClassFunction, i::Int) = QabElem(GAP.Globals.ValuesOfClassFunction(chi.values)[i])

# arithmetics with class functions
function Base.:(==)(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    chi.table === psi.table || error("character tables must be identical")
#T check_parent?
    return chi.values == psi.values
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
    val = QabElem(0)
    return group_class_function(chi.table, [val for i in 1:length(chi)])
end

Base.one(chi::GAPGroupClassFunction) = trivial_character(chi.table)

@doc Markdown.doc"""
    scalar_product(::Type{T} = fmpq, chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
                   where T <: Union{IntegerUnion, fmpz, fmpq, QabElem}

Return $\sum_{g \in G}$ `chi`($g$) `conj(psi)`($g$) / $|G|$,
where $G$ is the group of both `chi` and `psi`.
"""
scalar_product(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction) = scalar_product(fmpq, chi, psi)

function scalar_product(::Type{T}, chi::GAPGroupClassFunction, psi::GAPGroupClassFunction) where T <: Union{Integer, fmpz, fmpq, QabElem}
    chi.table === psi.table || error("character tables must be identical")
    return T(GAP.Globals.ScalarProduct(chi.values, psi.values))::T
end

function Base.:*(n::IntegerUnion, chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(chi.table, n * chi.values)
end

function Base.:^(chi::GAPGroupClassFunction, n::IntegerUnion)
    return GAPGroupClassFunction(chi.table, chi.values ^ n)
end

function Base.:^(chi::GAPGroupClassFunction, tbl::GAPGroupCharacterTable)
    return GAPGroupClassFunction(tbl, GAP.Globals.InducedClassFunction(chi.values, tbl.GAPTable))
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
    return GAPGroupClassFunction(chi.table, GAP.Globals.GaloisCyc(chi.values, -1))
end

@doc Markdown.doc"""
    (sigma::QabAutomorphism)(chi::GAPGroupClassFunction)

Return the class function whose values are the images of the values of `chi`
under `sigma`.
"""
function (sigma::QabAutomorphism)(chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(chi.table, GAP.Globals.GaloisCyc(chi.values, sigma.exp))
end

Base.:^(chi::Oscar.GAPGroupClassFunction, sigma::QabAutomorphism) = sigma(chi)

@doc Markdown.doc"""
    isirreducible(chi::GAPGroupClassFunction)

Return `true` if `chi` is an irreducible character, and `alse` otherwise.

A character is irreducible if it cannot be written as the sum of two
characters.
For ordinary characters this can be checked using the scalar product of
class functions (see [`scalar_product`](@ref).
For Brauer characters there is no generic method for checking irreducibility.
"""
function isirreducible(chi::GAPGroupClassFunction)
    return GAP.Globals.IsIrreducibleCharacter(chi.table, chi.values)
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
    return GAP.Globals.Indicator(chi.table.GAPTable, GAP.GapObj([chi.values]), n)[1]::Int
end

@doc Markdown.doc"""
    character_field(chi::GAPGroupClassFunction)

Return the pair `(F, phi)` where `F` is a number field that is generated
by the character values of `chi`, and `phi` is the embedding of `F` into
`abelian_closure(QQ)`.
"""
function character_field(chi::GAPGroupClassFunction)
    values = chi.values  # a list of GAP cyclotomics
    gapfield = GAP.Globals.Field(values)
    N = GAP.Globals.Conductor(gapfield)
    FF, = abelian_closure(QQ)
    if GAP.Globals.IsCyclotomicField(gapfield)
      # In this case, the want to return a field that knows to be cyclotomic
      # (and the embedding is easy).
      F, z = Oscar.AbelianClosure.cyclotomic_field(FF, N)
      f = x::nf_elem -> QabElem(x, N)
      finv = function(x::QabElem)
        g = gcd(x.c, N)
        K, = Oscar.AbelianClosure.cyclotomic_field(FF, g)
        x = Hecke.force_coerce_cyclo(K, x.data)
        x = Hecke.force_coerce_cyclo(F, x)
        return x
      end
    else
      # In the general case, we have to work for the embedding.
      gapgens = GAP.Globals.GeneratorsOfField(gapfield)
      @assert length(gapgens) == 1
      gappol = GAP.Globals.MinimalPolynomial(GAP.Globals.Rationals, gapgens[1])
      gapcoeffs = GAP.Globals.CoefficientsOfUnivariatePolynomial(gappol)
      coeffscyc = Vector{fmpq}(GAP.Globals.COEFFS_CYC(gapgens[1]))
      v = Vector{fmpq}(gapcoeffs)
      R, = PolynomialRing(QQ, "x")
      f = R(v)
      F, z = NumberField(f, "z"; cached = true, check = false)
      K, zz = Oscar.AbelianClosure.cyclotomic_field(FF, N)

      nfelm = QabElem(gapgens[1]).data

      # Compute the expression of powers of `z` as sums of roots of unity (once).
      powers = [coefficients(Hecke.force_coerce_cyclo(K, nfelm^i)) for i in 0:length(v)-2]
      c = transpose(matrix(QQ, powers))

      f = function(x::nf_elem)
        return QabElem(evaluate(R(x), nfelm), N)
      end

      finv = function(x::QabElem)
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

Return either the minimal integer `m` such that the character `m * chi`
is afforded by a representation over the character field of `chi`,
or `nothing`.

The latter happens if character theoretic criteria do not suffice for
computing `m`.
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
      gapfield = GAP.Globals.Field(values)
      N = GAP.Globals.Conductor(gapfield)
      for n in reverse(sort(divisors(N)))
        if GAP.Globals.E(n) in gapfield
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
      known, fus = known_class_fusion(s, tbl)
      @assert known "the class fusion is not stored"
      if length(class_positions_of_kernel(fus)) == 1
        psi = trivial_character(s)^(tbl)
        bound = gcd(bound, scalar_product(fmpz, chi, psi))
        bound == 1 && return 1
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

    # For the moment, we do not have more character theoretic criteria.
    return nothing
end
