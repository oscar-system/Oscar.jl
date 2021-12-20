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

import Base: getindex, length, mod, one, print, show, zero

import Oscar.AbelianClosure: QabElem, QabAutomorphism

import Nemo: degree

export
    character_field,
    character_table,
    decomposition_matrix,
    scalar_product,
    trivial_character


complex_conjugate(elm::QabElem) = elm^QabAutomorphism(-1)


#############################################################################
##
##  character tables
##
abstract type GroupCharacterTable end

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
Sym( [ 1 .. 3 ] )

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
    tbls = get_attribute(G, :character_tables)
    if tbls == nothing
      tbls = Dict()
      set_attribute!(G, :character_tables => tbls)
    end

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
            valbar = complex_conjugate(val)
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
    size = fmpz(GAP.Globals.Size(gaptbl))
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

    modtbls = get_attribute(tbl, :brauer_tables)
    if modtbls == nothing
      modtbls = Dict{Int,Any}()
      set_attribute!(tbl, :brauer_tables => modtbls)
    end
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


#############################################################################
##
##  class functions (and characters)
##
abstract type GroupClassFunction end
#TODO: support character rings and elements of it?
#      if yes then there is no need to have class other functions than these

struct GAPGroupClassFunction <: GroupClassFunction
    table::GAPGroupCharacterTable
    values::GAP.GapObj
end

function Base.show(io::IO, chi::GAPGroupClassFunction)
    print(io, "group_class_function($(chi.table), $(values(chi)))")
end

import Base.values

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

function trivial_character(tbl::GAPGroupCharacterTable)
    val = QabElem(1)
    return group_class_function(tbl, [val for i in 1:ncols(tbl)])
end

function trivial_character(G::GAPGroup)
    val = QabElem(1)
    return group_class_function(G, [val for i in 1:GAP.Globals.NrConjugacyClasses(G.X)])
end

Base.length(chi::GAPGroupClassFunction) = length(chi.values)

Base.iterate(chi::GAPGroupClassFunction, state = 1) = state > length(chi.values) ? nothing : (chi[state], state+1)

# the degree is an fmpq
# (for general class functions, denominators can occur)
function Nemo.degree(chi::GAPGroupClassFunction)
    val = values(chi)[1]
    return Nemo.coeff(val.data, 0)
end

# access character values
Base.getindex(chi::GAPGroupClassFunction, i::Int) = QabElem(GAP.Globals.ValuesOfClassFunction(chi.values)[i])

function Base.:(==)(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    chi.table === psi.table || error("character tables must be identical")
#T check_parent?
    return chi.values == psi.values
end

# aritmetics with class functions
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

function scalar_product(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    chi.table === psi.table || error("character tables must be identical")
    return Nemo.fmpz(GAP.Globals.ScalarProduct(chi.values, psi.values))
end

function Base.:*(n::IntegerUnion, chi::GAPGroupClassFunction)
    return GAPGroupClassFunction(chi.table, n * chi.values)
end

function Base.:^(chi::GAPGroupClassFunction, n::IntegerUnion)
    return GAPGroupClassFunction(chi.table, chi.values ^ n)
end

# apply a class function to a group element
function(chi::GAPGroupClassFunction)(g::BasicGAPGroupElem)
    # Identify the conjugacy class of `g`.
    ccl = GAP.Globals.ConjugacyClasses(GAP.Globals.UnderlyingGroup(chi.table.GAPTable))
    for i in 1:length(ccl)
      if g.X in ccl[i]
        return chi[i]
      end
    end
    error("$g is not an element in the underlying group")
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
      c = matrix(QQ, powers)'

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
        b = solve(c, matrix(QQ,length(a),1,a))'
        b = [b[i] for i in 1:length(b)]
        return F(b)
      end
    end

    return F, MapFromFunc(f, finv, F, FF)
end
