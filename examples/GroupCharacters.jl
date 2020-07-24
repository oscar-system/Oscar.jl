"""
This is a first attempt to implement group characters in Oscar.

The idea is that the available GAP objects (groups, character tables,
class functions) are used in a first step, and that access to character
values yields `QabElem` objects.

Once we agree on the functionality and the integration into Oscar,
this setup can in a second step be replaced by one that uses
native Julia objects for representing class functions,
but character tables and groups still have some counterpart in GAP.

In a third step, we replace the character table objects by native Julia
objects.

And here are a few examples how the current code is intended.

# Examples
```
include( "GroupCharacters.jl" )

using Oscar

# an ordinary character table from the GAP library
t = Main.GroupCharacters.character_table("A5");
print(t)   # short form
show(t)    # call GAP's 'Display'

chi = t[2]  # a character
t[2,4]  # a character value
chi[4]  # a character value

length(t) # number of conjugacy classes
Main.GroupCharacters.trivial_character(t)

values(chi)
chi == t[1]
chi == t[2]
chi + chi
chi * chi
- chi
chi - chi
2 * chi
chi ^ 2
one(chi)
zero(chi)
degree(chi)
length(chi)

Main.GroupCharacters.scalar_product(t[1], t[1])
Main.GroupCharacters.scalar_product(t[1], t[2])

# a Brauer character table
tmod2 = mod(t, 2);
print(tmod2)
show(tmod2)

tmod2[2]  # a character
tmod2[2,3]  # a character value

length(tmod2) # number of conjugacy classes
Main.GroupCharacters.trivial_character(tmod2)

Main.GroupCharacters.decomposition_matrix(tmod2)

# an ordinary character table computed from a group
G = Oscar.symmetric_group(4)
t = Main.GroupCharacters.character_table(G);
print(t)
show(t)

t[2]
t[2,3]

# the corresponding Brauer table (possible because the group is solvable)
tmod2 = mod(t, 2);
print(tmod2)
show(tmod2)

# in general, GAP cannot compute Brauer tables from groups
G = Oscar.symmetric_group(5)
t = Main.GroupCharacters.character_table(G);
mod(t, 2)
```
"""
module GroupCharacters

using Oscar

# character values are elements from QabField
#using QabModule
#import QabModule: QabElem
include("QabAndPChars.jl")

# functionality for displaying character tables
include("MatrixDisplay.jl")

import Base: getindex, length, mod, one, print, show, zero

import AbstractAlgebra: nrows, ncols

import Nemo: degree

# Load GAP's character table library (if it is installed)
GAP.Packages.load("ctbllib")


#############################################################################
##
##  conversion between Julia's QabElem and GAP's cyclotomics
##
function QabModule.QabElem(cyc::Union{GAP.GapObj,Int64})
    GAP.Globals.IsCyc(cyc) || error("cyc must be a GAP cyclotomic")
    denom = GAP.Globals.DenominatorCyc(cyc)
    n = GAP.Globals.Conductor(cyc)
    coeffs = GAP.Globals.ExtRepOfObj(cyc * denom)
    cycpol = GAP.Globals.CyclotomicPol(n)
    GAP.Globals.ReduceCoeffs(coeffs, cycpol)
    coeffs = GAP.gap_to_julia(coeffs, recursive = true)
    coeffs = coeffs[1:(length(cycpol)-1)]
    denom = GAP.gap_to_julia(denom)
    F, z = Nemo.CyclotomicField(n)
    val = Nemo.elem_from_mat_row(F, Nemo.matrix(Nemo.ZZ, [coeffs]), 1, fmpz(denom))
    return QabModule.QabElem(val,n)
end

#function QabModule.QabElem(cyc::GAP.GapObj, c::Int)
#end

function gap_cyclotomic(elm::QabModule.QabElem)
    coeffs = [Nemo.coeff(elm.data, i) for i in 0:(elm.c-1)]  # fmpq
    coeffs = [GAP.julia_to_gap(BigInt(numerator(x)) // BigInt(denominator(x)))
              for x in coeffs]
    return GAP.Globals.CycList(GAP.julia_to_gap(coeffs))
end


#############################################################################
##
##  character tables
##
abstract type GroupCharacterTable end

mutable struct GAPGroupCharacterTable <: GroupCharacterTable
    GAPGroup::Oscar.GAPGroup    # the underlying group, if any
    GAPTable::GAP.GapObj  # the character table object
    characteristic::Int
    AbstractAlgebra.@declare_other

    function GAPGroupCharacterTable(G::Oscar.GAPGroup, tab::GAP.GapObj, char::Int)
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

function character_table(G::Oscar.GAPGroup, p::Int = 0)
    tbl = AbstractAlgebra.get_special(G, :character_table)
    if tbl == nothing
      tbl = GAP.Globals.CharacterTable(G.X)
      AbstractAlgebra.set_special(G, :character_table => tbl)
    end

    if p != 0
      isprime(p) || error("p must be 0 or a prime integer")
      tbl = mod(tbl, p)
    end

    return GAPGroupCharacterTable(G, tbl, p)
end

function character_table(id::String, p::Int = 0)
    if GAP.Globals.IsBoundGlobal(GAP.julia_to_gap("CTblLib"))
      tbl = GAP.Globals.CharacterTable(GAP.julia_to_gap(id))
      if tbl == GAP.Globals.fail
        return nothing
      elseif p != 0
        isprime(p) || error("p must be 0 or a prime integer")
        tbl = GAP.Globals.mod(tbl, GAP.julia_to_gap(p))
        if tbl == GAP.Globals.fail
          return nothing
        end
      end
      return GAPGroupCharacterTable(tbl, p)
    else
      error("no character table library available")
    end
end

# utility:
# turn integer values to strings, but replace `0` by `"."`,
# print irrationalities via the `Hecke.math_html` method for `nf_elem`
function matrix_of_strings(tbl::GAPGroupCharacterTable)
  n = length(tbl)
  m = Array{String}(undef, n, n)
  buf = Base.IOBuffer()
  for i in 1:n
    for j in 1:n
      val = tbl[i,j].data
      if iszero(val)
        m[i, j] = "."
      else
        # irrationalities involving `\zeta`
        m[i,j] = sprint(Hecke.math_html, val)
      end
    end
  end
  return m
end

# Produce LaTeX output if `"text/html"` is prescribed,
# via the `:TeX` attribute of the io context.
function Base.show(io::IO, ::MIME"text/html", tbl::GAPGroupCharacterTable)
  print(io, "\$")
  show(IOContext(io, :TeX => true), tbl)
  print(io, "\$")
end

# Produce a screen format without LaTeX markup but with unicode characters
# and sub-/superscripts if LaTeX output is not requested..
function Base.show(io::IO, tbl::GAPGroupCharacterTable)
    n = nrows(tbl)
    size = GAP.gap_to_julia(GAP.Globals.Size(tbl.GAPTable))
    primes = [x[1] for x in collect(factor(size))]
    sort!(primes)

    # Compute the factored centralizer orders.
    cents = GAP.gap_to_julia(GAP.Globals.SizesCentralizers(tbl.GAPTable))
    fcents = [collect(factor(x)) for x in cents]
    d = Dict([p => fill(".", n) for p in primes]...)
    for i in 1:n
      for pair in fcents[i]
        d[pair[1]][i] = string(pair[2])
      end
    end
    cents_strings = [d[p] for p in primes]

    # Compute display format for power maps.
    names = Vector{String}(GAP.Globals.ClassNames(tbl.GAPTable))
    pmaps = Vector{Any}(GAP.Globals.ComputedPowerMaps(tbl.GAPTable))
    power_maps_primes = []
    power_maps_strings = []
    for i in 1:length(pmaps)
      map = pmaps[i]
      if map != nothing
        push!(power_maps_primes, string(i)*"P")
        push!(power_maps_strings, names[map])
      end
    end

    empty = ["" for i in 1:n]

    # Create the IO context.
    ioc = IOContext(io,
      # header (an array of strings):
      # name of the table and separating empty line
      :header => [String(GAP.Globals.Identifier(tbl.GAPTable)),
                  ""],
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
      :labels_row => ["\\chi_{" * string(i) * "}" for i in 1:n],

      # corner (a column vector is sufficient):
      # primes in the centralizer rows,
      # separating empty line,
      # separating empty line,
      # primes in the power map rows,
      # separating empty line,
      :corner => vcat( string.(primes), [ "", "" ], power_maps_primes, [ "" ] ),
    )

    # print the table
    labelled_matrix_formatted(ioc, matrix_of_strings(tbl))
end

# print: abbreviated form
function Base.print(io::IO, tbl::GAPGroupCharacterTable)
    gaptbl = tbl.GAPTable
    if GAP.Globals.HasUnderlyingGroup(gaptbl)
      id = string(GAP.Globals.UnderlyingGroup(gaptbl))
    else
      id = "\"" * String(GAP.Globals.Identifier(gaptbl)) * "\""
    end
    print(io, "character_table($id)")
end

Base.length(tbl::GAPGroupCharacterTable) = return length(GAP.Globals.Irr(tbl.GAPTable))
AbstractAlgebra.nrows(tbl::GAPGroupCharacterTable) = return length(GAP.Globals.Irr(tbl.GAPTable))
AbstractAlgebra.ncols(tbl::GAPGroupCharacterTable) = return length(GAP.Globals.Irr(tbl.GAPTable))

function Base.getindex(tbl::GAPGroupCharacterTable, i::Int)
    return group_class_function(tbl, GAP.Globals.Irr(tbl.GAPTable)[i])
end

function Base.getindex(tbl::GAPGroupCharacterTable, i::Int, j::Int)
    val = GAP.Globals.Irr(tbl.GAPTable)[i, j]
    return QabModule.QabElem(val)
end

function Base.mod(tbl::GAPGroupCharacterTable, p::Int)
    isprime(p) || error("p must be a prime integer")
    tbl.characteristic == 0 || error("tbl mod p only for ordinary table tbl")

    modtbls = AbstractAlgebra.get_special(tbl, :brauer_tables)
    if modtbls == nothing
      modtbls = Dict{Int,Any}()
      AbstractAlgebra.set_special(tbl, :brauer_tables => modtbls)
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

function decomposition_matrix(tbl::GAPGroupCharacterTable)
    isprime(tbl.characteristic) || error("characteristic of tbl must be a prime integer")
    return Matrix{Int}(GAP.Globals.DecompositionMatrix(tbl.GAPTable))
end


#############################################################################
##
##  characters
##
abstract type GroupClassFunction end

struct GAPGroupClassFunction <: GroupClassFunction
    table::GAPGroupCharacterTable
    values::GAP.GapObj
end

function Base.show(io::IO, chi::GAPGroupClassFunction)
    print(io, "group_class_function(" * string(chi.table) * ", " * string(values(chi)) * ")")
end

function values(chi::GAPGroupClassFunction)
    gapvalues = GAP.gap_to_julia(GAP.Globals.ValuesOfClassFunction(chi.values), recursive = false)
    return [QabModule.QabElem(x) for x in gapvalues]
end

function group_class_function(tbl::GAPGroupCharacterTable, values::GAP.GapObj)
    GAP.Globals.IsClassFunction(values) || error("values must be a class function")
    return GAPGroupClassFunction(tbl, values)
end

function group_class_function(tbl::GAPGroupCharacterTable, values::Vector{QabModule.QabElem})
    return GAPGroupClassFunction(tbl, GAP.Globals.ClassFunction(tbl.GAPTable, GAP.julia_to_gap([gap_cyclotomic(x) for x in values])))
end

#function group_class_function(G::Oscar.GAPGroup, values::Vector{QabModule.QabElem})
#end

function trivial_character(tbl::GAPGroupCharacterTable)
    val = QabModule.QabElem(1)
    return group_class_function(tbl, [val for i in 1:length(tbl)])
end

#function trivial_character(G::Oscar.GAPGroup)
#end

Base.length(chi::GAPGroupClassFunction) = return length(chi.values)

function Nemo.degree(chi::GAPGroupClassFunction)
    val = values(chi)[1]
    return Nemo.coeff(val.data, 0)
#TODO: make sure that the returned value is really equal to the 1st entry
#TODO: change result from fmpq to fmpz?
end

Base.getindex(chi::GAPGroupClassFunction, i::Int) = QabModule.QabElem(GAP.Globals.ValuesOfClassFunction(chi.values)[i])

function Base.:(==)(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    chi.table === psi.table || error("character tables must be identical")
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
    val = QabModule.QabElem(0)
    return group_class_function(chi.table, [val for i in 1:length(chi)])
end

Base.one(chi::GAPGroupClassFunction) = trivial_character(chi.table)

function scalar_product(chi::GAPGroupClassFunction, psi::GAPGroupClassFunction)
    chi.table === psi.table || error("the two class functions belong to different character tables")
    return Nemo.fmpz(GAP.gap_to_julia(GAP.Globals.ScalarProduct(chi.values, psi.values)))
end

#TODO: chi(g), ...

function Base.:*(n::T, chi::GAPGroupClassFunction) where T <: Integer
    return GAPGroupClassFunction(chi.table, n * chi.values)
end

function Base.:^(chi::GAPGroupClassFunction, n::Int)
    return GAPGroupClassFunction(chi.table, chi.values ^ n)
end

end  # module
