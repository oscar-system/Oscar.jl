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


Some open questions:

- Would it be advisable to introduce a `ClassFunctionSpace` object?
  Magma does this, and one benefit would be that this would be a natural
  object for constructing class functions, via `S(values)`, say.
  On the other hand, one can also use the character table object itself
  for this purpose.

- What is the recommended output for a function `character_field`?
  In general, we get an abelian number field, and one possibility would be
  to return an embedding into the enveloping cyclotomic field;
  or would an embedding into `QabField` be more natural?
  An embedding into the minimal splitting field might be too large
  for being considered.
  (And here the question about caching arises, one wants to reuse
  field objects which have already been created.)

- What is the recommended way to convert those `QabElem` or `nf_elem` values
  that lie in the base field to integers (or rationals)?
  For example, it makes sense that `degree(chi)` returns an integer
  (and not a `Qab_elem`) if `chi` is a (virtual) character.

- What is the recommended way to implement caching in Julia?
  The most critical point is that a character table of a group
  relies on a fixed ordering of conjugacy classes.
  (O.k., this is important enough such that we should provide a field
  in the table struct hat holds the array of conjugacy classes.)

  The next instance is the connection from the group to its character table.
  If one wants something like `trivial_character(G)` then the group must
  know how to fetch the same character table object in each call;
  perhaps we do not want to support such functions at all.

  Other information is optional, it concerns for example the question
  whether a given class function is a (virtual) character;
  sometimes this is known on construction (the sum of two characters is
  a character, inducing a character from a subgroup yields a character),
  but often this is not known a priori (difference of two characters).


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

import Base: getindex, length, mod, one, print, show, zero

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
#abstract type GroupCharacterTableOrdinary <: GroupCharacterTable end
#abstract type GroupCharacterTableBrauer <: GroupCharacterTable end

struct GAPGroupCharacterTable <: GroupCharacterTable
    GAPTable::GAP.GapObj  # the character table object
    characteristic::Int
end

struct GAPGroupCharacterTableWithGroup <: GroupCharacterTable
    GAPGroup::Oscar.GAPGroup    # the underlying group
    GAPTable::GAP.GapObj  # the character table object
    characteristic::Int
end

const GAPGroupCharacterTableWithOrWithoutGroup= Union{GAPGroupCharacterTable,GAPGroupCharacterTableWithGroup}

function character_table(G::Oscar.GAPGroup, p::Int = 0)
    tbl = GAP.Globals.CharacterTable(G.X)
    if p != 0
      isprime(p) || error("p must be 0 or a prime integer")
      tbl = mod(tbl, p)
    end
    return GAPGroupCharacterTableWithGroup(G, tbl, p)
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

function Base.show(io::IO, tbl::GAPGroupCharacterTableWithOrWithoutGroup)
    print(io, GAP.CSTR_STRING(GAP.Globals.StringDisplayObj(tbl.GAPTable)))
end

function Base.print(io::IO, tbl::GAPGroupCharacterTableWithOrWithoutGroup)
    gaptbl = tbl.GAPTable
    if GAP.Globals.HasUnderlyingGroup(gaptbl)
      id = string(GAP.Globals.UnderlyingGroup(gaptbl))
    else
      id = "\"" * String(GAP.Globals.Identifier(gaptbl)) * "\""
    end
    print(io, "character_table($id)")
end

Base.length(tbl::GAPGroupCharacterTableWithOrWithoutGroup) = return length(GAP.Globals.Irr(tbl.GAPTable))

function Base.getindex(tbl::GAPGroupCharacterTableWithOrWithoutGroup, i::Int)
    return group_class_function(tbl, GAP.Globals.Irr(tbl.GAPTable)[i])
end

function Base.getindex(tbl::GAPGroupCharacterTableWithOrWithoutGroup, i::Int, j::Int)
    val = GAP.Globals.Irr(tbl.GAPTable)[i, j]
    return QabModule.QabElem(val)
end

function Base.mod(tbl::GAPGroupCharacterTable, p::Int)
    isprime(p) || error("p must be a prime integer")
    tbl.characteristic == 0 || error("tbl mod p only for ordinary table tbl")
    return GAPGroupCharacterTable(mod(tbl.GAPTable, p), p)
end

function Base.mod(tbl::GAPGroupCharacterTableWithGroup, p::Int)
    isprime(p) || error("p must be a prime integer")
    tbl.characteristic == 0 || error("tbl mod p only for ordinary table tbl")
    return GAPGroupCharacterTableWithGroup(tbl.GAPGroup, mod(tbl.GAPTable, p), p)
end

function decomposition_matrix(tbl::GAPGroupCharacterTableWithOrWithoutGroup)
    isprime(tbl.characteristic) || error("characteristic of tbl must be a prime integer")
    return GAP.gap_to_julia(Matrix{Int}, GAP.Globals.DecompositionMatrix(tbl.GAPTable))
end


#############################################################################
##
##  characters
##
abstract type GroupClassFunction end

struct GAPGroupClassFunction <: GroupClassFunction
    table::GAPGroupCharacterTableWithOrWithoutGroup
    values::GAP.GapObj
end

function Base.show(io::IO, chi::GAPGroupClassFunction)
    print(io, "group_class_function(" * string(chi.table) * ", " * string(values(chi)) * ")")
end

function values(chi::GAPGroupClassFunction)
    gapvalues = GAP.gap_to_julia(GAP.Globals.ValuesOfClassFunction(chi.values), recursive = false)
    return [QabModule.QabElem(x) for x in gapvalues]
end

function group_class_function(tbl::GAPGroupCharacterTableWithOrWithoutGroup, values::GAP.GapObj)
    GAP.Globals.IsClassFunction(values) || error("values must be a class function")
    return GAPGroupClassFunction(tbl, values)
end

function group_class_function(tbl::GAPGroupCharacterTableWithOrWithoutGroup, values::Vector{QabModule.QabElem})
    return GAPGroupClassFunction(tbl, GAP.Globals.ClassFunction(tbl.GAPTable, GAP.julia_to_gap([gap_cyclotomic(x) for x in values])))
end

#function group_class_function(G::Oscar.GAPGroup, values::Vector{QabModule.QabElem})
#end

function trivial_character(tbl::GAPGroupCharacterTableWithOrWithoutGroup)
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

#TODO: chi(g), 2*chi, chi^3, ...

end
