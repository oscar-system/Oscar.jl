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
    possible_class_fusions,
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
           
χ₁