export
    atlas_group,
    number_atlas_groups

###################################################################
# Groups from the Atlas of Group Representations
###################################################################

"""
    atlas_group([::Type{T}, ]name::String) where T <: Union{PermGroup, MatrixGroup}

Return a group from the Atlas of Group Representations
whose isomorphism type is given by `name` and have the type `T`.

# Examples
```jldoctest
julia> atlas_group("A5")  # alternating group A5
Group([ (1,2)(3,4), (1,3,5) ])

julia> atlas_group("M11")  # Mathieu group M11
Group([ (2,10)(4,11)(5,7)(8,9), (1,4,3,8)(2,5,6,9) ])

```
"""
function atlas_group(name::String)
  G = GAP.Globals.AtlasGroup(GAP.GapObj(name))
  G === GAP.Globals.fail && error("the group atlas does not provide a representation for $name")
  T = _get_type(G)
  return T(G)
end

function atlas_group(::Type{T}, name::String) where T <: Union{PermGroup, MatrixGroup}
  if T === PermGroup
    G = GAP.Globals.AtlasGroup(GAP.GapObj(name), GAP.Globals.IsPermGroup, true)
  else
    G = GAP.Globals.AtlasGroup(GAP.GapObj(name), GAP.Globals.IsMatrixGroup, true)
  end
  G === GAP.Globals.fail && error("the group atlas does not provide a representation of type $T for $name")
  TT = _get_type(G)
  return TT(G)
end

"""
    number_atlas_groups([::Type{T}, ]name::String) where T <: Union{PermGroup, MatrixGroup}

Return the number of groups from the Atlas of Group Representations
whose isomorphism type is given by `name` and have the type `T`.

# Examples
```
julia> number_atlas_groups("A5")
18

julia> number_atlas_groups(PermGroup, "A5")
3

```
"""
function number_atlas_groups(name::String)
  return length(GAP.Globals.AllAtlasGeneratingSetInfos(GAP.GapObj(name)))
end

function number_atlas_groups(::Type{T}, name::String) where T <: Union{PermGroup, MatrixGroup}
  if T === PermGroup
    return length(GAP.Globals.AllAtlasGeneratingSetInfos(
                    GAP.GapObj(name), GAP.Globals.IsPermGroup, true))
  else
    return length(GAP.Globals.AllAtlasGeneratingSetInfos(
                    GAP.GapObj(name), GAP.Globals.IsMatrixGroup, true))
  end
end
