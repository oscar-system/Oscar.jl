export
    all_atlas_group_infos,
    atlas_group,
    atlas_program,
    atlas_subgroup,
    number_atlas_groups

###################################################################
# Groups from the Atlas of Group Representations
###################################################################

"""
    atlas_group([::Type{T}, ]name::String) where T <: Union{PermGroup, MatrixGroup}

Return a group from the Atlas of Group Representations
whose isomorphism type is given by `name` and have the type `T`.
If `T` is not given then `PermGroup` is chosen if a permutation group
for `name` is available, and `MatrixGroup` otherwise.

# Examples
```jldoctest
julia> atlas_group("A5")  # alternating group A5
Group([ (1,2)(3,4), (1,3,5) ])

julia> atlas_group(MatrixGroup, "A5")
Matrix group of degree 4 over Galois field with characteristic 2

julia> atlas_group("M11")  # Mathieu group M11
Group([ (2,10)(4,11)(5,7)(8,9), (1,4,3,8)(2,5,6,9) ])

julia> atlas_group("M")  # Monster group M
ERROR: the group atlas does not provide a representation for M
```
"""
function atlas_group(name::String)
  G = GAP.Globals.AtlasGroup(GapObj(name))
  G === GAP.Globals.fail && error("the group atlas does not provide a representation for $name")
  T = _get_type(G)
  return T(G)
end

function atlas_group(::Type{T}, name::String) where T <: Union{PermGroup, MatrixGroup}
  if T === PermGroup
    G = GAP.Globals.AtlasGroup(GapObj(name), GAP.Globals.IsPermGroup, true)::GapObj
  else
    G = GAP.Globals.AtlasGroup(GapObj(name), GAP.Globals.IsMatrixGroup, true)::GapObj
  end
  G === GAP.Globals.fail && error("the group atlas does not provide a representation of type $T for $name")
  TT = _get_type(G)
  return TT(G)
end


"""
    atlas_group(info::Dict)

Return the group from the Atlas of Group Representations
that is defined by `info`.
Typically, `info` is obtained from [`all_atlas_group_infos`](@ref).

# Examples
```jldoctest
julia> info = all_atlas_group_infos("A5", degree => 5)
1-element Vector{Dict{Symbol, Any}}:
 Dict(:repname => "A5G1-p5B0", :degree => 5, :name => "A5")

julia> atlas_group(info[1])
Group([ (1,2)(3,4), (1,3,5) ])

```
"""
function atlas_group(info::Dict)
  gapname = info[:name]
  l = GAP.Globals.AGR.MergedTableOfContents(GapObj("all"), GapObj(gapname))::GapObj
  pos = findfirst(r -> String(r.repname) == info[:repname], Vector{GAP.GapObj}(l))
  pos === nothing && error("no Atlas group for $info")
  G = GAP.Globals.AtlasGroup(l[pos])
  G === GAP.Globals.fail && error("the group atlas does not provide a representation for $info")

  if haskey(info, :base_ring_iso)
    # make sure that the given ring is used
    deg = GAP.Globals.DimensionOfMatrixGroup(G)
    iso = info[:base_ring_iso]
    ring = domain(iso)
    matgrp = MatrixGroup(deg, ring)
    matgrp.ring_iso = iso
    matgrp.X = G
    return matgrp
  else
    TT = _get_type(G)
    return TT(G)
  end
end


"""
    atlas_subgroup(G::GAPGroup, nr::Int)
    atlas_subgroup([::Type{T}, ]name::String, nr::Int) where T <: Union{PermGroup, MatrixGroup}
    atlas_subgroup(info::Dict, nr::Int)

Return a pair `(H, emb)` where `H` is a representative of the `nr`-th class
of maximal subgroups of the group `G`,
and `emb` is an embedding of `H` into `G`.

The group `G` can be given as the first argument,
in this case it is assumed that `G` has been created with
[`atlas_group`](@ref).
Otherwise `G` is the group obtained by calling [`atlas_group`](@ref)
with (`T` and) `name` or with `info`.

If the Atlas of Group Representations does not provide the information to
compute `G` or to compute generators of `H` from `G` then an exception is
thrown.

# Examples
```jldoctest
julia> g = atlas_group("M11");  # Mathieu group M11

julia> h1, emb = atlas_subgroup(g, 1);  h1
Group([ (1,4)(2,10)(3,7)(6,9), (1,6,10,7,11,3,9,2)(4,5) ])

julia> order(h1)  # largest maximal subgroup of M11
720

julia> h2, emb = atlas_subgroup("M11", 1);  h2
Group([ (1,4)(2,10)(3,7)(6,9), (1,6,10,7,11,3,9,2)(4,5) ])

julia> h3, emb = atlas_subgroup(MatrixGroup, "M11", 1 );  h3
Matrix group of degree 10 over Galois field with characteristic 2

julia> info = all_atlas_group_infos("M11", degree => 11);

julia> h4, emb = atlas_subgroup(info[1], 1);  h4
Group([ (1,4)(2,10)(3,7)(6,9), (1,6,10,7,11,3,9,2)(4,5) ])
```
"""
function atlas_subgroup(G::GAPGroup, nr::Int)
  GAP.Globals.HasAtlasRepInfoRecord(G.X) || error("$G was not constructed with atlas_group")
  info = GAP.Globals.AtlasRepInfoRecord(G.X)
  info.groupname == info.identifier[1] || error("$G was not constructed with atlas_group")
  H = GAP.Globals.AtlasSubgroup(G.X, nr)
  if H === GAP.Globals.fail
    name = string(info.groupname)
    error("the group atlas does not provide the restriction to the $nr-th class of maximal subgroups of $name")
  end
  return _as_subgroup(G, H)
end

atlas_subgroup(name::String, nr::Int) = atlas_subgroup(atlas_group(name), nr)

function atlas_subgroup(::Type{T}, name::String, nr::Int) where T <: Union{PermGroup, MatrixGroup}
  return atlas_subgroup(atlas_group(T, name), nr)
end

atlas_subgroup(info::Dict, nr::Int) = atlas_subgroup(atlas_group(info), nr)


"""
    all_atlas_group_infos(name::String, L...)

Return the vector of dictionaries that describe Atlas groups
whose isomorphism types are given by `name` and
which satisfy the conditions in `L`.
These conditions may be of one of the following forms:

- `func => intval` selects groups for which the function `func` returns `intval`
- `func => list` selects groups for which the function `func` returns any element inside `list`
- `func` selects groups for which the function `func` returns `true`
- `!func` selects groups for which the function `func` returns `false`

The following functions are currently supported as values for `func`:

For permutation groups

- `degree`
- `is_primitive`
- `is_transitive`
- `rank_action`
- `transitivity`

and for matrix groups

- `base_ring`
- `character`
- `characteristic`
- `dim`

# Examples
```jldoctest
julia> info = all_atlas_group_infos("A5", degree => [5, 6])
2-element Vector{Dict{Symbol, Any}}:
 Dict(:repname => "A5G1-p5B0", :degree => 5, :name => "A5")
 Dict(:repname => "A5G1-p6B0", :degree => 6, :name => "A5")

julia> atlas_group(info[1])
Group([ (1,2)(3,4), (1,3,5) ])

julia> info = all_atlas_group_infos("A5", dim => 4, characteristic => 3)
1-element Vector{Dict{Symbol, Any}}:
 Dict(:dim => 4, :repname => "A5G1-f3r4B0", :name => "A5")

julia> atlas_group(info[1])
Matrix group of degree 4 over Galois field with characteristic 3

```
"""
function all_atlas_group_infos(name::String, L...)
  iso = nothing

  # scan the given conditions
  gapargs = Any[GapObj(name)]
  for arg in L
    if arg isa Pair
      # handle e.g. `is_primitive => false`
      func = arg[1]
      data = arg[2]
      haskey(_atlas_group_filter_attrs, func) || throw(ArgumentError("Function not supported"))
      expected_type, gapfunc, _ = _atlas_group_filter_attrs[func]
      data isa expected_type || throw(ArgumentError("bad argument $(data) for function $(func)"))
      if func === base_ring
        # we will need the isomorphism later on
        iso = iso_oscar_gap(data)
        push!(gapargs, gapfunc, codomain(iso))
      elseif func === character
        push!(gapargs, gapfunc, data.values)
      else
        # we can translate `data` to GAP
        push!(gapargs, gapfunc, GAP.Obj(data))
      end
    elseif arg isa Function
      # handle e.g. `is_primitive` or `! is_primitive`
      func = arg
      haskey(_atlas_group_filter_attrs, func) || throw(ArgumentError("Function not supported"))
      expected_type, gapfunc, default = _atlas_group_filter_attrs[func]
      default !== nothing || throw(ArgumentError("missing argument for function $(func)"))
      push!(gapargs, gapfunc, default)
    else
      throw(ArgumentError("expected a function or a pair, got $arg"))
    end
  end

  # evaluate the conditions in GAP
  res_GAP = GAP.Globals.AllAtlasGeneratingSetInfos(gapargs...)::GapObj

  # translate the records to dictionaries
  res = Dict{Symbol, Any}[]
  for r in res_GAP
    # groupname and repname are always present
    d = Dict{Symbol, Any}(:name => string(r.groupname), :repname => string(r.repname))

    # permutation groups have a degree
    if hasproperty(r, :p)
      d[:degree] = r.p
    end

    # matrix groups have dim
    if hasproperty(r, :dim)
      d[:dim] = r.dim
    end

    # store a given ring
    if iso !== nothing
      d[:base_ring_iso] = iso
    end

    push!(res, d)
  end

  return res
end


"""
    number_atlas_groups([::Type{T}, ]name::String) where T <: Union{PermGroup, MatrixGroup}

Return the number of groups from the Atlas of Group Representations
whose isomorphism type is given by `name` and have the type `T`.

# Examples
```jldoctest
julia> number_atlas_groups("A5")
18

julia> number_atlas_groups(PermGroup, "A5")
3

julia> number_atlas_groups(MatrixGroup, "A5")
15

```
"""
function number_atlas_groups(name::String)
  return length(GAP.Globals.AllAtlasGeneratingSetInfos(GapObj(name))::GapObj)
end

function number_atlas_groups(::Type{T}, name::String) where T <: Union{PermGroup, MatrixGroup}
  if T === PermGroup
    return length(GAP.Globals.AllAtlasGeneratingSetInfos(
                    GapObj(name), GAP.Globals.IsPermGroup, true)::GapObj)
  else
    return length(GAP.Globals.AllAtlasGeneratingSetInfos(
                    GapObj(name), GAP.Globals.IsMatrixGroup, true)::GapObj)
  end
end

function atlas_program(name, paras...)
  if length(paras) == 1 && paras[1] == :classes
    slp = GAP.Globals.AtlasProgram(GapObj(name), GapObj("classes"))::GapObj
    slp === GAP.Globals.fail && return nothing
    gapcode = GAP.Globals.LinesOfStraightLineProgram(slp.program)::GapObj
    juliacode = GAP.gap_to_julia(gapcode, recursive = true)
    return Oscar.StraightLinePrograms.GAPSLProgram(juliacode)
  else
    error("not yet ...")
  end
end
