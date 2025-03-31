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
Permutation group of degree 5 and order 60

julia> atlas_group(MatrixGroup, "A5")
Matrix group of degree 4
  over prime field of characteristic 2

julia> atlas_group("M11")  # Mathieu group M11
Permutation group of degree 11 and order 7920

julia> atlas_group("M")  # Monster group M
ERROR: ArgumentError: the group atlas does not provide a representation for M
[...]
```
"""
function atlas_group(name::String)
  G = GAP.Globals.AtlasGroup(GapObj(name))
  if G === GAP.Globals.fail
    # Check whether the table of contents provides a representation.
    l = GAP.Globals.AGR.MergedTableOfContents(GapObj("all"), GapObj(gapname))::GapObj
    if length(l) != 0
      error("cannot access the group atlas, perhaps the download failed")
    else
      throw(ArgumentError("the group atlas does not provide a representation for $name"))
    end
  end
  return _oscar_group(G)
end

function atlas_group(::Type{T}, name::String) where T <: Union{PermGroup, MatrixGroup}
  if T === PermGroup
    G = GAP.Globals.AtlasGroup(GapObj(name), GAP.Globals.IsPermGroup, true)::GapObj
  else
    G = GAP.Globals.AtlasGroup(GapObj(name), GAP.Globals.IsMatrixGroup, true)::GapObj
  end

  if G === GAP.Globals.fail
    # Check whether the table of contents provides a representation.
    l = GAP.Globals.AGR.MergedTableOfContents(GapObj("all"), GapObj(gapname))::GapObj
    permtype = GapObj("perm")::GapObj
    if T === PermGroup
      l = filter(x -> x.type == permtype, l)
    else
      l = filter(x -> x.type != permtype, l)
    end
    if length(l) != 0
      error("cannot access the group atlas, perhaps the download failed")
    else
      throw(ArgumentError("the group atlas does not provide a representation of type $T for $name"))
    end
  end
  return _oscar_group(G)
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
 Dict(:constituents => [1, 4], :repname => "A5G1-p5B0", :degree => 5, :name => "A5")

julia> atlas_group(info[1])
Permutation group of degree 5 and order 60

```
"""
function atlas_group(info::Dict)
  gapname = info[:name]
  l = GAP.Globals.AGR.MergedTableOfContents(GapObj("all"), GapObj(gapname))::GapObj
  repname = info[:repname]
  pos = findfirst(r -> String(r.repname) == repname, Vector{GapObj}(l))
  @req (pos !== nothing) "no Atlas group for $repname"
  G = GAP.Globals.AtlasGroup(l[pos])
  # The table of contents knows about the requested representation.
  # If the result is `fail` then this is likely due to a download problem.
  @req (G !== GAP.Globals.fail) "cannot access the representation $repname from the group atlas, perhaps the download failed"

  if haskey(info, :base_ring_iso)
    # make sure that the given ring is used
    deg = GAP.Globals.DimensionOfMatrixGroup(G)
    iso = info[:base_ring_iso]
    ring = domain(iso)
    matgrp = matrix_group(ring, deg)
    matgrp.ring_iso = iso
    matgrp.X = G
    return matgrp
  else
    return _oscar_group(G)
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
Permutation group of degree 11 and order 720

julia> order(h1)  # largest maximal subgroup of M11
720

julia> h2, emb = atlas_subgroup("M11", 1);  h2
Permutation group of degree 11 and order 720

julia> h3, emb = atlas_subgroup(MatrixGroup, "M11", 1 );  h3
Matrix group of degree 10
  over prime field of characteristic 2

julia> info = all_atlas_group_infos("M11", degree => 11);

julia> h4, emb = atlas_subgroup(info[1], 1);  h4
Permutation group of degree 11 and order 720
```
"""
function atlas_subgroup(G::GAPGroup, nr::Int)
  @req GAP.Globals.HasAtlasRepInfoRecord(GapObj(G)) "$G was not constructed with atlas_group"
  info = GAP.Globals.AtlasRepInfoRecord(GapObj(G))
  @req (info.groupname == info.identifier[1]) "$G was not constructed with atlas_group"
  H = GAP.Globals.AtlasSubgroup(GapObj(G), nr)
  if H === GAP.Globals.fail
    # Check whether the table of contents provides the SLP in question.
    slp = GAP.Globals.AtlasProgramInfo(info.groupname, GapObj("maxes"), nr)
    name = string(info.groupname)
    if slp === GAP.Globals.fail
      error("the group atlas does not provide the restriction to the $nr-th class of maximal subgroups of $name")
    else
      # This is likely due to a download problem.
      error("cannot access the slp for the restriction to the $nr-th class of maximal subgroups of $name")
    end
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
 Dict(:constituents => [1, 4], :repname => "A5G1-p5B0", :degree => 5, :name => "A5")
 Dict(:constituents => [1, 5], :repname => "A5G1-p6B0", :degree => 6, :name => "A5")

julia> atlas_group(info[1])
Permutation group of degree 5 and order 60

julia> info = all_atlas_group_infos("A5", dim => 4, characteristic => 3)
1-element Vector{Dict{Symbol, Any}}:
 Dict(:dim => 4, :constituents => [4], :repname => "A5G1-f3r4B0", :name => "A5")

julia> atlas_group(info[1])
Matrix group of degree 4
  over prime field of characteristic 3

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
      @req haskey(_atlas_group_filter_attrs, func) "Function not supported"
      expected_type, gapfunc, _ = _atlas_group_filter_attrs[func]
      @req data isa expected_type "bad argument $(data) for function $(func)"
      if func === base_ring
        # we will need the isomorphism later on
        iso = iso_oscar_gap(data)
        push!(gapargs, gapfunc, codomain(iso))
      elseif func === character
        push!(gapargs, gapfunc, GapObj(data))
      else
        # we can translate `data` to GAP
        push!(gapargs, gapfunc, GAP.Obj(data))
      end
    elseif arg isa Function
      # handle e.g. `is_primitive` or `! is_primitive`
      func = arg
      @req haskey(_atlas_group_filter_attrs, func) "Function not supported"
      expected_type, gapfunc, default = _atlas_group_filter_attrs[func]
      @req default !== nothing "missing argument for function $(func)"
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

    # the character may be stored
    if hasproperty(r, :constituents)
      d[:constituents] = Vector{Int}(r.constituents)
    end

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
    number_of_atlas_groups([::Type{T}, ]name::String) where T <: Union{PermGroup, MatrixGroup}

Return the number of groups from the Atlas of Group Representations
whose isomorphism type is given by `name` and have the type `T`.

# Examples
```jldoctest
julia> number_of_atlas_groups("A5")
18

julia> number_of_atlas_groups(PermGroup, "A5")
3

julia> number_of_atlas_groups(MatrixGroup, "A5")
15

```
"""
function number_of_atlas_groups(name::String)
  return length(GAP.Globals.AllAtlasGeneratingSetInfos(GapObj(name))::GapObj)
end

function number_of_atlas_groups(::Type{T}, name::String) where T <: Union{PermGroup, MatrixGroup}
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
    return straight_line_program(slp.program)
  else
    error("not yet ...")
  end
end


###################################################################
# Show overviews of available information.
###################################################################

"""
    show_atlas_info()
    show_atlas_info(groupnames::Vector{String})
    show_atlas_info(groupname::String)

Print information concerning the group atlas about the groups
given by the arguments.

If the only argument is a vector `groupnames` of strings then
print one line about each Atlas group whose name is in `groupnames`,
showing the available information.

If there are no arguments then the same happens as if the vector of all
names of Atlas groups were entered as `groupnames`.

If the only argument is a string `groupname`, the name of an Atlas group `G`,
then print an overview of available representations and straight line
programs for `G`.

The overview for several groups is shown as a table with the following
columns.

- `group`:
  the name of `G`

- `#`:
  the number of faithful representations stored for `G`

- `maxes`:
  the number of available straight line programs for computing
  generators of maximal subgroups of `G`,

- `cl`:
  a `+` sign if at least one program for computing representatives of
  conjugacy classes of elements of `G` is stored,

- `cyc`:
  a `+` sign if at least one program for computing representatives of
  classes of maximally cyclic subgroups of `G` is stored,

- `out`:
  descriptions of outer automorphisms of `G` for which at least one
  program is stored,

- `fnd`:
  a `+` sign if at least one program is available for finding standard
  generators,

- `chk`:
  a `+` sign if at least one program is available for checking whether a
  set of generators is a set of standard generators, and

- `prs`:
  a `+` sign if at least one program is available that encodes a
  presentation.

The overview for a single group is shown in the form of two lists,
one for available representations, one for available straight line programs.

Each available representation is described by either `G <= Sym(nid)` or
`G <= GL(nid,descr)`, where the former means a permutation representation
of degree `n` and the latter means a matrix ring of dimension `n` over a
base ring described by `descr`; if `descr` is an integer then this means
a finite field with this cardinality.
In both cases, `id` is a (perhaps empty) identifier that distinguishes
different representations with the same `n`.
If known then information about transitivity, rank, point stabillizers of
permutation representations and about the character of matrix representations
is shown as well.

Below the representations, the programs available for `groupname` are listed.
"""
show_atlas_info() = show_atlas_info("all")

function show_atlas_info(groupnames::Vector{String})
  # Set a GAP user preference according to Oscar's unicode setting.
  AtlasRep = GapObj("AtlasRep")
  DisplayFunction = GapObj("DisplayFunction")
  if !is_unicode_allowed()
    oldDisplayFunction = GAP.Globals.UserPreference(AtlasRep, DisplayFunction)
    GAP.Globals.SetUserPreference(AtlasRep, DisplayFunction, GapObj("Print"));
  end

  # Remove the marker of non-core data.
  AtlasRepMarkNonCoreData = GapObj("AtlasRepMarkNonCoreData")
  oldmarker = GAP.Globals.UserPreference(AtlasRep, AtlasRepMarkNonCoreData)
  GAP.Globals.SetUserPreference(AtlasRep, AtlasRepMarkNonCoreData, GapObj(""))

  # Show the info
  info = GAP.Globals.AGR.StringAtlasInfoOverview(GapObj(groupnames, true),
             GapObj([])) # unconditional
  for line in info
    println(String(line));
  end

  # Reset the changed values.
  if !is_unicode_allowed()
    GAP.Globals.SetUserPreference(AtlasRep, DisplayFunction, oldDisplayFunction)
  end
  GAP.Globals.SetUserPreference(AtlasRep, AtlasRepMarkNonCoreData, oldmarker)

  return
end

function show_atlas_info(groupname::String)
  # Set a GAP user preference according to Oscar's unicode setting.
  AtlasRep = GapObj("AtlasRep")
  DisplayFunction = GapObj("DisplayFunction")
  if !is_unicode_allowed()
    oldvalue = GAP.Globals.UserPreference(AtlasRep, DisplayFunction)
    GAP.Globals.SetUserPreference(AtlasRep, DisplayFunction, GapObj("Print"));
  end

  # Remove the marker of non-core data.
  AtlasRepMarkNonCoreData = GapObj("AtlasRepMarkNonCoreData")
  oldmarker = GAP.Globals.UserPreference(AtlasRep, AtlasRepMarkNonCoreData)
  GAP.Globals.SetUserPreference(AtlasRep, AtlasRepMarkNonCoreData, GapObj(""))

  # Show the info
  if groupname == "all"
    # one-line overview of all supported groups
    info = GAP.Globals.AGR.StringAtlasInfoOverview(GapObj(groupnames, true),
               GapObj([])) # unconditional
  else
    # detailed overview for one group
    info = GAP.Globals.AGR.StringAtlasInfoGroup(GapObj([groupname], true))
  end
  for line in info
    println(String(line));
  end

  # Reset the changed values.
  if !is_unicode_allowed()
    GAP.Globals.SetUserPreference(AtlasRep, DisplayFunction, oldvalue)
  end
  GAP.Globals.SetUserPreference(AtlasRep, AtlasRepMarkNonCoreData, oldmarker)

  return
end
