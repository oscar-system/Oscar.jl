struct SmallGroupsLibrary{T} <: GroupLibrary{T} end

@doc raw"""
    small_groups_library(::Type{T} = GAPGroup)

Return a handle for the library of groups of small order,
see [`Oscar.GroupLibrary`](@ref) for what one can do with it.

The primary key is the `order`; the group `L[n, i]` is the `i`-th group of
order `n`, up to isomorphism.
The groups are provided by the GAP packages `SmallGrp` [SmallGrp](@cite),
`SOTGrps` [SOTGrps](@cite) and `SglPPow` [SglPPow](@cite).

A group is returned as `PcGroup` if it is solvable and as `PermGroup`
otherwise; if `T` is given, each group is converted to type `T`.

[`find(L::GroupLibrary, filters...)`](@ref) supports the functions
`exponent`, `is_abelian`, `is_almost_simple`, `is_cyclic`, `is_nilpotent`,
`is_perfect`, `is_quasisimple`, `is_simple`, `is_sporadic_simple`,
`is_solvable`, `is_supersolvable`, `number_of_conjugacy_classes`, `order`.
For most orders the library knows `is_abelian`, `is_nilpotent`, `is_solvable`
and `is_supersolvable` of its groups without constructing them.

# Examples
```jldoctest
julia> L = small_groups_library()
Library of small groups

julia> L[60, 5]
Permutation group of degree 5 and order 60

julia> length(find(L, order => 512, !is_abelian))
10494183

julia> small_groups_library(PermGroup)[8, 3]
Permutation group of degree 4 and order 8
```
"""
small_groups_library(::Type{T}=GAPGroup) where {T} = SmallGroupsLibrary{T}()

_name(::SmallGroupsLibrary) = "small groups"
_primary_key(::SmallGroupsLibrary) = order
_key_type(::SmallGroupsLibrary) = ZZRingElem
_filter_attrs(::SmallGroupsLibrary) = _group_filter_attrs
_number(::SmallGroupsLibrary, n) = number_of_small_groups(n)

function Base.getindex(::SmallGroupsLibrary{GAPGroup}, n::IntegerUnion, i::IntegerUnion)
  return small_group(n, i)
end

function Base.getindex(::SmallGroupsLibrary{T}, n::IntegerUnion, i::IntegerUnion) where {T}
  return small_group(T, n, i)
end

identify(::SmallGroupsLibrary, G::GAPGroup) = small_group_identification(G)

has_groups(::SmallGroupsLibrary, n::IntegerUnion) = has_small_groups(n)
has_number_of_groups(::SmallGroupsLibrary, n::IntegerUnion) = has_number_of_small_groups(n)
has_identification(::SmallGroupsLibrary, n::IntegerUnion) = has_small_group_identification(n)

function _pairs(S::GroupLibrarySelection{<:SmallGroupsLibrary}, n)
  _check_groups(S.library, n)
  return ((id, S.library[id]) for id in _ids(S, n))
end

function _ids(S::GroupLibrarySelection{<:SmallGroupsLibrary}, n)
  isempty(S.rest) && return ((n, ZZRingElem(i)) for i in 1:Int(_number(S.library, n)))

  _check_groups(S.library, n)
  ids = GAP.Globals.IdsOfAllSmallGroups(GAP.Globals.Size, GAP.Obj(n), S.gaprest...)
  return ((n, ZZRingElem(id[2])) for id in ids)
end

function _count(S::GroupLibrarySelection{<:SmallGroupsLibrary}, n)
  isempty(S.rest) && return _number(S.library, n)

  _check_groups(S.library, n)
  # answered from tables where the library indexes the given properties
  return GAP.Globals.NumberSmallGroups(GAP.Globals.Size, GAP.Obj(n), S.gaprest...)
end

function _first(S::GroupLibrarySelection{<:SmallGroupsLibrary}, n)
  _check_groups(S.library, n)
  G = GAP.Globals.OneSmallGroup(GAP.Globals.Size, GAP.Obj(n), S.gaprest...)
  return G === GAP.Globals.fail ? nothing : _wrap(S.library, G)
end

################################################################################
#
#  Information
#
################################################################################

# the keys of the bibliography of SmallGrp in the bibliography of Oscar
const _SMALLGRP_REFERENCES = Dict(
  "Hoelder1893" => "Hld93",
  "Hoelder1895" => "Hld95",
  "Bagnera1898" => "Bag98",
  "OBr90" => "OBr90",
  "JNO90" => "JNO90",
  "OBr91" => "OBr91",
  "BescheEick98" => "BE99",
  "BescheEick1000" => "BE99*1",
  "BescheEick768" => "BE01",
  "BEO01" => "BEO02",
  "EOB98" => "EO99*1",
  "EOB99" => "EO99",
  "NOV04" => "NOV04",
  "Gir03" => "Gir03",
  "DEi05" => "DE05",
  "OV05" => "OV05",
)

# the selection criteria indexed by SmallGrp that `find` accepts
const _SMALLGRP_FILTERS = Dict(
  "IsAbelian" => is_abelian,
  "IsNilpotentGroup" => is_nilpotent,
  "IsSupersolvableGroup" => is_supersolvable,
  "IsSolvableGroup" => is_solvable,
)

# the other indexed criteria, in words
const _SMALLGRP_CRITERIA = Dict(
  "RankPGroup" => "rank",
  "PClassPGroup" => "p-class",
  "LGLength" => "LG-length",
  "FrattinifactorSize" => "order of the Frattini quotient",
  "FrattinifactorId" => "Frattini quotient",
)

# the references for the layers that other GAP packages add
function _small_groups_layer_references(layer::String)
  layer == "SmallGrp" && return ["SmallGrp"]
  layer == "SOTGrps" && return ["SOTGrps", "DEP22"]
  startswith(layer, "SglPPow") && return ["SglPPow", "OBr90", "OV05"]
  return String[]
end

# what SmallGrp knows about the groups of order `n` (from version 1.8 on)
function _small_groups_record(n::IntegerUnion)
  GAP.Globals.IsBoundGlobal(GapObj("SmallGroupsInformationRecord"))::Bool || return nothing
  r = GAP.Globals.SmallGroupsInformationRecord(GAP.Obj(n))
  return r === GAP.Globals.fail ? nothing : r
end

function _small_groups_references(r::GapObj)
  refs = _small_groups_layer_references(String(r.layer))
  hasproperty(r, :sources) || return refs
  return [refs; [get(_SMALLGRP_REFERENCES, String(k), String(k)) for k in r.sources]]
end

function information(::SmallGroupsLibrary)
  lines = [
    "SmallGrp $(_gap_package_version("smallgrp")): order at most 2000 except 1024,",
    "  cubefree order at most 50000, p^n for n at most 6, p^7 for p = 3, 5, 7, 11,",
    "  q^n p for q^n dividing 2^8, 3^6, 5^5 or 7^4, squarefree order,",
    "  order with at most 3 prime factors",
    "SOTGrps $(_gap_package_version("sotgrps")): order with at most 4 prime factors, p^4 q",
    "SglPPow $(_gap_package_version("sglppow")): p^7 for p > 11, 3^8",
  ]
  references = ["SmallGrp", "BEO02", "SOTGrps", "DEP22", "SglPPow"]
  return LibraryInformation("Library of small groups"; lines, references)
end

function information(L::SmallGroupsLibrary, n::IntegerUnion)
  title = "Small groups of order $n"
  has_number_of_groups(L, n) || return LibraryInformation(title; lines=["not available"])
  count = _count_line(number_of_small_groups(n), "group")
  if !has_groups(L, n)
    return LibraryInformation(title; lines=[count * ", which are not available"])
  end

  r = _small_groups_record(n)
  if r === nothing
    lines = [count, "identify: " * _availability_line(has_identification(L, n))]
    return LibraryInformation(title; lines, references=["SmallGrp"])
  end

  ranges = [(x.from, x.to, String(x.text)) for x in r.ranges]
  hasproperty(r, :sortedBy) && (count *= ", sorted by $(String(r.sortedBy))")
  isempty(ranges) || (count *= ":")
  lines = [count; _range_lines(ranges)]
  hasproperty(r, :comment) && push!(lines, String(r.comment))

  indexed = [String(x) for x in r.indexed]
  filters = Any[_SMALLGRP_FILTERS[x] for x in indexed if haskey(_SMALLGRP_FILTERS, x)]
  others = [_SMALLGRP_CRITERIA[x] for x in indexed if haskey(_SMALLGRP_CRITERIA, x)]
  isempty(indexed) ||
    push!(lines, "decided from precomputed data: " * join([string.(filters); others], ", "))

  push!(lines, "identify: " * _availability_line(r.idAvailable))
  layer = String(r.layer)
  push!(lines, layer == "SmallGrp" ? "layer $(r.layerNumber) of SmallGrp" : "layer $layer")
  return LibraryInformation(title; lines, references=_small_groups_references(r), filters)
end

function information(L::SmallGroupsLibrary, id::Tuple{IntegerUnion,IntegerUnion})
  n, i = id
  @req has_number_of_groups(L, n) "small groups of order $n are not available"
  N = number_of_small_groups(n)
  @req 1 <= i <= N "there are only $N groups of order $n, not $i"
  title = "Small group ($n, $i)"
  lines = ["group $i of the $N groups of order $n"]

  r = _small_groups_record(n)
  r === nothing && return LibraryInformation(title; lines, references=["SmallGrp"])

  for x in r.ranges
    x.from <= i <= x.to || continue
    range = x.from == x.to ? "" : "in the range $(x.from) - $(x.to): "
    push!(lines, range * String(x.text))
  end
  return LibraryInformation(title; lines, references=_small_groups_references(r))
end
