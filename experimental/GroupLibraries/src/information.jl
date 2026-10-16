################################################################################
#
#  Information about the data of a library
#
################################################################################

# What `information` returns: a title, the lines of the overview, the keys of
# the references in Oscar's bibliography, and the filters the library decides
# from precomputed data.
struct LibraryInformation
  title::String
  lines::Vector{String}
  references::Vector{String}
  precomputed_filters::Vector{Any}
end

function LibraryInformation(
  title::String; lines=String[], references=String[], filters=Any[]
)
  return LibraryInformation(title, lines, unique(references), filters)
end

@doc raw"""
    information(L::GroupLibrary)
    information(L::GroupLibrary, k)
    information(L::GroupLibrary, id)

Return information about the data of the library `L`: about the whole
library, about the groups with primary key `k`, or about the group with
identifier `id`.

The result is displayed as an overview that says, for example, how many
groups there are, how they are arranged and which publications they come
from. [`references`](@ref) and [`precomputed_filters`](@ref) return parts of
it.

# Examples
```jldoctest
julia> I = information(transitive_groups_library(), 12)
Transitive groups of degree 12
  301 groups
  decided from precomputed data: order, transitivity, is_primitive
  identify: available
  references:
    GAP package TransGrp [TransGrp]
    Conway, Hulpke, McKay: On transitive permutation groups (1998) [CHM98]
    Royle: The transitive groups of degree twelve (1987) [Roy87]

julia> references(I)
3-element Vector{String}:
 "TransGrp"
 "CHM98"
 "Roy87"
```
"""
function information end

@doc raw"""
    references(I)

Return the keys of the publications behind the data that `I` describes,
where `I` is the result of [`information`](@ref).
The keys refer to the bibliography of OSCAR, the file
`docs/oscar_references.bib` in the OSCAR package.

# Examples
```jldoctest
julia> references(information(perfect_groups_library()))
3-element Vector{String}:
 "PerfGrp"
 "HP89"
 "Hul22"
```
"""
references(I::LibraryInformation) = I.references

@doc raw"""
    precomputed_filters(I)

Return the filters that the library decides from precomputed data for the
groups that `I` describes, where `I` is the result of
[`information`](@ref).
A selection with these filters constructs no group that fails them.

# Examples
```jldoctest
julia> string.(precomputed_filters(information(primitive_groups_library(), 8)))
4-element Vector{String}:
 "order"
 "transitivity"
 "is_simple"
 "is_solvable"
```
"""
precomputed_filters(I::LibraryInformation) = I.precomputed_filters

function Base.show(io::IO, ::MIME"text/plain", I::LibraryInformation)
  io = pretty(io)
  print(io, I.title, Indent())
  for line in I.lines, part in _wrap_line(line)
    print(io, "\n", part)
  end
  if !isempty(I.references)
    print(io, "\nreferences:", Indent())
    for key in I.references, part in _wrap_line(_reference_text(key))
      print(io, "\n", part)
    end
    print(io, Dedent())
  end
  print(io, Dedent())
end

Base.show(io::IO, I::LibraryInformation) = print(io, I.title)

################################################################################
#
#  Helpers for the libraries
#
################################################################################

_gap_package_version(name::String) =
  String(getproperty(GAP.Globals.GAPInfo.PackagesLoaded, Symbol(name))[2])

_availability_line(available::Bool) = available ? "available" : "not available"

function _filters_line(filters::Vector)
  isempty(filters) && return nothing
  return "decided from precomputed data: " * join(string.(filters), ", ")
end

function _count_line(number, name::String)
  number == 1 && return "1 $name"
  return "$number $(name)s"
end

# the ranges of `(from, to, text)` aligned in two columns
function _range_lines(ranges)
  labels = [from == to ? string(from) : "$from - $to" for (from, to, _) in ranges]
  width = maximum(length, labels; init=0)
  return ["  " * rpad(label, width + 2) * text for (label, (_, _, text)) in zip(labels, ranges)]
end

# `line` broken at spaces into parts of at most `width` characters, the parts
# after the first indented two more than `line`; the pretty printer would
# break it inside words
function _wrap_line(line::AbstractString, width::Int=70)
  length(line) <= width && return [line]
  indent = line[1:(findfirst(!isspace, line) - 1)]
  parts = _word_wrap(line, width - length(indent))
  return [i == 1 ? indent * part : indent * "  " * part for (i, part) in enumerate(parts)]
end

# `text` broken into lines of at most `width` characters
function _word_wrap(text::AbstractString, width::Int)
  lines = String[]
  line = ""
  for word in split(text)
    if !isempty(line) && length(line) + 1 + length(word) > width
      push!(lines, line)
      line = ""
    end
    line = isempty(line) ? String(word) : line * " " * word
  end
  isempty(line) || push!(lines, line)
  return lines
end

################################################################################
#
#  The bibliography of Oscar
#
################################################################################

const _BIBLIOGRAPHY = Ref{Union{Nothing,Dict{String,Dict{String,String}}}}(nothing)

# the fields of each entry of `docs/oscar_references.bib`, by key
function _bibliography()
  _BIBLIOGRAPHY[] === nothing || return _BIBLIOGRAPHY[]

  bib = Dict{String,Dict{String,String}}()
  text = read(joinpath(pkgdir(Oscar), "docs", "oscar_references.bib"), String)
  for entry in eachmatch(r"^@\w+\{([^,\s]+),\n(.*?)^\}"ms, text)
    fields = Dict{String,String}()
    for field in eachmatch(r"^\s+(\w+)\s+=\s+\{(.*?)\},?\n"ms, entry[2])
      fields[lowercase(field[1])] = join(split(field[2]), " ")
    end
    bib[entry[1]] = fields
  end
  _BIBLIOGRAPHY[] = bib
  return bib
end

_strip_tex(s::AbstractString) = replace(s, r"[{}$\\]" => "")

# a reference as "authors: title (year) [key]"; a GAP package by its name
function _reference_text(key::String)
  fields = get(_bibliography(), key, nothing)
  fields === nothing && return "[$key]"
  get(fields, "note", "") == "GAP package" && return "GAP package $key [$key]"

  authors = [strip(first(split(a, ","))) for a in split(get(fields, "author", ""), " and ")]
  title = _strip_tex(get(fields, "title", ""))
  year = get(fields, "year", "")
  return "$(join(authors, ", ")): $title ($year) [$key]"
end
