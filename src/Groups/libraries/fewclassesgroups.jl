import JSON

const VeraLopezTables = Vector[]

"""
    has_number_of_groups_with_class_number(n::IntegerUnion)

Return whether the number of groups with exactly `n` conjugacy classes is
available for use via [`number_of_groups_with_class_number`](@ref).

# Examples
```jldoctest
julia> has_number_of_groups_with_class_number(14)
true

julia> has_number_of_groups_with_class_number(15)
false
```
"""
has_number_of_groups_with_class_number(n::IntegerUnion) = has_groups_with_class_number(n)

"""
    has_groups_with_class_number(n::IntegerUnion)

Return whether the groups with exactly `n` conjugacy classes are available for
use via [`group_with_class_number`](@ref) and
[`all_groups_with_class_number`](@ref).
This function should be used to test for the scope of the library available.

# Examples
```jldoctest
julia> has_groups_with_class_number(14)
true

julia> has_groups_with_class_number(15)
false
```
"""
function has_groups_with_class_number(n::IntegerUnion)
  @req n >= 1 "class number must be positive, not $n"
  if length(VeraLopezTables) == 0
    # load the data file
    data = JSON.parsefile(joinpath(@__DIR__, "VeraLopez.json"))
    append!(VeraLopezTables, data[2])
  end
  return n <= length(VeraLopezTables)
end

"""
    number_of_groups_with_class_number(n::IntegerUnion)

Return the number of groups with exactly `n` conjugacy classes,
up to isomorphism.

# Examples
```jldoctest
julia> number_of_groups_with_class_number(8)
21

julia> number_of_groups_with_class_number(50)
ERROR: ArgumentError: the number of groups with 50 classes is not available
[...]
```
"""
function number_of_groups_with_class_number(n::IntegerUnion)
  @req has_number_of_groups_with_class_number(n) "the number of groups with $n classes is not available"
  return length(VeraLopezTables[n])
end

"""
    group_with_class_number(::Type{T}, n::IntegerUnion, i::IntegerUnion) where T
    group_with_class_number(n::IntegerUnion, i::IntegerUnion)

Return the `i`-th group in the list of groups with exactly `n`
conjugacy classes.
If a type `T` is specified then an attempt is made to return the result
with that type.
If `T` is omitted then the resulting group will have type `PcGroup`
if it is a solvable group that belongs to the library of small groups,
otherwise it will be of type `PermGroup`.

# Examples
```jldoctest
julia> describe(group_with_class_number(5, 4))
"D14"

julia> describe(group_with_class_number(5, 8))
"A5"

julia> group_with_class_number(5, 12)
ERROR: ArgumentError: there are only 8 groups with 5 classes, not 12
[...]
```
"""
function group_with_class_number(::Type{T}, n::IntegerUnion, i::IntegerUnion) where T
  return T(group_with_class_number(n, i))
end

function group_with_class_number(n::IntegerUnion, i::IntegerUnion)
  @req has_groups_with_class_number(n) "groups with $n classes are not available"
  desc = VeraLopezTables[n]
  @req i <= length(desc) "there are only $(length(desc)) groups with $n classes, not $i"
  return _group_with_class_number(desc[i])
end

function _group_with_class_number(entry::Vector)
  if length(entry) == 2
    return small_group(entry...)
  elseif length(entry) == 3 && entry[1] == "prim"
    return primitive_group(entry[2], entry[3])
  else
    error("not supported format $entry")
  end
end

"""
    all_groups_with_class_number(n::IntegerUnion)

Return a vector of all groups (up to isomorphism)
with exactly `n` conjugacy classes.

# Examples
```jldoctest
julia> map(describe, all_groups_with_class_number(4))
4-element Vector{String}:
 "C4"
 "C2 x C2"
 "D10"
 "A4"
```
"""
function all_groups_with_class_number(n::IntegerUnion)
  @req has_groups_with_class_number(n) "groups with $n classes are not available"
  return [_group_with_class_number(entry) for entry in VeraLopezTables[n]]
end
