# Using a library of groups

A library is represented by a handle, which is the first argument of the
functions below.
Which groups a library provides, and whether their number and their
identification are available, depends on the library.

| Expression | Meaning |
|:---|:---|
| `L[k, i]` | the group with identifier `(k, i)` |
| `find(L, filters...)` | a lazy selection `S` of groups |
| `for G in S`, `collect(S)`, `first(S)`, `isempty(S)` | the groups in `S` |
| `length(S)` | the number of groups in `S` |
| `keys(S)` | the identifiers of the groups in `S` |
| `identify(L, G)` | the identifier of the group `G` |

```@docs
Oscar.GroupLibrary
find(L::Oscar.GroupLibrary, filters...)
identify
has_groups
has_number_of_groups
has_identification
```

## Available libraries

```@docs
small_groups_library
transitive_groups_library
primitive_groups_library
perfect_groups_library
groups_with_class_number_library
```
