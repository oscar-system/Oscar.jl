```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Group libraries (variant 1)

Oscar provides access to several libraries of groups.
There is a Julia module for each such library,
the functions defined in this module are typically as follows,
see the library of [Transitive permutation groups of small degree](@ref transgrp_section)
for examples.
(Not all of these functions are available for all group libraries,
and some libraries provide also other functions.)

- `get` returns the group from the library that is defined by the arguments,
  if the library contains this group
  (which can be checked via the `has` function)
- `get_all` returns all groups from the library that have the properties
  given by the arguments,
  if the library contains these groups;
  it depends on the library which equivalence relation is expressed by
  "all" (up to abstract isomorphism, up to permutation isomorphism, etc.)
- `get_one` returns one group from the library that has the properties
  given by the arguments if the library contains the groups in question
  and if such a group exists, and `nothing` otherwise
- `has` returns `true` if the library contains all groups
  given by the arguments
- `has_count` returns `true` if the number of groups given by the arguments
  is known, and `false` otherwise
- `count` returns the number of groups given by the arguments,
  if the library contains this information
- `identification` takes a group `G` and returns a value
  such that calling `get` with this value yields a group from the library
  that corresponds to `G`,
  if the library contains this information
- `has_identification` returns `true` if there is such an `identification`
  function, and `false` otherwise

## [Transitive permutation groups of small degree](@id transgrp_section)

```@docs
TransitiveGroups
```

```@docs
TransitiveGroups.get
TransitiveGroups.get_all
TransitiveGroups.has
TransitiveGroups.has_count
TransitiveGroups.count
TransitiveGroups.has_identification
TransitiveGroups.identification
```

## Perfect groups of small order

```@docs
PerfectGroups
```

```@docs
PerfectGroups.get
PerfectGroups.get_all
PerfectGroups.has
PerfectGroups.has_count
PerfectGroups.count
PerfectGroups.has_identification
PerfectGroups.identification
PerfectGroups.orders
```
