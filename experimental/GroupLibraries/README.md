# Group libraries behind one interface

## Aims

Each library of groups in OSCAR (small groups, transitive groups, primitive
groups, perfect groups, groups with few conjugacy classes) has its own set of
functions, e.g. `small_group`, `all_small_groups`, `number_of_small_groups`.
This project offers the same libraries through one interface: a handle
object per library, and generic functions that take the handle.
Selections are lazy, and counting the groups with given properties does not
construct them where the library can avoid it.

See oscar-system/Oscar.jl#1168 for the discussion that led to this design.

## Status

The names and the behaviour of the functions may still change.
The functions `small_group` etc. are not affected.
