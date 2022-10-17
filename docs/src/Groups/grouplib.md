```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["grouplib.md"]
```

# Group libraries

## Transitive permutation groups of small degree

TODO: explain about the scope of this.

TODO: give proper attribution to the transgrp package (in particular, cite it)

The arrangement and the names of the groups of degree up to 15 is the same as given in
[CHM98](@cite). With the exception of the symmetric and alternating group (which are represented
as `symmetric_group` and `alternating_group`) the generators for these groups also conform to this paper with
the only difference that 0 (which is not permitted in GAP for permutations to act on) is always replaced by
the degree.

The arrangement for all degrees is intended to be equal to the arrangement within the systems GAP and Magma, thus it
should be safe to refer to particular (classes of) groups by their index numbers.


```@docs
all_transitive_groups
has_number_transitive_groups
has_transitive_group_identification
has_transitive_groups
number_transitive_groups
transitive_group
transitive_group_identification
```

## Primitive permutation groups of small degree

TODO: explain about the scope of this.

TODO: give proper attribution to the primitive groups library (in particular, cite it)

```@docs
all_primitive_groups
has_number_primitive_groups
has_primitive_group_identification
has_primitive_groups
number_primitive_groups
primitive_group
primitive_group_identification
```

## Perfect groups of small order

The functions in this section are wrappers for the GAP library of finite perfect
groups which provides, up to isomorphism, a list of all perfect groups whose
sizes are less than $2\cdot 10^6$. The groups of most orders up to $10^6$ have been
enumerated by Derek Holt and Wilhelm Plesken, see [HP89](@cite). For orders
$n = 86016$, 368640, or 737280 this work only counted the groups (but did not
explicitly list them), the groups of orders $n = 61440$, 122880, 172032,
245760, 344064, 491520, 688128, or 983040 were omitted.

Several additional groups omitted from the book [HP89](@cite) have also
been included. Two groups -- one of order 450000 with a factor group of
type $A_6$ and the one of order 962280 -- were found by Jack Schmidt in
2005. Two groups of order 243000 and one each of orders 729000, 871200, 878460
were found in 2020 by Alexander Hulpke.

The perfect groups of size less than $2\cdot 10^6$ which had not been
classified in the work of Holt and Plesken have been enumerated by Alexander
Hulpke. They are stored directly and provide less construction information
in their names.

As all groups are stored by presentations, a permutation representation
is obtained by coset enumeration. Note that some of the library groups do
not have a faithful permutation representation of small degree.
Computations in these groups may be rather time consuming.

```@docs
has_number_perfect_groups
has_perfect_group_identification
has_perfect_groups
number_perfect_groups
orders_perfect_groups
perfect_group
perfect_group_identification
```

## Groups of small order

TODO: explain about the scope of this.

TODO: give proper attribution to the smallgrp package and other things used (in particular, cite it)

```@docs
all_small_groups
has_number_small_groups
has_small_group_identification
has_small_groups
number_small_groups
small_group
small_group_identification
```

## Atlas of Group Representations

The functions in this section give access to data in the
Atlas of Group Representations [WWTSPNNLBA](@cite).
The isomorphism types of the groups in question are specified via
names for the groups, which coincide with the names of the
corresponding character tables in the library of character tables,
see [`character_table(id::String, p::Int = 0)`](@ref).

```@docs
number_atlas_groups
all_atlas_group_infos
atlas_group
atlas_subgroup
```
