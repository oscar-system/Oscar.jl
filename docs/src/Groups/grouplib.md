```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Group libraries

## Transitive permutation groups of small degree

The functions in this section are wrappers for the GAP library of
transitive permutation groups up to degree 48,
via the GAP package `TransGrp` [Hul23](@cite).

(The groups of degrees 32 and 48 are currently not automatically
available in Oscar,
one has to install additional data in order to access them.)

The arrangement and the names of the groups of degree up to 15 is the same as given in
[CHM98](@cite). With the exception of the symmetric and alternating group (which are represented
as `symmetric_group` and `alternating_group`) the generators for these groups also conform to this paper with
the only difference that 0 (which is not permitted in GAP for permutations to act on) is always replaced by
the degree.

The arrangement for all degrees is intended to be equal to the arrangement within the systems GAP and Magma, thus it
should be safe to refer to particular (classes of) groups by their index numbers.


```@docs
all_transitive_groups
has_number_of_transitive_groups
has_transitive_group_identification
has_transitive_groups
number_of_transitive_groups
transitive_group
transitive_group_identification
```

## Primitive permutation groups of small degree

The functions in this section are wrappers for the GAP library of
primitive permutation groups up to degree 8191,
via the GAP package `PrimGrp` [HRR23](@cite).
See the documentation of this package for more information about
the source of the data.

```@docs
all_primitive_groups
has_number_of_primitive_groups
has_primitive_group_identification
has_primitive_groups
number_of_primitive_groups
primitive_group
primitive_group_identification
```

## Perfect groups of small order

The functions in this section are wrappers for the GAP library of finite perfect
groups which provides, up to isomorphism, a list of all perfect groups whose
sizes are less than $2\cdot 10^6$. The groups of most orders up to $10^6$ have been
enumerated by Derek Holt and Wilhelm Plesken, see [HP89](@cite). For orders
86016, 368640, or 737280 this work only counted the groups (but did not
explicitly list them), the groups of orders 61440, 122880, 172032,
245760, 344064, 491520, 688128, or 983040 were omitted.

Several additional groups omitted from the book [HP89](@cite) have also
been included. Two groups -- one of order 450000 with a factor group of
type $A_6$ and the one of order 962280 -- were found in 2005 by Jack Schmidt.
Two groups of order 243000 and one each of orders 729000, 871200, 878460
were found in 2020 by Alexander Hulpke.

The perfect groups of size less than $2\cdot 10^6$ which had not been
classified in the work of Holt and Plesken have been enumerated by Alexander
Hulpke, see [Hul22](@cite).
They are stored directly and provide less construction information
in their names.

As all groups are stored by presentations, a permutation representation
is obtained by coset enumeration. Note that some of the library groups do
not have a faithful permutation representation of small degree.
Computations in these groups may be rather time consuming.

```@docs
all_perfect_groups
has_number_of_perfect_groups
has_perfect_group_identification
has_perfect_groups
number_of_perfect_groups
orders_perfect_groups
perfect_group
perfect_group_identification
```

## Groups of small order

The functions in this section are wrappers for the GAP library of
the following groups.

The GAP package `SmallGrp` [BEO23](@cite) provides

- those of order at most 2000 (except those of order 1024),
- those of cubefree order at most 50000,
- those of order $p^7$ for the primes $p = 3, 5, 7, 11$,
- those of order $p^n$ for $n \leq 6$ and all primes $p$,
- those of order $q^n p$ where $q^n$ divides $2^8$, $3^6$, $5^5$
  or $7^4$ and $p$ is an arbitrary prime not equal to $q$,
- those of squarefree order,
- those whose order factorises into at most 3 primes.

The GAP package `SOTGrps` [Pan23](@cite) provides

- those whose order factorises into at most 4 primes,
- those of order $p^4 q$ where $p$ and $q$ are distinct primes.

The GAP package `SglPPow` [VE22](@cite)  provides

- those of order $p^7$ for primes $p > 11$,
- those of order $3^8$.

```@docs
all_small_groups
has_number_of_small_groups
has_small_group_identification
has_small_groups
number_of_small_groups
small_group
small_group_identification
```

## Groups with a small number of conjugacy classes

The functions in this section give access to the groups with
up to 14 conjugacy classes.
These groups have been classified in [V-LV-L85](@cite), [V-LV-L86](@cite),
[VS07](@cite).

```@docs
all_groups_with_class_number
has_number_of_groups_with_class_number
has_groups_with_class_number
number_of_groups_with_class_number
group_with_class_number
```

## Atlas of Group Representations

The functions in this section give access to data in the
Atlas of Group Representations [ATLAS](@cite).
The isomorphism types of the groups in question are specified via
names for the groups, which coincide with the names of the
corresponding character tables in the library of character tables,
see [`character_table(id::String, p::Int = 0)`](@ref).

```@docs
number_of_atlas_groups
all_atlas_group_infos
atlas_group
atlas_subgroup
```
