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

```@docs
number_transitive_groups
transitive_group
transitive_identification
all_transitive_groups
```

## Primitive permutation groups of small degree

```@docs
number_primitive_groups
primitive_group
all_primitive_groups
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
orders_perfect_groups
number_perfect_groups
perfect_group
perfect_group_identification
```

## Groups of small order

```@docs
number_small_groups
small_group
small_group_identification
all_small_groups
```
