```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Partitioned Permutations

Partitioned Permutations are used in the context of Free Probability Theory to model higher order freeness and higher order free cumulants, see e.g. [CMS07](@cite).

We provide basic functions for working with partitioned permutations. The focus is on factorizing partitioned permutations.

# Basics

Formally, a partitioned permutation $(V, \pi)$ consists of a permutation $\pi$ and a partition $V$ of the set $\{1, ..., n\}$ such that the partition dominates the permutation in the sense that every cycle of $\pi$ is contained in one block of $V$. We call $n$ the length of $(V, \pi)$. Mathematically, another length function is more important. It is given by
$$|(V, \pi)| := n - ( 2 \cdot \text{number of blocks of } V - \text{number of cycles of } \pi),$$
and we call this the adjusted length of $(V, \pi)$. Note that this terminology is not used in the literature.

```@docs
PartitionedPermutation
length(pp::PartitionedPermutation)
adjusted_length(pp::PartitionedPermutation)
```

# Products of Partitioned Permutations

For two partitioned permutations $(V, \pi)$ and $(W, \sigma)$ one defines their product as
$$(V, \pi) \cdot (W, \sigma) = (V \vee W, \pi \sigma)$$
if $|(V, \pi)| + |(W, \sigma)| = |(V \vee W, \pi \sigma)|$. Otherwise, one sets $(V, \pi) \cdot (W, \sigma) = (O, \mathrm{id})$. Here, $O$ is the partition where every block consists of exactly one element, and $V \vee W$ denotes the join of the partitions $V$ and $W$.

A major problem is the factorization of a partitioned permutation $(V, \pi)$. This involves finding all pairs $(W_1, \sigma_1), (W_2, \sigma_2)$ of partitioned permutations with $(V, \pi) = (W_1, \sigma_1) \cdot (W_2, \sigma_2)$.

```@docs
*(pp_1::PartitionedPermutation, pp_2::PartitionedPermutation)
enumerate_partitioned_permutations(n::Int)
factor(pp::PartitionedPermutation)
```
