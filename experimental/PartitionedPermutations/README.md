# Partitioned Permutations

This project implements algorithms for working with partitioned permutations. These are used in the context of Free Probability Theory to model higher order freeness and higher order free cumulants. The central reference is the paper [Second Order Freeness and Fluctuations of Random Matrices III. Higher order freeness and free cumulants (2007)](https://arxiv.org/pdf/math/0606431.pdf).

## Mathematical Background

A partitioned permutation $(V, \pi)$ consists of a permutation $\pi$ and a partition $V$ of the set of the set {1, ..., n\} such that the partition dominates the permutation in the sense that every cycle of $\pi$ is contained in one block of $V$. One sets
$$|(V, \pi)| := n - ( 2 \cdot \text{number of blocks of } V - \text{number of cycles of } \pi).$$

For two partitioned permutations $(V, \pi)$ and $(W, \sigma)$ one defines their product as
$$(V, \pi) \cdot (W, \sigma) = (V \vee W, \pi \sigma)$$
if $|(V, \pi)| + |(W, \sigma)| = |(V \vee W, \pi \sigma)|$. Otherwise one sets $(V, \pi) \cdot (W, \sigma) = (O, \mathrm{id})$. Here, $O$ is the partition where every block consists of exactly one element.

## Status

We implemeted the type `PartitionedPermutation` together with the following methods.
- a function `join` for computing the join of two set partitions (presented as objects of type `SetPartition`
- a comparison `<=` for set partitions
- a function `cycle_partition` that returns the cycle partition of a permutation as set partition
- functions `length` `length2` for computing the number $n$ of underlying elements of a partitioned permutation $(V, \pi)$ and the number $|(V, \pi)|$, respectively
- a function that enumerates all partitioned permutations of a fixed length $n$
- a function `*` that returns the product of two partitioned permutations

We plan to implement
- a function that computes all factorizations of a given partitioned permutation

## Contact

Please direct questions about this part of OSCAR to the following people:
* Björn Schäfer
* Sebastian Volz

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).
