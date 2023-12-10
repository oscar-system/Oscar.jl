# Set-Partitions

## Aims

The goal of this project is to
bring set-partitions to OSCAR.
These are usually depicted as string diagrams and can be used to define partition algebras like the Temperley–Lieb algebra. 
More generally, one can use set-partitions to construct tensor categories, which appear for example as representation categories of so-called easy quantum groups [BS09](@cite).

## Status

We implemented data structures for partitions and corresponding algorithms as described in [Vol23](@cite). 
This includes:
* basic set-partitions and operations on them (e.g. composition, tensor product, involution)
* variations like colored partitions [TW18](@cite) and spatial partitions [CW16](@cite)
* enumeration of partitions which can be constructed from a set of generators

In the future, we plan to implement:
* linear combinations of partitions
* examples of tensor categories in the framework of [TensorCategories.jl](https://github.com/FabianMaeurer/TensorCategories.jl)
