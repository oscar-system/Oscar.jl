```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Phylogenetic Trees

## Introduction

Phylogenetic trees represent the evolutionary history of some species of consideration.
Here we consider phylogenetic (rooted) trees with branch lengths as defined in [SS03](@cite).

The construction of tropical median consensus trees is one of the nontrivial algorithm here.

## Construction

```@docs
phylogenetic_tree
tropical_median_consensus
```

## Some Helpful Functions

```@docs
adjacency_tree
is_equidistant
cophenetic_matrix
taxa
newick
leaves
n_leaves
```
