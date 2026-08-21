```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Operators for Digraphs

The operations in this chapter mirror Chapter 4 of the GAP Digraphs manual. In particular, digraphs can be compared for equality and ordered, and the subdigraph relations can be tested.

## Operators for digraphs

```@docs
==(D1::Digraph, D2::Digraph)
<(D1::Digraph, D2::Digraph)
is_subdigraph(super::Digraph, sub::Digraph)
is_undirected_spanning_tree(super::Digraph, sub::Digraph)
is_undirected_spanning_forest(super::Digraph, sub::Digraph)
```
