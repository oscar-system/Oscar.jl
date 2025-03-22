```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Group recognition

The idea of constructive group recognition is to compute a *recognition tree*
for a given (permutation or matrix) group, which describes the structure
of this group in a recursive way:
Each non-leaf node of the tree describes an epimorphism such that
the kernel and the image belong to the two subtrees of the node.
Each leaf node describes a group for which efficient methods are available
that allow one to decide whether a group element is an element of this group,
and if yes to write the element as a word in terms of suitable generators.

The recognition tree has enough information to decide whether a group element
is an element of the given group,
and if yes to write the element as a word in terms of suitable generators
of the given group.

```@docs
recognize
is_ready
nice_gens
straight_line_program(tree::GroupRecognitionTree, g::GAPGroupElem)
```
