```@meta
CurrentModule = Oscar
```

# Partially Ordered Sets

## Introduction

Partially ordered sets (a.k.a. *posets*) describe a finite set of elements with a partial ordering.
All such objects in OSCAR need to have a rank function but it is not required
that it is graded, i.e., adjacent nodes may have a rank difference larger than one.
There must also be a unique minimal element in the set.

Internally, posets are encoded in terms of their *Hasse diagrams*.
The latter object is the directed graph whose nodes are the elements, and the oriented edges are given by the covering relations.
Usually, the direction is from bottom to top; but there are exceptions, e.g., due to lazy evaluation.

Posets are static in the sense that they are born once, and then they remain immutable.
This is different from [Graphs](@ref).

!!! note
    A unique maximal element is also required but in some cases this is realized by adding an
    extra artificial top node which covers all maximal elements. This node is suppressed in most
    functions except for the visualization.


## Construction

```@docs
partially_ordered_set(covrels::Matrix{Int})
partially_ordered_set(g::Graph{Directed})
partially_ordered_set(g::Graph{Directed}, node_ranks::Dict{Int,Int})
partially_ordered_set_from_inclusions(I::IncidenceMatrix)
face_poset(p::Union{Polyhedron,Cone,PolyhedralFan,PolyhedralComplex,SimplicialComplex})
maximal_ranked_poset(v::AbstractVector{Int})
lattice_of_flats(m::Matroid)
lattice_of_cyclic_flats(m::Matroid)
```


## Auxiliary functions

```@docs
rank(p::PartiallyOrderedSet)
comparability_graph(pos::PartiallyOrderedSet)
maximal_chains(p::PartiallyOrderedSet)
order_polytope(p::PartiallyOrderedSet)
chain_polytope(p::PartiallyOrderedSet)
graph(p::PartiallyOrderedSet)
rank(pe::PartiallyOrderedSetElement)
node_ranks(p::PartiallyOrderedSet)
n_atoms(p::PartiallyOrderedSet)
n_coatoms(p::PartiallyOrderedSet)
atoms(p::PartiallyOrderedSet)
coatoms(p::PartiallyOrderedSet)
least_element(p::PartiallyOrderedSet)
greatest_element(p::PartiallyOrderedSet)
element(p::PartiallyOrderedSet, i::Int)
elements(p::PartiallyOrderedSet)
elements_of_rank(p::PartiallyOrderedSet, rk::Int)
```

## Visualization

```@docs
visualize(p::PartiallyOrderedSet; AtomLabels=[], filename=nothing)
```
