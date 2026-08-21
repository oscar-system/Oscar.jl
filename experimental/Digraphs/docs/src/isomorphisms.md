```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Isomorphisms and Homomorphisms

The operations in this chapter mirror Chapter 7 of the GAP Digraphs manual: acting on digraphs, isomorphisms and canonical labellings, and graph homomorphisms.

## Acting on digraphs

```@docs
on_digraphs(D::Digraph, p::PermGroupElem)
Base.:^(D::Digraph, p::PermGroupElem)
on_multi_digraphs(D::Digraph, pair)
on_tuples_digraphs
on_sets_digraphs
```

## Isomorphisms and canonical labellings

```@docs
digraphs_use_nauty()
digraphs_use_bliss()
automorphism_group(D::Digraph)
bliss_automorphism_group(D::Digraph)
nauty_automorphism_group(D::Digraph)
bliss_canonical_labelling(D::Digraph)
nauty_canonical_labelling(D::Digraph)
bliss_canonical_digraph(D::Digraph)
nauty_canonical_digraph(D::Digraph)
digraph_group(D::Digraph)
digraph_orbits(D::Digraph)
digraph_orbit_reps(D::Digraph)
digraph_schreier_vector(D::Digraph)
digraph_stabilizer(D::Digraph, v::Int)
is_isomorphic(D1::Digraph, D2::Digraph)
isomorphism_digraphs(D1::Digraph, D2::Digraph)
representative_out_neighbours(D::Digraph)
is_digraph_isomorphism
is_digraph_automorphism
is_digraph_colouring(D::Digraph, list::Vector{<:Integer})
maximal_common_subdigraph(D1::Digraph, D2::Digraph)
minimal_common_superdigraph(D1::Digraph, D2::Digraph)
```

## Homomorphisms of digraphs

```@docs
homomorphism_digraphs_finder(D1::Digraph, D2::Digraph, hook, user_param, max_results, hint, injective, image, partial_map, colors1, colors2; order=nothing, aut_grp=nothing)
digraph_homomorphism(D1::Digraph, D2::Digraph)
homomorphisms_digraphs(D1::Digraph, D2::Digraph)
homomorphisms_digraphs_representatives(D1::Digraph, D2::Digraph)
digraph_monomorphism(D1::Digraph, D2::Digraph)
monomorphisms_digraphs(D1::Digraph, D2::Digraph)
monomorphisms_digraphs_representatives(D1::Digraph, D2::Digraph)
digraph_epimorphism(D1::Digraph, D2::Digraph)
epimorphisms_digraphs(D1::Digraph, D2::Digraph)
epimorphisms_digraphs_representatives(D1::Digraph, D2::Digraph)
digraph_embedding(D1::Digraph, D2::Digraph)
embeddings_digraphs(D1::Digraph, D2::Digraph)
embeddings_digraphs_representatives(D1::Digraph, D2::Digraph)
is_digraph_homomorphism
is_digraph_epimorphism
is_digraph_monomorphism
is_digraph_endomorphism
is_digraph_embedding
subdigraphs_monomorphisms(D1::Digraph, D2::Digraph)
subdigraphs_monomorphisms_representatives(D1::Digraph, D2::Digraph)
digraphs_respects_colouring(src::Digraph, ran::Digraph, x, col1, col2)
generators_of_endomorphism_monoid(D::Digraph)
generators_of_endomorphism_monoid_attr(D::Digraph)
digraph_colouring(D::Digraph, n::Int)
digraph_greedy_colouring(D::Digraph)
digraph_welsh_powell_order(D::Digraph)
chromatic_number(D::Digraph)
digraph_core(D::Digraph)
lattice_digraph_embedding(L1::Digraph, L2::Digraph)
is_lattice_homomorphism
is_lattice_epimorphism
is_lattice_embedding
is_lattice_monomorphism
is_lattice_endomorphism
```
