```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [The Quadrillion F-Theory Standard Models (QSMs)](@id qsm_models)

A yet more special instance of literature models are the Quadrillion F-theory Standard Models
(F-theory QSMs) [CHLLT19](@cite). Those hypersurface models come in 708 different families.

The base geometry of an F-theory QSM is obtained from triangulating one of 708 reflexive 3-dimensional
polytopes. The models, whose bases are obtained from triangulations of the same polytope form a family.
The following information on the polytope in question and its triangulations is available within our database:

```@docs
vertices(m::AbstractFTheoryModel)
polytope_index(m::AbstractFTheoryModel)
has_quick_triangulation(m::AbstractFTheoryModel)
max_lattice_pts_in_facet(m::AbstractFTheoryModel)
estimated_number_of_triangulations(m::AbstractFTheoryModel)
```

Beyond the polytope and its triangulations, a number of other integers are of key importance. The following
are supported in our database.

```@docs
kbar3(m::AbstractFTheoryModel)
hodge_h11(m::AbstractFTheoryModel)
hodge_h12(m::AbstractFTheoryModel)
hodge_h13(m::AbstractFTheoryModel)
hodge_h22(m::AbstractFTheoryModel)
```

More recently, a research program estimated the exact massless spectra of the F-theory QSMs
(cf. [Bie24](@cite)). These studies require yet more information about the F-theory QSM geometries,
which are supported by our database.

First, recall that (currently), the base of an F-theory QSM is a 3-dimensional toric variety B3.
Let s in H^0(B3, Kbar_B3), then V(s) is a K3-surface. Moreover, let xi be the coordinates
of the Cox ring of B3. Then V(xi) is a divisor in B3. Consequently, Ci = V(xi) cap V(s)
is a divisor in the K3-surface V(s). For the root bundle counting program, these curves Ci are
of ample importance (cf. [Bie24](@cite)). We support the following information on these curves:

```@docs
genera_of_ci_curves(m::AbstractFTheoryModel)
degrees_of_kbar_restrictions_to_ci_curves(m::AbstractFTheoryModel)
topological_intersection_numbers_among_ci_curves(m::AbstractFTheoryModel)
indices_of_trivial_ci_curves(m::AbstractFTheoryModel)
topological_intersection_numbers_among_nontrivial_ci_curves(m::AbstractFTheoryModel)
```

The collection of the Ci-curves form a nodal curve. To every nodal curve one can associate a
(dual) graph. In this graph, every irreducible component of the nodal curve becomes a node/vertex
of the dual graph, and every nodal singularity of the nodal curve turns into an edge of the dual
graph. In the case at hand, this is rather simple.

The Ci-curves turn into the irreducible components of the nodel curve. Certainly, we only need
to focus on the non-trivial Ci-curves. A non-trivial Ci-curve can split into multiple irreducible
components. This is taken into account when the nodes/vertices of the dual graph are constructed.

The topological intersection numbers among the Ci-curves (or rather, their irreducible components)
tells us how many nodal singularities link the Ci-curves (or rather, their irreducible components)
in question. Hence, if the topological intersection numbers is zero, there is no edge between the
corresponding nodes. Otherwise, if the topological intersection number is positive - say n -, then
there are exactly n edges between the nodes in question.

The following functions access/create the so-obtained dual graph:

```@docs
dual_graph(m::AbstractFTheoryModel)
components_of_dual_graph(m::AbstractFTheoryModel)
degrees_of_kbar_restrictions_to_components_of_dual_graph(m::AbstractFTheoryModel)
genera_of_components_of_dual_graph(m::AbstractFTheoryModel)
```

The dual graph is essential in counting root bundles (cf. [BCL21](@cite)). It turns out, that one
can simplify this graph so that the computations at hand can be conducted on a simpler graph
instead. The following functionality exists to access this simplified dual graph.

```@docs
simplified_dual_graph(m::AbstractFTheoryModel)
components_of_simplified_dual_graph(m::AbstractFTheoryModel)
degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph(m::AbstractFTheoryModel)
genera_of_components_of_simplified_dual_graph(m::AbstractFTheoryModel)
```
