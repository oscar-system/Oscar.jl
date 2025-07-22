```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [The F-Theory QSMs](@id qsm_models)

The **Quadrillion F-Theory Standard Models (QSMs)** were introduced in [CHLLT19](@cite)
as a large class of compactifications within F-theory that exhibit promising features for
particle physics model building. Due to their richness and physical relevance, these models
have been extensively studied; see for example [Bie24](@cite) and the references therein.

The detailed geometric and physical properties of the QSMs are captured within `FTheoryTools`
through the framework of [Literature Models](@ref literature_models), enabling computations
and explorations of their characteristics.

---

## What are the QSM?

The QSMs correspond to a vast ensemble of F-theory compactifications constructed on elliptically
fibered Calabi–Yau fourfolds with carefully chosen base geometries. Their defining feature is that
they yield low-energy effective theories closely resembling the Standard Model of particle physics,
including realistic gauge groups and chiral matter spectra.

The geometries arise from systematically classifying certain favorable base threefolds, often toric,
and then carefully engineering elliptic fibrations that produce the Standard Model gauge group and
chiral matter spectrum.

In `FTheoryTools`, we focus on the QSMs built over toric base spaces. It turned out that there are
exactly 708 families of such QSM geometries. Each family corresponds to the fine, regular, star
triangulation of one reflexive 3-dimensional polytope. Since the number of such triangulations can
be astronomical — more specifically, the largest such number is a quadrillion (``10^{15}``), which gives
this family of F-theory models its name — in `FTheoryTools` one such triangulation is chosen. In
other words, we could read this as a chosen representative for each family, that can be constructed
via the framework of [Literature Models](@ref literature_models).

---

## Constructing QSMs

The construction of QSMs happens via the literature model constructor. Here is an example:

```julia
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base
```

The value assigned to the model parameter `k` must be altered to construct different QSMs.
Recall that there are exactly 708 distinct classes of QSMs, each family corresponding to the
fine, regular, star triangulations of the k-th reflexive polytope in the Kreuzer-Skarke list
(cf. [KS98](@cite)).

Accordingly, `k` corresponds exactly to the index of the polytopes as they appear in the
Kreuzer-Skarke list. The "first" QSM is associated with the 4th polytope, and the last one with
index `k = 4291`. Between these, 706 other values correspond to QSMs. If a given `k` does not
correspond to a QSM, the constructor raises an error.

---

## [The Underlying Polytope](@id qsm_polytope)

Each QSM corresponds to a family of models associated with the fine regular star triangulations
of a reflexive 3-dimensional polytope. The following methods provide detailed information about
the underlying polytope and its triangulations.

### Vertices

Return the vertices of the underlying polytope. Vertices are normalized to rational numbers following the Polymake standard.

```@docs
vertices(m::AbstractFTheoryModel)
```

### Polytope Index

Return the index of the underlying polytope within the Kreuzer-Skarke list ([KS98](@cite)).

```@docs
polytope_index(m::AbstractFTheoryModel)
```

### Quick Triangulation Feasibility

Enumerating all full, regular, star triangulations of a 3-dimensional reflexive polytope can be extremely time-consuming. This method returns `true` if enumerating all such triangulations is expected to complete in a reasonable time (e.g., within 5 minutes on a personal computer), and `false` otherwise.

```@docs
has_quick_triangulation(m::AbstractFTheoryModel)
```

### Maximal Number of Lattice Points in a Facet

To enumerate all full, regular, star triangulations of the underlying polytope, one can first work out the corresponding triangulations of all facets [HT17](@cite). A rough estimate for the complexity of finding all triangulations of a facet is its number of lattice points. This method returns the maximal number of lattice points over all facets of the underlying polytope.

```@docs
max_lattice_pts_in_facet(m::AbstractFTheoryModel)
```

### Estimated Number of Triangulations

Provide an estimate of the total number of full, regular, star triangulations of the underlying polytope.

```@docs
estimated_number_of_triangulations(m::AbstractFTheoryModel)
```

---

## [Topological Data of a QSM](@id qsm_top_data)

In addition to the polytope data, several topological invariants play a central role in analyzing
the geometry of a QSM.

The following method returns the triple self-intersection number of the anticanonical class
``\overline{K}_{B_3}`` of the 3-dimensional base:

```@docs
kbar3(m::AbstractFTheoryModel)
```

The Hodge numbers of the elliptically fibered Calabi–Yau 4-fold are also available. For details, see
[Functionality for all F-theory models](@ref functionality_for_all_f_theory_models).

---

## [The Nodal Curve](@id qsm_nodal_curve)

Recent studies (see [Bie24](@cite)) have refined the exact massless spectrum estimates for F-theory QSMs.
These developments rely on deeper geometric data, much of which is accessible via the `FTheoryTools` database.

Recall that for `FTheoryTools`, the base of a QSM is a 3-dimensional toric variety ``B_3``. Let
``s \in H^0(B_3, \bar{K}_{B_3})`` be a generic section of the anticanonical bundle. Then ``V(s) \subset B_3``
defines a K3 surface. Additionally, each homogeneous coordinate ``x_i`` of the Cox ring of ``B_3`` defines
a divisor ``V(x_i) \subset B_3``, and we define the curve:

$C_i := V(x_i) \cap V(s)\,.$

Collectively, the curves ``C_i`` form a nodal curve, which plays a crucial role in the exact massless spectrum
estimates (see [Bie24](@cite) and references therein).

The following quantities associated to the curves ``C_i`` are supported.

### Genera of the Curves ``C_i``

Return the genera of all ``C_i = V(x_i, s)``, labeled by their respective Cox ring coordinate indices.

```@docs
genera_of_ci_curves(m::AbstractFTheoryModel)
```

### Degrees of ``\bar{K}_{B_3}|_{C_i}``

Compute the degree of the restriction of the anticanonical bundle ``\bar{K}_{B_3}`` to each curve ``C_i``.

```@docs
degrees_of_kbar_restrictions_to_ci_curves(m::AbstractFTheoryModel)
```

### Topological Intersection Numbers Among All ``C_i``

Returns the intersection matrix among the ``C_i`` curves.

```@docs
topological_intersection_numbers_among_ci_curves(m::AbstractFTheoryModel)
```

### Indices of Trivial Curves

Some intersections ``V(x_i, s)`` may be empty. This function returns the list of all indices ``i`` such that ``C_i = \emptyset``.

```@docs
indices_of_trivial_ci_curves(m::AbstractFTheoryModel)
```

### Intersection Matrix of Non-trivial Curves

Returns the intersection matrix among the non-trivial ``C_i``.

```@docs
topological_intersection_numbers_among_nontrivial_ci_curves(m::AbstractFTheoryModel)
```

---

## [The Dual Graph](@id qsm_dual_graph)

The collection of ``C_i``-curves described above forms a nodal curve. To such a curve, one can associate a **dual graph**:

- Each **irreducible component** of the nodal curve corresponds to a **vertex** (or node) in the graph.
- Each **nodal intersection** (i.e. a singular point where two components meet) becomes an **edge**.

Only **non-trivial** ``C_i = V(x_i, s)`` curves are included in this graph. If a curve ``C_i`` is reducible,
its irreducible components each get their own vertex. The number of edges between any two vertices is determined
by the **topological intersection number** of their corresponding components.

The following methods allow access to this graph and its associated data.

### The Dual Graph

Returns the unlabelled, undirected dual graph for the given QSM. Each node corresponds to an irreducible component of a non-trivial ``C_i`` curve.

```@docs
dual_graph(m::AbstractFTheoryModel)
```

### Component Labels

Returns the labels — each a string — for each node in the dual graph, indicating the geometric origin of the corresponding irreducible component. If ``C_i`` is irreducible, the label is `"Ci"`. If ``C_i`` is reducible, its components are labeled as `"Ci-0"`, `"Ci-1"`, etc.

```@docs
components_of_dual_graph(m::AbstractFTheoryModel)
```

### Degrees of the restriction of ``\bar{K}_{B_3}`` to the Components

Computes the degree of the anticanonical bundle ``\bar{K}_{B_3}`` when restricted to each component (node) of the dual graph.

```@docs
degrees_of_kbar_restrictions_to_components_of_dual_graph(m::AbstractFTheoryModel)
```

### Genera of the Components

Returns the genus of each irreducible component represented by a node in the dual graph.

```@docs
genera_of_components_of_dual_graph(m::AbstractFTheoryModel)
```

---

## [The Simplified Dual Graph](@id qsm_simple_dual_graph)

The estimates for the exact massless spectrum of the QSMs are based on root bundle counting
(cf. [BCL21](@cite)). This in turn employs the **dual graph** described in the previous section.
However, it turns out, that the **dual graph** can often be simplified without loss of relevant
information for the root bundle counting. This leads to a **reduced** version of the graph on
which the root bundle computations are easier to perform. The following functions provide access
to this **simplified dual graph**.

```@docs
simplified_dual_graph(m::AbstractFTheoryModel)
```

### Component Labels

Returns the labels — each a string — for each vertex in the simplified dual graph. These labels correspond to irreducible components of non-trivial ``C_i = V(x_i, s)`` curves.

```@docs
components_of_simplified_dual_graph(m::AbstractFTheoryModel)
```

### Degrees of ``\bar{K}_{B_3}`` on the Components

Returns the degrees of the anticanonical bundle ``\bar{K}_{B_3}`` restricted to the components (nodes) of the simplified dual graph.

```@docs
degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph(m::AbstractFTheoryModel)
```

### Genera of the Components

Returns the genus of each irreducible component corresponding to a node in the simplified dual graph.

```@docs
genera_of_components_of_simplified_dual_graph(m::AbstractFTheoryModel)
```
