# Literature constructions

Certain models have been studied in the physics literature over and over again.
Thereby, these constructions became famous and some were given special names. We
aim to provide support for such standard constructions. An example of such a model is
the following:
```@docs
su5_tate_model_over_arbitrary_3d_base()
```
More generally, we support literature constructions.
```@docs
literature_model(; doi::String="", arxiv_id::String="", version::String="", equation::String="", model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String,Any}(), defining_classes::Dict{String, <:Any} = Dict{String,Any}(), completeness_check::Bool = true)
```

## Attributes

For literature models, we provide the following attributes referencing meta data:
```@docs
arxiv_id(m::AbstractFTheoryModel)
arxiv_doi(m::AbstractFTheoryModel)
arxiv_link(m::AbstractFTheoryModel)
arxiv_model_equation_number(m::AbstractFTheoryModel)
arxiv_model_page(m::AbstractFTheoryModel)
arxiv_model_section(m::AbstractFTheoryModel)
arxiv_version(m::AbstractFTheoryModel)
associated_literature_models(m::AbstractFTheoryModel)
generating_sections(m::AbstractFTheoryModel)
journal_doi(m::AbstractFTheoryModel)
journal_link(m::AbstractFTheoryModel)
journal_model_equation_number(m::AbstractFTheoryModel)
journal_model_page(m::AbstractFTheoryModel)
journal_model_section(m::AbstractFTheoryModel)
journal_name(m::AbstractFTheoryModel)
journal_pages(m::AbstractFTheoryModel)
journal_report_numbers(m::AbstractFTheoryModel)
journal_volume(m::AbstractFTheoryModel)
journal_year(m::AbstractFTheoryModel)
literature_identifier(m::AbstractFTheoryModel)
model_parameters(m::AbstractFTheoryModel)
paper_authors(m::AbstractFTheoryModel)
paper_buzzwords(m::AbstractFTheoryModel)
paper_description(m::AbstractFTheoryModel)
paper_title(m::AbstractFTheoryModel)
birational_literature_models(m::AbstractFTheoryModel)
```
Such meta data can be modified with setters. For instance, there is a function
`set_description(m::AbstractFTheoryModel, description::String)`, which takes the
model in question as the first argument and the desired description - provided as string -
as the second argument. Such a setter function exists for all of the above. If appropriate,
we also offer a method that adds a new value. For instance, we have a function
`add_paper_buzzword(m::AbstractFTheoryModel, addition::String)`.

In addition, the following attributes are available to access advanced model information:
```@docs
resolutions(m::AbstractFTheoryModel)
resolution_generating_sections(m::AbstractFTheoryModel)
resolution_zero_sections(m::AbstractFTheoryModel)
weighted_resolutions(m::AbstractFTheoryModel)
weighted_resolution_generating_sections(m::AbstractFTheoryModel)
weighted_resolution_zero_sections(m::AbstractFTheoryModel)
zero_section(m::AbstractFTheoryModel)
zero_section_class(m::AbstractFTheoryModel)
zero_section_index(m::AbstractFTheoryModel)
exceptional_classes(m::AbstractFTheoryModel)
exceptional_divisor_indices(m::AbstractFTheoryModel)
torsion_sections(m::AbstractFTheoryModel)
```

One can check if a model has a particular set of information. This is achieved with the
following methods:
* `has_arxiv_id(m::AbstractFTheoryModel)`,
* `has_arxiv_doi(m::AbstractFTheoryModel)`,
* `has_arxiv_link(m::AbstractFTheoryModel)`,
* `has_arxiv_model_equation_number(m::AbstractFTheoryModel)`,
* `has_arxiv_model_page(m::AbstractFTheoryModel)`,
* `has_arxiv_model_section(m::AbstractFTheoryModel)`,
* `has_arxiv_version(m::AbstractFTheoryModel)`,
* `has_associated_literature_models(m::AbstractFTheoryModel)`,
* `has_generating_sections(m::AbstractFTheoryModel)`,
* `has_journal_doi(m::AbstractFTheoryModel)`,
* `has_journal_link(m::AbstractFTheoryModel)`,
* `has_journal_model_equation_number(m::AbstractFTheoryModel)`,
* `has_journal_model_page(m::AbstractFTheoryModel)`,
* `has_journal_model_section(m::AbstractFTheoryModel)`,
* `has_journal_name(m::AbstractFTheoryModel)`,
* `has_journal_pages(m::AbstractFTheoryModel)`,
* `has_journal_report_numbers(m::AbstractFTheoryModel)`,
* `has_journal_volume(m::AbstractFTheoryModel)`,
* `has_journal_year(m::AbstractFTheoryModel)`,
* `has_literature_identifier(m::AbstractFTheoryModel)`,
* `has_model_description(m::AbstractFTheoryModel)`,
* `has_model_parameters(m::AbstractFTheoryModel)`,
* `has_paper_authors(m::AbstractFTheoryModel)`,
* `has_paper_buzzwords(m::AbstractFTheoryModel)`,
* `has_paper_description(m::AbstractFTheoryModel)`,
* `has_paper_title(m::AbstractFTheoryModel)`,
* `has_birational_literature_models(m::AbstractFTheoryModel)`,
* `has_resolutions(m::AbstractFTheoryModel)`,
* `has_resolution_generating_sections(m::AbstractFTheoryModel)`,
* `has_resolution_zero_sections(m::AbstractFTheoryModel)`,
* `has_weighted_resolutions(m::AbstractFTheoryModel)`,
* `has_weighted_resolution_generating_sections(m::AbstractFTheoryModel)`,
* `has_weighted_resolution_zero_sections(m::AbstractFTheoryModel)`,
* `has_zero_section(m::AbstractFTheoryModel)`,
* `has_zero_section_class(m::AbstractFTheoryModel)`,
* `has_torsion_sections(m::AbstractFTheoryModel)`,
* `has_gauge_algebra(m::AbstractFTheoryModel)`,
* `has_global_gauge_group_quotient(m::AbstractFTheoryModel)`.


## Methods

### Resolution(s) of a singular model

A central task in F-theory is to resolve a singular model.
For literature models, we have stored resolutions in our data base.
Upon construction of a literature model, we load these known resolutions.

In addition to listing the known resolutions with `resolutions(m::AbstractFTheoryModel)`,
the user might want to add a resolution. This can be achieved with the following method:
```@docs
add_resolution(m::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})
```
Provided that a resolution for a model is known, we can (attempt to) resolve the model.
```@docs
resolve(m::AbstractFTheoryModel, index::Int)
```


## The Quadrillion F-Theory Standard Models

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
