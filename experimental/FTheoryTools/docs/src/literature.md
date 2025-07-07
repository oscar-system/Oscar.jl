```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Literature Models

In the landscape of F-theory model building, many constructions introduced over the years remain relevant
and influential. However, revisiting or building upon these models can be challenging—especially when prior
computations must be repeated manually due to evolving mathematical techniques or notation. This not only
hinders efficiency but can also deter newcomers from engaging deeply with the field.

To overcome these obstacles, we introduce the concept of a **literature model**: a well-established F-theory
construction from the research literature, made accessible in a fully computational format.

---

## What is a Literature Model?

### A Computational Interface to the Literature

The `LiteratureModels` framework provides a curated database of such models, drawn from influential F-theory papers.
Each entry is preprocessed to automate core calculations, store topological and geometric data, and support further
exploration. Rather than repeating extensive algebraic or topological derivations, researchers can immediately load
a literature model and apply the full power of `FTheoryTools` to study or extend it.

This offers several advantages:
- **Rapid access to existing constructions**: Retrieve previously published models without manually reproducing their setup.
- **Integrated results**: Access known data, such as intersection numbers, crepant resolutions, or gauge groups.
- **Targeted search**: Filter models by desired features (e.g., gauge symmetry, fiber type).
- **Seamless analysis**: Load a model and perform further computations using the broader functionality of `FTheoryTools`.

The current database includes a wide range of examples, including models:

- [Krause, Mayrhofer, Weigand 2011](@cite KMW12),
- [Morrison, Park 2012](@cite MP12),
- [Lawrie, Schafer-Nameki 2013](@cite LS13),
- [Klevers, Mayorga, Damian, Oehlmann, Piragua, Reuter 2015](@cite KM-POPR15),
- [Cvetič, Klevers, Piragua, Taylor 2015](@cite CKPT15),
- [Taylor, Wang 2015](@cite TW15),
- [Cvetič, Halverson, Ling, Liu, Tian 2019](@cite CHLLT19).

This collection is actively expanding.

### Interoperable and Reproducible Model Storage

All literature models are stored in a structured and platform-independent format based on JSON. The framework is transitioning
toward the **MaRDI file format**, a standardized data format developed as part of the
[Mathematics Research Data Initiative (MaRDI)](https://www.mardi4nfdi.de). This format ensures that models built in `FTheoryTools`
can be used across a wide range of computer algebra systems, including `Oscar`, `Sage`, `Macaulay`, and others.

The benefits of this include:
- **Interoperability**: Models can be exported, shared, and loaded across different systems.
- **Data compression**: Efficient storage of complex models (via subtree reduction).
- **Machine-readability**: Makes it easy to contribute new models or verify existing ones.
- **Transparency**: Published models can be peer-reviewed in both mathematical and computational terms.

We envision literature models not only as a computational resource, but also as a **new standard for reproducibility** in F-theory
research. Authors can supplement future publications with accompanying data files, making it easier for others to explore, verify,
and build upon their work. This aligns with growing community efforts—such as those within [MaRDI](https://www.mardi4nfdi.de)—to
ensure that mathematical software and data are as verifiable and accessible as the results they support.

---

## Constructing Literature Models

The constructor for literature models is necessarily complex, reflecting the wide range of constructions explored in the F-theory
literature. While some keyword arguments are self-explanatory—for example, the string `doi` uniquely identifies a publication—others
are more technical or use terminology that may be non-standard. A prominent example is `model_sections`, which plays a key role in
defining the geometry of the model.

We defer a detailed discussion of these technical inputs (and outputs) to the next section. The examples in the constructor below serve as
a high-level overview of how literature models are initialized in practice:

```@docs
literature_model(; doi::String="", arxiv_id::String="", version::String="", equation::String="", model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String,Any}(), defining_classes::Dict{String, <:Any} = Dict{String,Any}(), completeness_check::Bool = true)
```

---

## Understanding Defining Classes, Model Parameters and Their Structure

The construction of any F-theory model in `FTheoryTools`—especially those inspired by the literature—relies
on a systematic set of global sections of line bundles over the base of the elliptic fibration. These sections
can be grouped according to their role. For instance, there are **tunable sections**, which play a crucial role
when modifying the singularity structure of the model by tuning it.

Some models also require additional input parameters. For instance, the constructions in [Lawrie, Schafer-Nameki 2013](@cite LS13)
describe a family of global Tate models with gauge group ``\mathrm{SU}(2k)``, where the integer parameter ``k``
must be specified to define a concrete model. Another example appears in [Cvetič, Halverson, Ling, Liu, Tian 2019](@cite CHLLT19),
where models are classified by three-dimensional reflexive polytopes indexed by an integer. A unique base geometry
is then obtained by choosing a fine, regular, star triangulation of the selected polytope. We refer to such structural or discrete
data (e.g., integers fixing singularity types, choices of triangulations, etc.), as model parameters.

In general, the `literature_model` constructor must support a broad spectrum of such inputs from model parameters to divisor classes.

Although most users will engage primarily with the high-level geometry or physics of a model, understanding the underlying
structure of these sections and parameters is essential for advanced tasks such as extending the database, computing fluxes,
or constructing custom geometries. The following sections detail the core components of this infrastructure.

---

### Defining Classes

The **defining classes** of an F-theory model constitute the minimal set of divisor classes used to disambiguate the model
from its most generic form, i.e. generic Weierstrass, Tate, or hypersurface model.

As an example, consider the global Tate model introduced in [Krause, Mayrhofer, Weigand 2011](@cite KMW12). In that construction,
the Tate coefficients are factored as follows:

```julia
a1 = a10 * w
a2 = a21 * w
a3 = a32 * w^2
a4 = a43 * w^3
a6 = 0
```

Here, ``w`` is a section of a line bundle associated to a divisor class ``W`` on the base. The class ``W`` is the sole defining
class for this model—it uniquely determines how all the Tate sections are twisted.

The `literature_model` constructor accepts an optional `defining_classes` argument, which should be a dictionary mapping
parameter names (such as `"w"`) to their corresponding divisor classes. To enable this, the base space must also be provided
explicitly and must be a concrete variety—typically a toric variety.

The following example demonstrates how to construct the model of [Krause, Mayrhofer, Weigand 2011](@cite KMW12) over a toric base
with a chosen defining class ``W``:

```@docs
defining_classes(::AbstractFTheoryModel)
```

---

### Tunable Sections

The hypersurface equation of an F-theory model—whether a Weierstrass, Tate, or hypersurface model—is built from a number of global
sections of line bundles over the base. These sections fall into two conceptual categories:

- **Structural sections**, such as the Weierstrass or Tate coefficients, are determined by the type of model being constructed. To preserve the model's nature—e.g., that it remains a Weierstrass or global Tate model—these sections must retain their prescribed divisor classes. Altering them would fundamentally change the class of the model.

- **Tunable sections** are those sections in the defining hypersurface equation that can be freely varied or further factorized, without changing the qualitative structure of the model. In particular, sections of the defining classes are tunable sections.

For example, in the global Tate model from [Krause, Mayrhofer, Weigand 2011](@cite KMW12), the Tate coefficients are factorized as:

```julia
a1 = a10 * w
a2 = a21 * w
a3 = a32 * w^2
a4 = a43 * w^3
a6 = 0
```

This model is, by definition, a global Tate model. To preserve this, the classes of the `a1`, `a2`, etc., must remain fixed. However,
the section ``w`` of the defining class ``W`` as well as the auxiliary sections ``aij`` introduced via the factorization can be adjusted
freely. Therefore the tunable sections of this model are `w`, `a10`, `a21`, `a32` and `a43`.

```@docs
tunable_sections(::AbstractFTheoryModel)
```

---

### Classes of Tunable Sections

The function `classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes` returns a matrix encoding the divisor
classes of all **tunable sections** of an F-theory model. These classes are expressed in a basis formed by the
anticanonical class ``\overline{K}_B`` of the base space and the defining classes of the model.

Each **column** of the returned matrix corresponds to a tunable section, and each **row** corresponds to a basis element—i.e.,
the first row gives the coefficient with respect to ``\overline{K}_B``, and the subsequent rows represent contributions
from the defining classes.

For an example, consider the global Tate model presented in [Krause, Mayrhofer, Weigand 2011](@cite KMW12), where the
Tate sections are factorized as follows:

```julia
a1 = a10 * w
a2 = a21 * w
a3 = a32 * w^2
a4 = a43 * w^3
a6 = 0
```

The tunable sections in this model are `w`, `a10`, `a21`, `a32` and `a43`. The only defining class is `W`. Therefore, the
divisor classes of the tunable sections in the basis (``\overline{K}_B``, `W`) are given by the matrix:

$$
\begin{bmatrix}
0 & 1 & 2 & 3 & 4 \\
1 & 0 & -1 & -2 & -3
\end{bmatrix}
$$

Each column in this matrix corresponds to one of the tunable sections, listed in the same order as above.

```@docs
classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes(::AbstractFTheoryModel)
```

---

### Model Section Parametrization

```@docs
model_section_parametrization(::AbstractFTheoryModel)
```

This dictionary explains how the standard sections used in the model (e.g., `a2`, `a3`, `a4`, etc.) are written in terms of the tunable sections.

Using the earlier example:

```julia
"a2" => a21 * w
```

means that the section `a2` is not an independent section—it is constructed as a product of tunable ones. This information is key to understanding how the model was derived from a more generic form.

---

### Model Sections

```@docs
model_sections(::AbstractFTheoryModel)
```

This function returns the complete list of all named sections that appear in the model. It includes:
- The tunable sections (like `w`, `a21`, `a32`, ...)
- Any derived sections that are constructed from them (like `a2`, `a3`, ...)

It is the union of the keys from `model_section_parametrization` and `tunable_sections`.

---

### Explicit Model Sections

```@docs
explicit_model_sections(::AbstractFTheoryModel)
```

This dictionary gives the actual expressions for all model sections—written as polynomials in the coordinates of the base.

Each entry in the dictionary corresponds to a name returned by `model_sections`.

---

### Classes of Model Sections

```@docs
classes_of_model_sections(::AbstractFTheoryModel)
```

This dictionary assigns a divisor class to every model section in the model. Again, the keys match those in `model_sections`.

---

### Hypersurface Equation Parametrization

```@docs
hypersurface_equation_parametrization(::HypersurfaceModel)
```

This function ties all the above together by giving the hypersurface equation that defines the geometry of the model, written in terms of the symbolic model sections.

---

Together, these tools provide fine-grained control over the algebraic structure of F-theory models and allow users to trace back each part of the geometry to its symbolic origin.














---

## Meta Data Attributes for Liteature Models

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

!!! warning
    Calling `set_description(m::AbstractFTheoryModel, description::String)` overwrites the existing description. This applies similarly to all other setter functions. Use with care, as existing data will be replaced without warning.

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





## Attributes for Liteature Models

In addition, the following attributes are available to access advanced model information:
```@docs
resolutions(m::AbstractFTheoryModel)
resolution_generating_sections(m::AbstractFTheoryModel)
resolution_zero_sections(m::AbstractFTheoryModel)
torsion_sections(m::AbstractFTheoryModel)
weighted_resolutions(m::AbstractFTheoryModel)
weighted_resolution_generating_sections(m::AbstractFTheoryModel)
weighted_resolution_zero_sections(m::AbstractFTheoryModel)
zero_section(m::AbstractFTheoryModel)
zero_section_class(m::AbstractFTheoryModel)
zero_section_index(m::AbstractFTheoryModel)
exceptional_classes(m::AbstractFTheoryModel)
exceptional_divisor_indices(m::AbstractFTheoryModel)
```




* `has_resolutions(m::AbstractFTheoryModel)`,
* `has_resolution_generating_sections(m::AbstractFTheoryModel)`,
* `has_resolution_zero_sections(m::AbstractFTheoryModel)`,
* `has_torsion_sections(m::AbstractFTheoryModel)`,
* `has_weighted_resolutions(m::AbstractFTheoryModel)`,
* `has_weighted_resolution_generating_sections(m::AbstractFTheoryModel)`,
* `has_weighted_resolution_zero_sections(m::AbstractFTheoryModel)`,
* `has_zero_section(m::AbstractFTheoryModel)`,
* `has_zero_section_class(m::AbstractFTheoryModel)`,
* `has_zero_section_index(m::AbstractFTheoryModel)`, DOES THAT WORK????

For the following, the corresponding methods are not listed above.
Do they exist? Do they have a doc string?

* `has_gauge_algebra(m::AbstractFTheoryModel)`,
* `has_global_gauge_group_quotient(m::AbstractFTheoryModel)`.
* `has_birational_literature_models(m::AbstractFTheoryModel)`,




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
