```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Literature Models](@id literature_models)

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

## [Model Parameters, Defining Classes, and Related Concepts](@id model_section_explanation_section)

### General Overview

The `literature_model` constructor in `FTheoryTools` is designed to support a broad spectrum of input parameters,
ranging from discrete numerical values to divisor classes. While most users will interact primarily with the geometric
or physical properties of a model, understanding its internal structure is essential for advanced applications such
as extending the model database, computing fluxes, or constructing custom geometries.

Some literature models require discrete input parameters. For instance, the constructions in [Lawrie, Schafer-Nameki 2013](@cite LS13)
describe a family of global Tate models with gauge group ``\mathrm{SU}(2k)``, where the integer parameter ``k`` must be
specified to define a concrete model. Similarly, in [Cvetič, Halverson, Ling, Liu, Tian 2019](@cite CHLLT19), models are
indexed by three-dimensional reflexive polytopes, each identified by an integer. A specific base geometry is selected by
choosing a fine, regular, star triangulation of the chosen polytope. We refer to such discrete data—e.g., integers fixing
singularity types, triangulation choices—as **model parameters**.

Let us now focus on the hypersurface equation of a model. Regardless of whether the model is a Weierstrass, global Tate, or
hypersurface model, the hypersurface equation is constructed from global sections of line bundles over the base. These sections
fall into two conceptual categories:

- **Structural sections**, such as the Weierstrass or Tate coefficients, are determined by the type of model being constructed. To preserve the model's nature, these sections must retain their prescribed divisor classes. Altering them would fundamentally change the type of the model, e.g. turning a global Tate model into a Weierstrass model.

- **Tunable sections**: These are sections in the defining equation that can be freely varied or further factorized without altering the qualitative structure of the model. In other words, they preserve the type of the model.

These two categories are not mutually exclusive. Unless factorized at construction time, a structural section can also be a tunable section.

Together, structural and tunable sections form the set of **model sections**. Since each model section is a global section of a line bundle,
it has an associated divisor class. However, not all of these classes are needed to construct the model. Only a subset—the classes of certain
tunable sections—are generally required. These are called the **defining classes**.

The defining classes play a central role in specifying literature models. The `literature_model` constructor accepts an optional `defining_classes`
argument, which should be a dictionary mapping class symbols to their corresponding divisor classes. Note that because divisor classes are
only supported over concrete base spaces, this assignment can only be made for explicitly chosen bases. For unspecified base families, the
model is constructed to be compatible with a more general divisor class structure.

In summary, the **defining classes** form the minimal set of divisor classes needed to distinguish a given F-theory model from its most generic
Weierstrass, Tate, or hypersurface formulation.

#### Example: Global Tate Model

Consider the global Tate model introduced in [Krause, Mayrhofer, Weigand 2011](@cite KMW12). The Tate coefficients are factorized as:

```julia
a1 = a10

a2 = a21 * w

a3 = a32 * w^2

a4 = a43 * w^3

a6 = 0
```

From this, we find:

- **Structural sections**: `a1`, `a2`, `a3`, `a4`, `a6`
- **Tunable sections**: `w`, `a10`, `a21`, `a32`, `a43`
- **Model sections**: `w`, `a10`, `a21`, `a32`, `a43`, `a1`, `a2`, `a3`, `a4`, `a6`
- **Defining classes**: `W`, the divisor class of `w`

The next section introduces accessor functions that allow users to retrieve defining classes, model sections, tunable sections, and
related properties. Several examples are also included to demonstrate how to construct literature models by specifying defining classes.

### Accessor Functions

Retrieve the defining classes:

```@docs
defining_classes(::AbstractFTheoryModel)
```

Access all model and tunable sections:

```@docs
model_sections(::AbstractFTheoryModel)
tunable_sections(::AbstractFTheoryModel)
```

The tunable sections parametrize structural sections. The following method expresses each model section in terms of the tunable sections:

```@docs
model_section_parametrization(::AbstractFTheoryModel)
```

The hypersurface equation can also be written explicitly in terms of structural and tunable sections. This is accessible via:

!!! warning
    Currently supported only for [Hypersurface Models](@ref hypersurface_models).

```@docs
hypersurface_equation_parametrization(::HypersurfaceModel)
```

To study the divisor classes of the tunable sections in more detail, the following method provides their expression in a basis consisting
of the anticanonical divisor ``\overline{K}_B`` of the base and the model's defining classes. Each **column** of the returned matrix
corresponds to a tunable section, while each **row** represents a basis element:

```@docs
classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes(::AbstractFTheoryModel)
```

Once the base of the elliptic fibration is a concrete toric variety, we can obtain further explicit information about the model sections.

Retrieve the toric divisor class of each model section:

```@docs
classes_of_model_sections(::AbstractFTheoryModel)
```

Express model sections explicitly as polynomials in the Cox ring of the base space:

!!! warning
    If the base is a family of spaces, the polynomial expressions are equivalent to the return values of
    [`model_section_parametrization(::AbstractFTheoryModel)`](@ref).

```@docs
explicit_model_sections(::AbstractFTheoryModel)
```

---

## [Metadata Attributes](@id meta_data_attributes)

We provide the following metadata attributes for literature models, each of which returns
the corresponding attribute if it exists and otherwise throws an error:

```@docs
associated_literature_models(m::AbstractFTheoryModel)
arxiv_doi(m::AbstractFTheoryModel)
arxiv_id(m::AbstractFTheoryModel)
arxiv_link(m::AbstractFTheoryModel)
arxiv_model_equation_number(m::AbstractFTheoryModel)
arxiv_model_page(m::AbstractFTheoryModel)
arxiv_model_section(m::AbstractFTheoryModel)
arxiv_version(m::AbstractFTheoryModel)
birational_literature_models(m::AbstractFTheoryModel)
journal_doi(m::AbstractFTheoryModel)
journal_link(m::AbstractFTheoryModel)
journal_model_equation_number(m::AbstractFTheoryModel)
journal_model_page(m::AbstractFTheoryModel)
journal_model_section(m::AbstractFTheoryModel)
journal_pages(m::AbstractFTheoryModel)
journal_report_numbers(m::AbstractFTheoryModel)
journal_volume(m::AbstractFTheoryModel)
journal_name(m::AbstractFTheoryModel)
journal_year(m::AbstractFTheoryModel)
literature_identifier(m::AbstractFTheoryModel)
model_description(m::AbstractFTheoryModel)
model_parameters(m::AbstractFTheoryModel)
paper_authors(m::AbstractFTheoryModel)
paper_buzzwords(m::AbstractFTheoryModel)
paper_description(m::AbstractFTheoryModel)
paper_title(m::AbstractFTheoryModel)
```

For metadata fields that are collections, we provide helper functions to add new entries:

```@docs
add_associated_literature_model!(m::AbstractFTheoryModel, addition::String)
add_birational_literature_model!(m::AbstractFTheoryModel, addition::String)
add_journal_report_number!(m::AbstractFTheoryModel, addition::String)
add_paper_author!(m::AbstractFTheoryModel, addition::String)
add_paper_buzzword!(m::AbstractFTheoryModel, addition::String)
```

We do **not** provide `set_*` functions to overwrite attributes, in order to reduce the risk of
accidental data loss. If you truly need to replace a value, use `set_attribute!(m, :attribute_name, value)`.
