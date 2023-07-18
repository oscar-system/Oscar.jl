```@meta
CurrentModule = Oscar
```

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
literature_model(; doi::String="", arxiv_id::String="", version::String="", equation::String="")
```

## Attributes

For literature models, we provide the following attributes:
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
journal_pages(m::AbstractFTheoryModel)
journal_report_numbers(m::AbstractFTheoryModel)
journal_volume(m::AbstractFTheoryModel)
journal_year(m::AbstractFTheoryModel)
literature_identifier(m::AbstractFTheoryModel)
model_description(m::AbstractFTheoryModel)
model_parameters(m::AbstractFTheoryModel)
paper_authors(m::AbstractFTheoryModel)
paper_buzzwords(m::AbstractFTheoryModel)
paper_description(m::AbstractFTheoryModel)
paper_title(m::AbstractFTheoryModel)
related_literature_models(m::AbstractFTheoryModel)
resolutions(m::AbstractFTheoryModel)
resolution_generating_sections(m::AbstractFTheoryModel)
resolution_zero_sections(m::AbstractFTheoryModel)
weighted_resolutions(m::AbstractFTheoryModel)
weighted_resolution_generating_sections(m::AbstractFTheoryModel)
weighted_resolution_zero_sections(m::AbstractFTheoryModel)
```
One can add this information for a model that does
not have it:
```@docs
set_description(m::AbstractFTheoryModel, description::String)
```
Note however, that these changes will (currently) not be stored
in our data base. One can also check if a model has a particular
set of information. This is achieved with the following methods:
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
* `has_related_literature_models(m::AbstractFTheoryModel)`,
* `has_resolutions(m::AbstractFTheoryModel)`,
* `has_resolution_generating_sections(m::AbstractFTheoryModel)`,
* `has_resolution_zero_sections(m::AbstractFTheoryModel)`,
* `has_weighted_resolutions(m::AbstractFTheoryModel)`,
* `has_weighted_resolution_generating_sections(m::AbstractFTheoryModel)`,
* `has_weighted_resolution_zero_sections(m::AbstractFTheoryModel)`,
* `has_zero_section(m::AbstractFTheoryModel)`.


## Methods

### Resolution(s) of a singular model

A central task in F-theory is to resolve a singular model.
For literature models, we have stored resolutions in our data base.
Upon construction of a literature model, we load these known resolutions.

In addition to listing the known resolutions with `resolutions(m::AbstractFTheoryModel)`, the user might want to add a resolution. This can be achieved with
the following method:
```@docs
add_resolution(m::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})
```
Provided that a resolution for a model is known, we can (attempt to) resolve the model.
```@docs
resolve(m::AbstractFTheoryModel, index::Int)
```
