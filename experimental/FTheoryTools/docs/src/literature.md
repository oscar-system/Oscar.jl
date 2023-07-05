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
doi(t::AbstractFTheoryModel)
arxiv_id(t::AbstractFTheoryModel)
version(t::AbstractFTheoryModel)
equation_number(t::AbstractFTheoryModel)
description(t::AbstractFTheoryModel)
link(t::AbstractFTheoryModel)
```
One can add this information for a model that does
not have it:
```@docs
set_description(t::AbstractFTheoryModel, description::String)
```
Note however, that these changes will (currently) not be stored
in our data base. One can also check if a model has a particular
set of information. This is achieved with the following methods:
* `has_doi(t::AbstractFTheoryModel)`,
* `has_arxiv_id(t::AbstractFTheoryModel)`,
* `has_version(t::AbstractFTheoryModel)`,
* `has_equation_number(t::AbstractFTheoryModel)`,
* `has_description(t::AbstractFTheoryModel)`,
* `has_link(t::AbstractFTheoryModel)`.


## Methods

### Resolution(s) of a singular model

A central task in F-theory is to resolve a singular model.
For literature models, we have stored resolutions in our data base.
Upon construction of a literature model, we load these known resolutions.
The user can list all the resolutions in our database as follows:
```@docs
resolutions(t::AbstractFTheoryModel)
```
In addition, the user might want to add a resolution. This can be achieved with
the following method:
```@docs
add_resolution(t::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})
```
Provided that a resolution for a model is known, we can (attempt to) resolve the model.
```@docs
resolve(t::AbstractFTheoryModel, index::Int)
```
