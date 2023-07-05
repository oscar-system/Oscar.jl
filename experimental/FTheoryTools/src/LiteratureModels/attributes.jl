#####################################################
# 1: Known resolutions
#####################################################

@doc raw"""
    resolutions(t::AbstractFTheoryModel)

Return the list of all known resolutions for the given model.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> length(resolutions(t))
1
```
"""
function resolutions(t::AbstractFTheoryModel)
  @req has_attribute(t, :resolutions) "No resolutions known for this model"
  return get_attribute(t, :resolutions)
end


#####################################################
# 2: Attributes for literature models
#####################################################

@doc raw"""
    doi(t::AbstractFTheoryModel)

Return the `doi` of the publication that introduced
the given model. If no `doi` is
known, an error is raised.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> doi(t)
"10.1016/j.nuclphysb.2011.12.013"
```
"""
function doi(t::AbstractFTheoryModel)
  @req has_attribute(t, :doi) "No doi known for this model"
  return get_attribute(t, :doi)
end

@doc raw"""
    arxiv_id(t::AbstractFTheoryModel)

Return the `arxiv_id` of the preprint that introduced
the given model. If no `arxiv_id` is
known, an error is raised.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> arxiv_id(t)
"1109.3454"
```
"""
function arxiv_id(t::AbstractFTheoryModel)
  @req has_attribute(t, :arxiv_id) "No arxiv identifier known for this model"
  return get_attribute(t, :arxiv_id)
end

@doc raw"""
    version(t::AbstractFTheoryModel)

Return the `version` of the arXiv preprint that
introduced the given model. If no
`version` is known, an error is raised.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> version(t)
"2"
```
"""
function version(t::AbstractFTheoryModel)
  @req has_attribute(t, :version) "No version known for this model"
  return get_attribute(t, :version)
end

@doc raw"""
    equation_number(t::AbstractFTheoryModel)

Return the `equation_number` in which the given model was introduced
in the arXiv preprint in our record. If no `equation_number`
is known, an error is raised.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> equation_number(t)
"3.1"
```
"""
function equation_number(t::AbstractFTheoryModel)
  @req has_attribute(t, :equation_number) "No equation number known for this model"
  return get_attribute(t, :equation_number)
end

@doc raw"""
    description(t::AbstractFTheoryModel)

Return the `description` of the given model.
If no `description` is known, an error is raised.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> description(t)
"SU(5)xU(1) restricted Tate model"
```
"""
function description(t::AbstractFTheoryModel)
  @req has_attribute(t, :description) "No description known for this model"
  return get_attribute(t, :description)
end

@doc raw"""
    link(t::AbstractFTheoryModel)

Return the `link` (formatted as string) to the online
version of the paper that introduced the given model.
If no `link` is known, an error is raised.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> link(t)
"https://arxiv.org/abs/1109.3454v2"
```
"""
function link(t::AbstractFTheoryModel)
  @req has_attribute(t, :link) "No link known for this model"
  return get_attribute(t, :link)
end
