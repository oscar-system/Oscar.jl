#####################################################
# 1: Attributes for literature models
#####################################################

@doc raw"""
    arxiv_id(m::AbstractFTheoryModel)

Return the `arxiv_id` of the preprint that introduced
the given model. If no `arxiv_id` is
known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> arxiv_id(m)
"1109.3454"
```
"""
function arxiv_id(m::AbstractFTheoryModel)
  @req has_arxiv_id(m) "No arxiv identifier known for this model"
  return get_attribute(m, :arxiv_id)
end

@doc raw"""
    arxiv_doi(m::AbstractFTheoryModel)

Return the `arxiv_doi` of the preprint that introduced
the given model. If no `arxiv_doi` is
known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> arxiv_doi(m)
"10.48550/arXiv.1109.3454"
```
"""
function arxiv_doi(m::AbstractFTheoryModel)
  @req has_arxiv_doi(m) "No arXiv DOI known for this model"
  return get_attribute(m, :arxiv_doi)
end

@doc raw"""
    arxiv_link(m::AbstractFTheoryModel)

Return the `arxiv_link` (formatted as string) to the arXiv
version of the paper that introduced the given model.
If no `arxiv_link` is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> arxiv_link(m)
"https://arxiv.org/abs/1109.3454v2"
```
"""
function arxiv_link(m::AbstractFTheoryModel)
  @req has_arxiv_link(m) "No arXiv link known for this model"
  return get_attribute(m, :arxiv_link)
end

@doc raw"""
    arxiv_model_equation_number(m::AbstractFTheoryModel)

Return the `arxiv_model_equation_number` in which the given model was introduced
in the arXiv preprint in our record. If no `arxiv_model_equation_number`
is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> arxiv_model_equation_number(m)
"3.1"
```
"""
function arxiv_model_equation_number(m::AbstractFTheoryModel)
  @req has_arxiv_model_equation_number(m) "No arXiv equation number known for this model"
  return get_attribute(m, :arxiv_model_equation_number)
end

@doc raw"""
    arxiv_model_page(m::AbstractFTheoryModel)

Return the `arxiv_model_page` on which the given model was introduced
in the arXiv preprint in our record. If no `arxiv_model_page`
is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> arxiv_model_page(m)
"10"
```
"""
function arxiv_model_page(m::AbstractFTheoryModel)
  @req has_arxiv_model_page(m) "No arXiv page number known for this model"
  return get_attribute(m, :arxiv_model_page)
end

@doc raw"""
    arxiv_model_section(m::AbstractFTheoryModel)

Return the `arxiv_model_section` in which the given model was introduced
in the arXiv preprint in our record. If no `arxiv_model_section`
is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> arxiv_model_section(m)
"3"
```
"""
function arxiv_model_section(m::AbstractFTheoryModel)
  @req has_arxiv_model_section(m) "No arXiv section number known for this model"
  return get_attribute(m, :arxiv_model_section)
end

@doc raw"""
    arxiv_version(m::AbstractFTheoryModel)

Return the `arxiv_version` of the arXiv preprint that
introduced the given model. If no
`arxiv_version` is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> arxiv_version(m)
"2"
```
"""
function arxiv_version(m::AbstractFTheoryModel)
  @req has_arxiv_version(m) "No arXiv version known for this model"
  return get_attribute(m, :arxiv_version)
end

@doc raw"""
    associated_literature_models(m::AbstractFTheoryModel)

Return a list of the unique identifiers any `associated_literature_models` of
the given model. These are either other presentations (Weierstrass, Tate, ...)
of the given model, or other version of the same model from a different paper
in the literature. If no `associated_literature_models` are known,
an error is raised.
"""
function associated_literature_models(m::AbstractFTheoryModel)
  @req has_associated_literature_models(m) "No associated models known for this model"
  return get_attribute(m, :associated_literature_models)
end

@doc raw"""
    generating_sections(m::AbstractFTheoryModel)

Return a list of the known Mordell–Weil generating sections of
the given model.  If no generating sections are known,
an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> generating_sections(m)
1-element Vector{Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 [0, 0, 1]
```
"""
function generating_sections(m::AbstractFTheoryModel)
  @req has_generating_sections(m) "No generating sections known for this model"
  return get_attribute(m, :generating_sections)
end

@doc raw"""
    journal_doi(m::AbstractFTheoryModel)

Return the `journal_doi` of the publication that introduced
the given model. If no `journal_doi` is
known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> journal_doi(m)
"10.1016/j.nuclphysb.2011.12.013"
```
"""
function journal_doi(m::AbstractFTheoryModel)
  @req has_journal_doi(m) "No journal DOI known for this model"
  return get_attribute(m, :journal_doi)
end

@doc raw"""
    journal_link(m::AbstractFTheoryModel)

Return the `journal_link` (formatted as string) to the published
version of the paper that introduced the given model.
If no `journal_link` is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> journal_link(m)
"https://www.sciencedirect.com/science/article/pii/S0550321311007115"
```
"""
function journal_link(m::AbstractFTheoryModel)
  @req has_journal_link(m) "No journal link known for this model"
  return get_attribute(m, :journal_link)
end

@doc raw"""
    journal_model_equation_number(m::AbstractFTheoryModel)

Return the `journal_model_equation_number` in which the given model was introduced
in the published paper in our record. If no `journal_model_equation_number`
is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> journal_model_equation_number(m)
"3.1"
```
"""
function journal_model_equation_number(m::AbstractFTheoryModel)
  @req has_journal_model_equation_number(m) "No journal equation number known for this model"
  return get_attribute(m, :journal_model_equation_number)
end

@doc raw"""
    journal_model_page(m::AbstractFTheoryModel)

Return the `journal_model_page` on which the given model was introduced
in the published paper in our record. If no `journal_model_page`
is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> journal_model_page(m)
"9"
```
"""
function journal_model_page(m::AbstractFTheoryModel)
  @req has_journal_model_page(m) "No journal page number known for this model"
  return get_attribute(m, :journal_model_page)
end

@doc raw"""
    journal_model_section(m::AbstractFTheoryModel)

Return the `journal_model_section` in which the given model was introduced
in the published paper in our record. If no `journal_model_section`
is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> journal_model_section(m)
"3"
```
"""
function journal_model_section(m::AbstractFTheoryModel)
  @req has_journal_model_section(m) "No journal section number known for this model"
  return get_attribute(m, :journal_model_section)
end

@doc raw"""
    journal_pages(m::AbstractFTheoryModel)

Return the `journal_pages` of the published paper in which the given model was introduced.
If no `journal_pages` are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> journal_pages(m)
"1–47"
```
"""
function journal_pages(m::AbstractFTheoryModel)
  @req has_journal_pages(m) "No journal pages known for this model"
  return get_attribute(m, :journal_pages)
end

@doc raw"""
    journal_report_numbers(m::AbstractFTheoryModel)

Return the `journal_report_numbers` of the published paper in which the given model was introduced.
If no `journal_report_numbers` is known, an error is raised.
"""
function journal_report_numbers(m::AbstractFTheoryModel)
  @req has_journal_report_numbers(m) "No journal pages known for this model"
  return get_attribute(m, :journal_report_numbers)
end

@doc raw"""
    journal_volume(m::AbstractFTheoryModel)

Return the `journal_volume` of the published paper in which the given model was introduced.
If no `journal_volume` are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> journal_volume(m)
"858"
```
"""
function journal_volume(m::AbstractFTheoryModel)
  @req has_journal_volume(m) "No journal volume known for this model"
  return get_attribute(m, :journal_volume)
end

@doc raw"""
    journal_year(m::AbstractFTheoryModel)

Return the `journal_year` of the published paper in which the given model was introduced.
If no `journal_year` is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> journal_year(m)
"2012"
```
"""
function journal_year(m::AbstractFTheoryModel)
  @req has_journal_year(m) "No journal year known for this model"
  return get_attribute(m, :journal_year)
end

@doc raw"""
    literature_identifier(m::AbstractFTheoryModel)

Return the `literature_identifier` of the given mode, which is a unique string
that distinguishes the model from all others in the literature model database.
If no `literature_identifier` is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> literature_identifier(m)
"1109_3454"
```
"""
function literature_identifier(m::AbstractFTheoryModel)
  @req has_literature_identifier(m) "No literature identifier known for this model"
  return get_attribute(m, :literature_identifier)
end

@doc raw"""
    model_description(m::AbstractFTheoryModel)

Return the `model_description` of the given model.
If no `model_description` is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> model_description(m)
"SU(5)xU(1) restricted Tate model"
```
"""
function model_description(m::AbstractFTheoryModel)
  @req has_model_description(m) "No model description known for this model"
  return get_attribute(m, :model_description)
end

@doc raw"""
    model_parameters(m::AbstractFTheoryModel)

Return the `model_parameters` of the given model.
If no `model_parameters` are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1212.2949", equation = "3.2", model_parameters = Dict("k" => 5))
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(2k+1) Tate model with parameter values (k = 5) based on arXiv paper 1212.2949 Eq. (3.2)

julia> model_parameters(m)
Dict{String, Int64} with 1 entry:
  "k" => 5
```
"""
function model_parameters(m::AbstractFTheoryModel)
  @req has_model_parameters(m) "No model parameters known for this model"
  return get_attribute(m, :model_parameters)
end

@doc raw"""
    paper_authors(m::AbstractFTheoryModel)

Return the `paper_authors` of the paper that introduced the given model.
If no `paper_authors` are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> paper_authors(m)
3-element Vector{String}:
 "Sven Krause"
 "Christoph Mayrhofer"
 "Timo Weigand"
```
"""
function paper_authors(m::AbstractFTheoryModel)
  @req has_paper_authors(m) "No paper authors known for this model"
  return get_attribute(m, :paper_authors)
end

@doc raw"""
    paper_buzzwords(m::AbstractFTheoryModel)

Return the `paper_buzzwords` of the paper that introduced the given model.
If no `paper_buzzwords` are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> paper_buzzwords(m)
4-element Vector{String}:
 "GUT model"
 "Tate"
 "U(1)"
 "SU(5)"
```
"""
function paper_buzzwords(m::AbstractFTheoryModel)
  @req has_paper_buzzwords(m) "No paper buzzwords known for this model"
  return get_attribute(m, :paper_buzzwords)
end

@doc raw"""
    paper_description(m::AbstractFTheoryModel)

Return the `paper_description` of the paper that introduced the given model.
If no `paper_description` is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> paper_description(m)
"SU(5)xU(1) restricted Tate model"
```
"""
function paper_description(m::AbstractFTheoryModel)
  @req has_paper_description(m) "No paper description known for this model"
  return get_attribute(m, :paper_description)
end

@doc raw"""
    paper_title(m::AbstractFTheoryModel)

Return the `paper_title` of the arXiv preprint that introduced the given model.
If no `paper_title` is known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> paper_title(m)
"\$G_4\$ flux, chiral matter and singularity resolution in F-theory compactifications"
```
"""
function paper_title(m::AbstractFTheoryModel)
  @req has_paper_title(m) "No paper title known for this model"
  return get_attribute(m, :paper_title)
end

@doc raw"""
    related_literature_models(m::AbstractFTheoryModel)

Return a list of the unique identifiers of any `related_literature_models` of
the given model. These are models that are introduced in the same paper as
the given model, but that are distinct from the given model. If no
`related_literature_models` are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1212.2949", equation = "3.2", model_parameters = Dict("k" => 5))
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(2k+1) Tate model with parameter values (k = 5) based on arXiv paper 1212.2949 Eq. (3.2)

julia> related_literature_models(m)
6-element Vector{String}:
 "1212_2949-2"
 "1212_2949-3"
 "1212_2949-4"
 "1212_2949-5"
 "1212_2949-6"
 "1212_2949-7"
```
"""
function related_literature_models(m::AbstractFTheoryModel)
  @req has_related_literature_models(m) "No associated models known for this model"
  return get_attribute(m, :related_literature_models)
end

@doc raw"""
    resolutions(m::AbstractFTheoryModel)

Return the list of all known resolutions for the given model.
If no resolutions are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> resolutions(m)
1-element Vector{Vector{Vector}}:
 [[["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]], ["e1", "e4", "e2", "e3", "s"]]
```
"""
function resolutions(m::AbstractFTheoryModel)
  @req has_resolutions(m) "No resolutions known for this model"
  return get_attribute(m, :resolutions)
end

@doc raw"""
    resolution_generating_sections(m::AbstractFTheoryModel)

Return a list of lists of known Mordell–Weil generating sections
for the given model after each known resolution. Each element of
the outer list corresponds to a known resolution (in the same order),
and each element of the list associated to a given resolution
corresponds to a known generating section (in the same order).
If no resolution generating sections are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> resolution_generating_sections(m)
1-element Vector{Vector{Vector{Vector{String}}}}:
 [[["0", "0", "1"], ["0", "0", "1"], ["0", "1"], ["0", "1"], ["0", "1"], ["a32", "-a43"]]]
```
"""
function resolution_generating_sections(m::AbstractFTheoryModel)
  @req has_resolution_generating_sections(m) "No resolution generating sections known for this model"
  return get_attribute(m, :resolution_generating_sections)
end

@doc raw"""
    resolution_zero_sections(m::AbstractFTheoryModel)

Return a list of known Mordell–Weil zero sections for the
given model after each known resolution. Each element of the
list corresponds to a known resolution (in the same order).
If no resolution zero sections are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> resolution_zero_sections(m)
1-element Vector{Vector{Vector{String}}}:
 [["1", "1", "0"], ["1", "1", "w"], ["1", "1"], ["1", "1"], ["1", "1"], ["1", "1"]]
```
"""
function resolution_zero_sections(m::AbstractFTheoryModel)
  @req has_resolution_zero_sections(m) "No resolution zero sections known for this model"
  return get_attribute(m, :resolution_zero_sections)
end

@doc raw"""
    weighted_resolutions(m::AbstractFTheoryModel)

Return the list of all known weighted resolutions for the given model.
If no weighted resolutions are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> weighted_resolutions(m)
1-element Vector{Vector{Vector}}:
 [Vector{Vector{Any}}[[["x", "y", "w"], [1, 1, 1]], [["x", "y", "w"], [1, 2, 1]], [["x", "y", "w"], [2, 2, 1]], [["x", "y", "w"], [2, 3, 1]], [["x", "y"], [1, 1]]], ["e1", "e4", "e2", "e3", "s"]]
```
"""
function weighted_resolutions(m::AbstractFTheoryModel)
  @req has_weighted_resolutions(m) "No weighted resolutions known for this model"
  return get_attribute(m, :weighted_resolutions)
end

@doc raw"""
    weighted_resolution_generating_sections(m::AbstractFTheoryModel)

Return a list of lists of known Mordell–Weil generating sections
for the given model after each known weighted resolution. Each element of
the outer list corresponds to a known weighted resolution (in the same order),
and each element of the list associated to a given weighted resolution
corresponds to a known generating section (in the same order).
If no weighted resolution generating sections are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> weighted_resolution_generating_sections(m)
1-element Vector{Vector{Vector{Vector{String}}}}:
 [[["0", "0", "1"], ["0", "0", "1"], ["0", "0", "1"], ["0", "0", "1"], ["0", "0", "1"], ["a32", "-a43"]]]
```
"""
function weighted_resolution_generating_sections(m::AbstractFTheoryModel)
  @req has_weighted_resolution_generating_sections(m) "No weighted resolution generating sections known for this model"
  return get_attribute(m, :weighted_resolution_generating_sections)
end

@doc raw"""
    weighted_resolution_zero_sections(m::AbstractFTheoryModel)

Return a list of known Mordell–Weil zero sections for the
given model after each known weighted resolution. Each element of the
list corresponds to a known weighted resolution (in the same order).
If no weighted resolution zero sections are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> weighted_resolution_zero_sections(m)
1-element Vector{Vector{Vector{String}}}:
 [["1", "1", "0"], ["1", "1", "w"], ["1", "1", "w"], ["1", "1", "w"], ["1", "1", "w"], ["1", "1"]]
```
"""
function weighted_resolution_zero_sections(m::AbstractFTheoryModel)
  @req has_weighted_resolution_zero_sections(m) "No weighted resolution zero sections known for this model"
  return get_attribute(m, :weighted_resolution_zero_sections)
end

# This example cannot be used until we support literature hypersurface models
# At that point, we should simply be able to uncomment this block

# @doc raw"""
#     zero_section(m::AbstractFTheoryModel)

# Return the zero section of the given model.
# If no zero section is known, an error is raised.
# This information is not typically stored as an attribute for
# Weierstrass and global Tate models, whose zero sections are known.

# ```jldoctest
# julia> m = literature_model(arxiv_id = "1507.05954", equation = "3.4")
# Hyoersurface model over a not fully specified base -- U(1)xU(1) hypersurface model based on arXiv paper 1507.05954 Eq. (A.1)

# julia> zero_section(m)
# 1-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
#  0
#  -b1
#  a1
# ```
# """
# function zero_section(m::AbstractFTheoryModel)
#   @req has_zero_section(m) "No zero section stored for this model"
#   return get_attribute(m, :zero_section)
# end
