##########################################
### (1) General attributes
##########################################

@doc raw"""
    base_space(m::AbstractFTheoryModel)

Return the base space of the F-theory model.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> base_space(m)
A family of spaces of dimension d = 3
```
"""
function base_space(m::AbstractFTheoryModel)
  is_base_space_fully_specified(m) || @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning AUXILIARY base space.\n"
  return m.base_space
end


@doc raw"""
    ambient_space(m::AbstractFTheoryModel)

Return the ambient space of the F-theory model.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> ambient_space(m)
A family of spaces of dimension d = 5
```
"""
function ambient_space(m::AbstractFTheoryModel)
  is_base_space_fully_specified(m) || @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning AUXILIARY ambient space.\n"
  return m.ambient_space
end


@doc raw"""
    fiber_ambient_space(m::AbstractFTheoryModel)

Return the fiber ambient space of an F-theory model.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base

julia> fiber_ambient_space(t)
Normal toric variety
```
"""
function fiber_ambient_space(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :fiber_ambient_space) "fiber_ambient_space not supported for this F-theory model"
  return m.fiber_ambient_space
end


@doc raw"""
    explicit_model_sections(m::AbstractFTheoryModel)

Return the model sections in explicit form, that
is as polynomials of the base space coordinates.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> explicit_model_sections(t)
Dict{String, QQMPolyRingElem} with 9 entries:
  "a6"  => 0
  "a21" => a21
  "a3"  => w^2*a32
  "w"   => w
  "a2"  => w*a21
  "a1"  => a1
  "a4"  => w^3*a43
  "a43" => a43
  "a32" => a32
```
"""
function explicit_model_sections(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :explicit_model_sections) "explicit_model_sections not supported for this F-theory model"
  return m.explicit_model_sections
end


@doc raw"""
    defining_section_parametrization(m::AbstractFTheoryModel)

Return the model sections in explicit form, that
is as polynomials of the base space coordinates.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> defining_section_parametrization(t)
Dict{String, MPolyRingElem} with 4 entries:
  "a6" => 0
  "a3" => w^2*a32
  "a2" => w*a21
  "a4" => w^3*a43
```
"""
function defining_section_parametrization(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :defining_section_parametrization) "defining_section_parametrization not supported for this F-theory model"
  return m.defining_section_parametrization
end


@doc raw"""
    classes_of_model_sections(m::AbstractFTheoryModel)

Return the divisor classes of all model sections.

```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> classes_of_model_sections(t)
Dict{String, ToricDivisorClass} with 9 entries:
  "a21" => Divisor class on a normal toric variety
  "a6"  => Divisor class on a normal toric variety
  "a3"  => Divisor class on a normal toric variety
  "w"   => Divisor class on a normal toric variety
  "a2"  => Divisor class on a normal toric variety
  "a1"  => Divisor class on a normal toric variety
  "a43" => Divisor class on a normal toric variety
  "a4"  => Divisor class on a normal toric variety
  "a32" => Divisor class on a normal toric variety
```
"""
@attr Dict{String, ToricDivisorClass} function classes_of_model_sections(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :explicit_model_sections) "explicit_model_sections not supported for this F-theory model"
  @req base_space(m) isa NormalToricVariety "classes_of_model_sections only supported for models over a toric base"
  my_dict = Dict{String, ToricDivisorClass}()
  for (key, value) in explicit_model_sections(m)
    sec = value
    @req is_homogeneous(sec) "Encountered a non-homogeneous model section"
    if is_zero(sec)
      my_dict[key] = trivial_divisor_class(base_space(m))
    else
      deg = collect(keys(homogeneous_components(sec)))[1]
      my_dict[key] = toric_divisor_class(base_space(m), deg)
    end
  end
  return my_dict
end


@doc raw"""
    defining_classes(m::AbstractFTheoryModel)

Return the defining divisor classes of the model in question.

```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> defining_classes(t)
Dict{String, ToricDivisorClass} with 2 entries:
  "w"    => Divisor class on a normal toric variety
  "Kbar" => Divisor class on a normal toric variety
```
"""
function defining_classes(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :defining_classes) "defining_classes not supported for this F-theory model"
  return m.defining_classes
end


##########################################
### (2) Meta data attributes
##########################################

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
    birational_literature_models(m::AbstractFTheoryModel)

Return a list of the unique identifiers of `birational_literature_models` of
the given model. These are either other presentations (Weierstrass, Tate, ...)
of the given model, or other version of the same model from a different paper
in the literature. If no `birational_literature_models` are known,
an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1507.05954", equation = "A.1")
Assuming that the first row of the given grading is the grading under Kbar

Weierstrass model over a not fully specified base -- U(1)xU(1) Weierstrass model based on arXiv paper 1507.05954 Eq. (A.1)

julia> birational_literature_models(m)
1-element Vector{String}:
 "1507_05954-1"
```
"""
function birational_literature_models(m::AbstractFTheoryModel)
  @req has_birational_literature_models(m) "No birationally equivalent models known for this model"
  return get_attribute(m, :birational_literature_models)
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

```jldoctest
julia> m = literature_model(arxiv_id = "1507.05954", equation = "A.1")
Assuming that the first row of the given grading is the grading under Kbar

Weierstrass model over a not fully specified base -- U(1)xU(1) Weierstrass model based on arXiv paper 1507.05954 Eq. (A.1)

julia> journal_report_numbers(m)
3-element Vector{String}:
 "UPR-1274-T"
 "CERN-PH-TH-2015-157"
 "MIT-CTP-4678"
```
"""
function journal_report_numbers(m::AbstractFTheoryModel)
  @req has_journal_report_numbers(m) "No journal report numbers known for this model"
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
    journal_name(m::AbstractFTheoryModel)

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
function journal_name(m::AbstractFTheoryModel)
@req has_journal_name(m) "No journal volume known for this model"
  return get_attribute(m, :journal_name)
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

Global Tate model over a not fully specified base -- SU(11) Tate model with parameter values (k = 5) based on arXiv paper 1212.2949 Eq. (3.2)

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
    associated_literature_models(m::AbstractFTheoryModel)

Return a list of the unique identifiers of any `associated_literature_models` of
the given model. These are models that are introduced in the same paper as
the given model, but that are distinct from the given model. If no
`associated_literature_models` are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1212.2949", equation = "3.2", model_parameters = Dict("k" => 5))
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(11) Tate model with parameter values (k = 5) based on arXiv paper 1212.2949 Eq. (3.2)

julia> associated_literature_models(m)
6-element Vector{String}:
 "1212_2949-2"
 "1212_2949-3"
 "1212_2949-4"
 "1212_2949-5"
 "1212_2949-6"
 "1212_2949-7"
```
"""
function associated_literature_models(m::AbstractFTheoryModel)
  @req has_associated_literature_models(m) "No associated models known for this model"
  return get_attribute(m, :associated_literature_models)
end


@doc raw"""
    model_index(m::AbstractFTheoryModel)
Return database index of a literature model. This index is a unique identifier that can be used to more conveniently construct the model. 
All models have a model_index and these will not change in the future.

```jldoctest
julia> t = literature_model(31)
Assuming that the first row of the given grading is the grading under Kbar

Weierstrass model over a not fully specified base -- F-theory weierstrass model dual to hypersurface model with fiber ambient space F_10 based on arXiv paper 1408.4808 Eq. (3.130)

julia> model_index(t)
31
```
"""
function model_index(m::AbstractFTheoryModel)
  directory = joinpath(dirname(@__DIR__), "LiteratureModels/")
  model_indices = JSON.parsefile(directory * "model_indices.json")
  return parse(Int, model_indices["model" * literature_identifier(m) * ".json"])
end



##########################################
### (3) Specialized model attributes
##########################################

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
1-element Vector{Vector{QQMPolyRingElem}}:
 [0, 0, 1]
```
"""
function generating_sections(m::AbstractFTheoryModel)
  @req has_generating_sections(m) "No generating sections known for this model"
  return get_attribute(m, :generating_sections)
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
1-element Vector{Vector{Vector{Vector{QQMPolyRingElem}}}}:
 [[[0, 0, 1], [0, 0, 1], [0, 1], [0, 1], [0, 1], [a32, -a43]]]
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
1-element Vector{Vector{Vector{QQMPolyRingElem}}}:
 [[1, 1, 0], [1, 1, w], [1, 1], [1, 1], [1, 1], [1, 1]]
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
1-element Vector{Vector{Vector{Vector{QQMPolyRingElem}}}}:
 [[[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [a32, -a43]]]
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
1-element Vector{Vector{Vector{QQMPolyRingElem}}}:
 [[1, 1, 0], [1, 1, w], [1, 1, w], [1, 1, w], [1, 1, w], [1, 1]]
```
"""
function weighted_resolution_zero_sections(m::AbstractFTheoryModel)
  @req has_weighted_resolution_zero_sections(m) "No weighted resolution zero sections known for this model"
  return get_attribute(m, :weighted_resolution_zero_sections)
end


@doc raw"""
    zero_section(m::AbstractFTheoryModel)

Return the zero section of the given model.
If no zero section is known, an error is raised.
This information is not typically stored as an attribute for
Weierstrass and global Tate models, whose zero sections are known.

```jldoctest
julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5")
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base

julia> zero_section(h)
3-element Vector{QQMPolyRingElem}:
 0
 1
 0
```
"""
function zero_section(m::AbstractFTheoryModel)
  @req has_zero_section(m) "No zero section stored for this model"
  return get_attribute(m, :zero_section)
end


@doc raw"""
    gauge_algebra(m::AbstractFTheoryModel)

Return the gauge algebra of the given model.
If no gauge algebra is known, an error is raised.
This information is typically available for all models, however.

```jldoctest
julia> t = literature_model(arxiv_id = "1408.4808", equation = "3.190", type = "hypersurface")
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base

julia> gauge_algebra(t)
Direct sum Lie algebra
  of dimension 13
with summands
  sl_2
  sl_2
  sl_2
  sl_2
  linear Lie algebra
over field of algebraic numbers
```
"""
function gauge_algebra(m::AbstractFTheoryModel)
  @req has_gauge_algebra(m) "No gauge algebra stored for this model"
  return get_attribute(m, :gauge_algebra)
end


@doc raw"""
    global_gauge_quotients(m::AbstractFTheoryModel)

Return list of lists of matrices, where each list of matrices corresponds to a gauge factor of the same index given by `gauge_algebra(m)`.
These matrices are elements of the center of the corresponding gauge factor and quotienting by them replicates the action of some discrete group on the center of the lie algebra.
This list combined with `gauge_algebra(m)` completely determines the gauge group of the model.
If no gauge quotients are known, an error is raised.

```jldoctest
julia> t = literature_model(arxiv_id = "1408.4808", equation = "3.190", type = "hypersurface")
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base

julia> global_gauge_quotients(t)
5-element Vector{Vector{String}}:
 ["-identity_matrix(C,2)", "-identity_matrix(C,2)"]
 ["-identity_matrix(C,2)"]
 ["-identity_matrix(C,2)"]
 ["-identity_matrix(C,2)", "-identity_matrix(C,2)"]
 ["-identity_matrix(C,1)"]
```
"""
function global_gauge_quotients(m::AbstractFTheoryModel)
  @req has_global_gauge_quotients(m) "No gauge quotients stored for this model"
  return get_attribute(m, :global_gauge_quotients)
end


@doc raw"""
    chern_class_c1(m::AbstractFTheoryModel)

If the elliptically fibered n-fold $Y_n$ underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute a cohomology class $h$ on the
toric ambient space $X_\Sigma$, such that its restriction to $Y_n$ is the first Chern class
$c_1$ of the tangent bundle of $Y_n$. If those assumptions are satisfied, this method returns
this very cohomology class $h$, otherwise it raises and error.

Note that $Y_n$ is a Calabi-Yau variety precisely if $c_1$ is trivial. Thus, the restriction
of $h$ to the hypersurface $Y_n$ must be trivial. Upon closer inspection, in the given toric
setting, this is equivalent to $h$ being trivial.

The computation of the cohomology ring verifies if the toric variety is simplicial and
complete. The check for it to be complete can be very time consuming. This can be switched
off by setting the optional argument `check` to the value `false`, as in the example below.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> h = chern_class_c1(qsm_model; check = false)
Cohomology class on a normal toric variety given by 0

julia> is_trivial(h)
true
```
"""
function chern_class_c1(m::AbstractFTheoryModel; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "First Chern class of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "First Chern class of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "First Chern class of F-theory model currently supported only for toric ambient space"

  # Check if the answer is known
  if has_attribute(m, :chern_class_c1)
    return get_attribute(m, :chern_class_c1)
  end

  # Trigger potential short-cut computation of cohomology ring
  cohomology_ring(ambient_space(m); check)

  # Compute the cohomology class corresponding to the hypersurface equation
  if m isa WeierstrassModel
    cl = toric_divisor_class(ambient_space(m), degree(weierstrass_polynomial(m)))
  end
  if m isa GlobalTateModel
    cl = toric_divisor_class(ambient_space(m), degree(tate_polynomial(m)))
  end
  if m isa HypersurfaceModel
    cl = toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))
  end
  cy = cohomology_class(cl)

  # Compute and set h
  c1_t_ambient = cohomology_class(anticanonical_divisor(ambient_space(m)))
  set_attribute!(m, :chern_class_c1, c1_t_ambient - cy)
  return get_attribute(m, :chern_class_c1)
end


@doc raw"""
    chern_class_c2(m::AbstractFTheoryModel)

If the elliptically fibered n-fold $Y_n$ underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute a cohomology class $h$ on the
toric ambient space $X_\Sigma$, such that its restriction to $Y_n$ is the 2nd Chern class
$c_2$ of the tangent bundle of $Y_n$. If those assumptions are satisfied, this method returns
this very cohomology class $h$, otherwise it raises and error.

As of right now, this method is computationally quite expensive for involved toric varieties,
such as in the example below. Therefore, think carefully if you truly want to compute this quantity.

The computation of the cohomology ring verifies if the toric variety is simplicial and
complete. The check for it to be complete can be very time consuming. This can be switched
off by setting the optional argument `check` to the value `false`, as in the example below.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> h = chern_class_c2(qsm_model; check = false);

julia> is_trivial(h)
false
```
"""
function chern_class_c2(m::AbstractFTheoryModel; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Second Chern class of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Second Chern class of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Second Chern class of F-theory model currently supported only for toric ambient space"

  # Check if the answer is known
  if has_attribute(m, :chern_class_c2)
    return get_attribute(m, :chern_class_c2)
  end

  # Trigger potential short-cut computation of cohomology ring
  cohomology_ring(ambient_space(m); check)

  # Compute the cohomology class corresponding to the hypersurface equation
  if m isa WeierstrassModel
    cl = toric_divisor_class(ambient_space(m), degree(weierstrass_polynomial(m)))
  end
  if m isa GlobalTateModel
    cl = toric_divisor_class(ambient_space(m), degree(tate_polynomial(m)))
  end
  if m isa HypersurfaceModel
    cl = toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))
  end
  cy = cohomology_class(cl)

  # Compute sum of products of cohomolgy classes of torus invariant prime divisors
  c_ds = [polynomial(cohomology_class(d)) for d in torusinvariant_prime_divisors(ambient_space(m))]
  c2_t_ambient = sum(c_ds[i]*c_ds[j] for i in 1:length(c_ds)-1 for j in i+1:length(c_ds))
  c2_t_ambient = cohomology_class(ambient_space(m), c2_t_ambient)

  # Compute and set h
  h = c2_t_ambient - cy * chern_class_c1(m; check = check)
  set_attribute!(m, :chern_class_c2, h)
  return h
end


##########################################
### (4) Attributes specially for the QSMs
##########################################

### (4.1) Attributes regarding the polytope in the Kreuzer-Skarke database

@doc raw"""
    vertices(m::AbstractFTheoryModel)

This method returns the vertices of the polytope the the base of the F-theory QSM is
build from. Note that those vertices are normalized according to the Polymake standard
to rational numbers.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> vertices(qsm_model)
4-element Vector{Vector{QQFieldElem}}:
 [-1, -1, -1]
 [1, -1//2, -1//2]
 [-1, 2, -1]
 [-1, -1, 5]
```
"""
function vertices(m::AbstractFTheoryModel)
  @req has_attribute(m, :vertices) "No vertices known for this model"
  return get_attribute(m, :vertices)
end


@doc raw"""
    polytope_index(m::AbstractFTheoryModel)

Of the 3-dimensional reflexive polytope that the base of this F-theory model is build from,
this method returns the index within the Kreuzer-Skarke list.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> polytope_index(qsm_model)
4
```
"""
function polytope_index(m::AbstractFTheoryModel)
  @req has_attribute(m, :poly_index) "No polytope index known for this model"
  return get_attribute(m, :poly_index)
end


@doc raw"""
    has_quick_triangulation(m::AbstractFTheoryModel)

For a 3-dimensional reflexive polytope in the Kreuzer-Skarke list, the list of
full (sometimes also called fine), regular, star triangulations can be extremely
large. Consequently, one may wonder if the triangulations can be enumerated in a
somewhat reasonable time (say 5 minutes on a personal computer). This method tries
to provide an answer to this. It returns `true` if one should expect a timly response
to the atttempt to enumerate all (full, regular, star) triangulations. Otherwise, this
method returns `false`.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> has_quick_triangulation(qsm_model)
true
```
"""
function has_quick_triangulation(m::AbstractFTheoryModel)
  @req has_attribute(m, :triang_quick) "It is not known if the base of this model can be triangulated quickly"
  return get_attribute(m, :triang_quick)::Bool
end


@doc raw"""
    max_lattice_pts_in_facet(m::AbstractFTheoryModel)

In order to enumerate the number of full, regular, star triangulations of a
3-dimensional reflexive polytope, it is possible to first find the corresponding
triangulations of all facets of the polytope [HT17](@cite). A first indication for
the complexity of this triangulation task is the maximum number of lattice points
in a facet of the polytope in question. This method returns this maximal number of
lattice points.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> max_lattice_pts_in_facet(qsm_model)
16
```
"""
function max_lattice_pts_in_facet(m::AbstractFTheoryModel)
  @req has_attribute(m, :max_lattice_pts_in_facet) "Maximal number of lattice points of facets not known for this model"
  return get_attribute(m, :max_lattice_pts_in_facet)
end


@doc raw"""
    estimated_number_of_triangulations(m::AbstractFTheoryModel)

This method returns an estimate for the number of full, regular, star triangulations
of the 3-dimensional reflexive polytope, those triangulations define the possible base
spaces of the F-theory QSM in question.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> estimated_number_of_triangulations(qsm_model)
212533333333
```
"""
function estimated_number_of_triangulations(m::AbstractFTheoryModel)
  @req has_attribute(m, :estimated_number_of_triangulations) "Estimated number of (full, regular, star) triangulation not known for this model"
  return get_attribute(m, :estimated_number_of_triangulations)
end


### (4.2) Attributes regarding the polytope in the Kreuzer-Skarke database


@doc raw"""
    kbar3(m::AbstractFTheoryModel)

Let Kbar denote the anticanonical class of the 3-dimensional base space of the F-theory QSM.
Of ample importance is the triple intersection number of Kbar, i.e. Kbar * Kbar * Kbar.
This method returns this intersection number.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> kbar3(qsm_model)
6
```
"""
function kbar3(m::AbstractFTheoryModel)
  @req has_attribute(m, :Kbar3) "Kbar3 not known for this model"
  return get_attribute(m, :Kbar3)
end


@doc raw"""
    hodge_h11(m::AbstractFTheoryModel)

This methods return the Hodge number h11 of the elliptically
fibered 4-fold that defined the F-theory QSM in question.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> hodge_h11(qsm_model)
31
```
"""
function hodge_h11(m::AbstractFTheoryModel)
  @req has_attribute(m, :h11) "Hodge number h11 of ambient space not known for this model"
  return get_attribute(m, :h11)
end


@doc raw"""
    hodge_h12(m::AbstractFTheoryModel)

This methods return the Hodge number h12 of the elliptically
fibered 4-fold that defined the F-theory QSM in question.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> hodge_h12(qsm_model)
10
```
"""
function hodge_h12(m::AbstractFTheoryModel)
  @req has_attribute(m, :h12) "Hodge number h12 of ambient space not known for this model"
  return get_attribute(m, :h12)
end


@doc raw"""
    hodge_h13(m::AbstractFTheoryModel)

This methods return the Hodge number h13 of the elliptically
fibered 4-fold that defined the F-theory QSM in question.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> hodge_h13(qsm_model)
34
```
"""
function hodge_h13(m::AbstractFTheoryModel)
  @req has_attribute(m, :h13) "Hodge number h13 of ambient space not known for this model"
  return get_attribute(m, :h13)
end


@doc raw"""
    hodge_h22(m::AbstractFTheoryModel)

This methods return the Hodge number h22 of the elliptically
fibered 4-fold that defined the F-theory QSM in question.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> hodge_h22(qsm_model)
284
```
"""
function hodge_h22(m::AbstractFTheoryModel)
  @req has_attribute(m, :h22) "Hodge number h22 of ambient space not known for this model"
  return get_attribute(m, :h22)
end


### (4.3) Attributes regarding the Ci-curves


@doc raw"""
    genera_of_ci_curves(m::AbstractFTheoryModel)

This methods return the genera of the Ci curves.
Recall that Ci = V(xi, s), where xi is a homogeneous
coordinate of the 3-dimensional toric base space B3 of the
QSM hypersurface model in question, and s is a generic
section of the anticanonical bundle of B3. Consequently,
we may use the coordinates xi as labels for the curves Ci.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> my_key = collect(keys(genera_of_ci_curves(qsm_model)))[1]
x7

julia> genera_of_ci_curves(qsm_model)[my_key]
0
```
"""
function genera_of_ci_curves(m::AbstractFTheoryModel)
  @req has_attribute(m, :genus_ci) "Genera of Ci curves not known for this model"
  return get_attribute(m, :genus_ci)
end


@doc raw"""
    degrees_of_kbar_restrictions_to_ci_curves(m::AbstractFTheoryModel)

The anticanonical divisor of the 3-dimensional toric base space B3 of the
QSM hypersurface model in question can be restricted to the Ci curves. The
result of this operation is a line bundle. This method returns the degree of
this line bundle for every Ci curve.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> my_key = collect(keys(degrees_of_kbar_restrictions_to_ci_curves(qsm_model)))[1]
x7

julia> degrees_of_kbar_restrictions_to_ci_curves(qsm_model)[my_key]
0
```
"""
function degrees_of_kbar_restrictions_to_ci_curves(m::AbstractFTheoryModel)
  @req has_attribute(m, :degree_of_Kbar_of_tv_restricted_to_ci) "Degree of Kbar restriction to Ci curves not known for this model"
  return get_attribute(m, :degree_of_Kbar_of_tv_restricted_to_ci)
end


@doc raw"""
    topological_intersection_numbers_among_ci_curves(m::AbstractFTheoryModel)

The topological intersection numbers among Ci curves are also of ample importance.
This method returns those intersection numbers in the form of a matrix.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> n_rows(topological_intersection_numbers_among_ci_curves(qsm_model))
29

julia> n_columns(topological_intersection_numbers_among_ci_curves(qsm_model))
29
```
"""
function topological_intersection_numbers_among_ci_curves(m::AbstractFTheoryModel)
  @req has_attribute(m, :intersection_number_among_ci_cj) "Topological intersection numbers among Ci curves not known for this model"
  return get_attribute(m, :intersection_number_among_ci_cj)
end


@doc raw"""
    indices_of_trivial_ci_curves(m::AbstractFTheoryModel)

Some of the Ci curves are trivial, in that V(xi, s) is the empty set.
This method returns the vector of all indices of trivial Ci curves.
That is, should V(x23, s) be the empty set, then 23 will be included in
the returned list.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> indices_of_trivial_ci_curves(qsm_model)
10-element Vector{Int64}:
 23
 22
 18
 19
 20
 26
 10
 11
 12
 15
```
"""
function indices_of_trivial_ci_curves(m::AbstractFTheoryModel)
  @req has_attribute(m, :index_facet_interior_divisors) "Degree of Kbar restriction to Ci curves not known for this model"
  return get_attribute(m, :index_facet_interior_divisors)
end


@doc raw"""
    topological_intersection_numbers_among_nontrivial_ci_curves(m::AbstractFTheoryModel)

The topological intersection numbers among the non-trivial Ci curves are used
frequently. This method returns those intersection numbers in the form of a matrix.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> n_rows(topological_intersection_numbers_among_nontrivial_ci_curves(qsm_model))
19

julia> n_columns(topological_intersection_numbers_among_nontrivial_ci_curves(qsm_model))
19
```
"""
function topological_intersection_numbers_among_nontrivial_ci_curves(m::AbstractFTheoryModel)
  @req has_attribute(m, :intersection_number_among_nontrivial_ci_cj) "Topological intersection numbers among non-trivial Ci curves not known for this model"
  return get_attribute(m, :intersection_number_among_nontrivial_ci_cj)
end


### (4.4) Attributes regarding the dual graph


@doc raw"""
    dual_graph(m::AbstractFTheoryModel)

This method returns the dual graph of the QSM model in question.
Note that no labels are (currently) attached to the vertices/nodes or edges.
To understand/read this graph correctly, please use the methods listed below.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> dual_graph(qsm_model)
Undirected graph with 21 nodes and the following edges:
(5, 1)(6, 5)(7, 6)(8, 7)(9, 4)(9, 8)(10, 1)(11, 4)(12, 3)(12, 10)(13, 3)(13, 11)(14, 1)(15, 4)(16, 3)(17, 3)(18, 2)(18, 14)(19, 2)(19, 15)(20, 2)(20, 16)(21, 2)(21, 17)
```
"""
function dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :dual_graph) "Dual graph not known for this model"
  return get_attribute(m, :dual_graph)
end


@doc raw"""
    components_of_dual_graph(m::AbstractFTheoryModel)

This method returns a vector with labels for each node/vertex of the dual graph of the QSM
model in question. Those labels allow to understand the geometric origin of the node/vertex.

Specifically, recall that those nodes are associated to the Ci-curves, which are in turn
given by Ci = V(xi, s). xi is a homogenous coordinate of the 3-dimensional toric base space
B3 of the QSM in question, and s is a generic section of the anticanonical bundle of B3.

Only non-trivial Ci = V(xi, s) correspond to vertices/nodes of the dual graph.

If Ci = V(xi, s) is irreducible and corresponds to the k-th component, then the label "Ci"
appears at position k of the vector returned by this method. However, if Ci = V(xi, s) is
reducible, then we introduce the labels Ci-0, Ci-1, Ci-2 etc. for those irreducible
components of Ci. If Ci-0 corresponds to the k-th components of the dual graph,
then the label "Ci-0" appears at position k of the vector returned by this method.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> components_of_dual_graph(qsm_model)
21-element Vector{String}:
 "C0"
 "C1"
 "C2"
 "C3"
 "C4"
 "C5"
 "C6"
 "C7"
 "C8"
 "C9"
 ⋮
 "C16"
 "C17"
 "C21"
 "C24-0"
 "C24-1"
 "C25"
 "C27"
 "C28-0"
 "C28-1"
```
"""
function components_of_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :components_of_dual_graph) "Components of dual graph not known for this model"
  return get_attribute(m, :components_of_dual_graph)
end


@doc raw"""
    degrees_of_kbar_restrictions_to_components_of_dual_graph(m::AbstractFTheoryModel)

The anticanonical bundle of the toric 3-dimensional base space of the F-theory QSM in
question can be restricted to the (geometric counterparts of the) nodes/vertices of
the dual graph. The result is a line bundle for each node/vertex. This method returns
a vector with the degrees of these line bundles.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> degrees_of_kbar_restrictions_to_components_of_dual_graph(qsm_model)["C28-1"]
0
```
"""
function degrees_of_kbar_restrictions_to_components_of_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph) "Degree of Kbar restricted to components of dual graph not known for this model"
  return get_attribute(m, :degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph)
end


@doc raw"""
    genera_of_components_of_dual_graph(m::AbstractFTheoryModel)

This methods returns a vector with the genera of the (geometric
counterparts of the) nodes/vertices of the dual graph.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> genera_of_components_of_dual_graph(qsm_model)["C28-1"]
0
```
"""
function genera_of_components_of_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :genus_of_components_of_dual_graph) "Genera of components of dual graph not known for this model"
  return get_attribute(m, :genus_of_components_of_dual_graph)
end


### (4.5) Attributes regarding the simplified dual graph


@doc raw"""
    simplified dual_graph(m::AbstractFTheoryModel)

This method returns the simplified dual graph of the QSM model in question.
Note that no labels are (currently) attached to the vertices/nodes or edges.
To understand/read this graph correctly, please use the methods listed below.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> simplified_dual_graph(qsm_model)
Undirected graph with 4 nodes and the following edges:
(2, 1)(3, 1)(3, 2)(4, 1)(4, 2)(4, 3)
```
"""
function simplified_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :simplified_dual_graph) "Simplified dual graph not known for this model"
  return get_attribute(m, :simplified_dual_graph)
end


@doc raw"""
    components_of_simplified_dual_graph(m::AbstractFTheoryModel)

This method returns a vector with labels for each node/vertex of the simplified dual graph.
Otherwise, works identical to `components_of_dual_graph(m::AbstractFTheoryModel)`.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> components_of_simplified_dual_graph(qsm_model)
4-element Vector{String}:
 "C0"
 "C1"
 "C2"
 "C3"
```
"""
function components_of_simplified_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :components_of_simplified_dual_graph) "Components of simplified dual graph not known for this model"
  return get_attribute(m, :components_of_simplified_dual_graph)
end


@doc raw"""
    degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph(m::AbstractFTheoryModel)

Same as `degrees_of_kbar_restrictions_to_components_of_dual_graph(m::AbstractFTheoryModel)`,
but for the simplified dual graph.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph(qsm_model)["C2"]
2
```
"""
function degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph) "Degree of Kbar restricted to components of simplified dual graph not known for this model"
  return get_attribute(m, :degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph)
end


@doc raw"""
    genera_of_components_of_simplified_dual_graph(m::AbstractFTheoryModel)

This methods returns a vector with the genera of the (geometric
counterparts of the) nodes/vertices of the dual graph.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> genera_of_components_of_simplified_dual_graph(qsm_model)["C2"]
0
```
"""
function genera_of_components_of_simplified_dual_graph(m::AbstractFTheoryModel)
  @req has_attribute(m, :genus_of_components_of_simplified_dual_graph) "Genera of components of simplified dual graph not known for this model"
  return get_attribute(m, :genus_of_components_of_simplified_dual_graph)
end
