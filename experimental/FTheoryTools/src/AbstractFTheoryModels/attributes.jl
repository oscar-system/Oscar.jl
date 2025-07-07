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
Family of spaces of dimension d = 3
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
Family of spaces of dimension d = 5
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
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

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

Return the model sections in explicit form, that is as polynomials
of the base space coordinates. The set of keys of the returned
dictionary matches the output of `model_sections`. 

```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
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
    model_section_parametrization(m::AbstractFTheoryModel)

Return a dictionary that defines how the "default" parameters of
the given model type are defined in terms of the tunable sections
(those returned by `tunable_sections`).

```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> model_section_parametrization(t)
Dict{String, MPolyRingElem} with 4 entries:
  "a6" => 0
  "a3" => w^2*a32
  "a2" => w*a21
  "a4" => w^3*a43
```
"""
function model_section_parametrization(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :model_section_parametrization) "model_section_parametrization not supported for this F-theory model"
  return m.model_section_parametrization
end


@doc raw"""
    classes_of_model_sections(m::AbstractFTheoryModel)

Return the divisor classes of all model sections. The set
of keys of the returned dictionary matches the output
of `model_sections`.

```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
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
This is the minimum set of information required to specify a
family of models. 

```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> defining_classes(t)
Dict{String, ToricDivisorClass} with 1 entry:
  "w" => Divisor class on a normal toric variety
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
  return get_attribute(m, :arxiv_id)::String
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
  return get_attribute(m, :arxiv_doi)::String
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
  return get_attribute(m, :arxiv_link)::String
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
  return get_attribute(m, :arxiv_model_equation_number)::String
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
  return get_attribute(m, :arxiv_model_page)::String
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
  return get_attribute(m, :arxiv_model_section)::String
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
  return get_attribute(m, :arxiv_version)::String
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
  return get_attribute(m, :birational_literature_models)::Vector{String}
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
  return get_attribute(m, :journal_doi)::String
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
  return get_attribute(m, :journal_link)::String
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
  return get_attribute(m, :journal_model_equation_number)::String
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
  return get_attribute(m, :journal_model_page)::String
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
  return get_attribute(m, :journal_model_section)::String
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
  return get_attribute(m, :journal_pages)::String
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
  return get_attribute(m, :journal_report_numbers)::Vector{String}
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
  return get_attribute(m, :journal_volume)::String
end


@doc raw"""
    journal_name(m::AbstractFTheoryModel)

Return the `journal_name` of the published paper in which the given model was introduced.
If no `journal_name` are known, an error is raised.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> journal_name(m)
"Nucl. Phys. B"
```
"""
function journal_name(m::AbstractFTheoryModel)
@req has_journal_name(m) "No journal volume known for this model"
  return get_attribute(m, :journal_name)::String
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
  return get_attribute(m, :journal_year)::String
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
  return get_attribute(m, :literature_identifier)::String
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
  return get_attribute(m, :model_description)::String
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
  return get_attribute(m, :model_parameters)::Dict{String, Int64}
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
  return get_attribute(m, :paper_authors)::Vector{String}
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
  return get_attribute(m, :paper_buzzwords)::Vector{String}
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
  return get_attribute(m, :paper_description)::String
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
  return get_attribute(m, :paper_title)::String
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
  return get_attribute(m, :associated_literature_models)::Vector{String}
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
1-element Vector{Tuple{Vector{Vector{String}}, Vector{String}}}:
 ([["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]], ["e1", "e4", "e2", "e3", "s"])
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
1-element Vector{Tuple{Vector{Tuple{Vector{String}, Vector{Int64}}}, Vector{String}}}:
 ([(["x", "y", "w"], [1, 1, 1]), (["x", "y", "w"], [1, 2, 1]), (["x", "y", "w"], [2, 2, 1]), (["x", "y", "w"], [2, 3, 1]), (["x", "y"], [1, 1])], ["e1", "e4", "e2", "e3", "s"])
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
    zero_section_class(m::AbstractFTheoryModel)

Return the zero section class of a model as a cohomology class in the toric ambient space.
If no zero section class is known, an error is raised.
This information is always available for
Weierstrass and global Tate models, whose zero section classes are known.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> zero_section_class(qsm_model)
Cohomology class on a normal toric variety given by e2 + 2*u + 3*e4 + e1 - w
```
"""
function zero_section_class(m::AbstractFTheoryModel)
  @req has_zero_section_class(m) "No zero section class stored for this model"
  return get_attribute(m, :zero_section_class)
end


@doc raw"""
    zero_section_index(m::AbstractFTheoryModel)

Return the index of the generator of the Cox ring of the ambient space, whose corresponding vanishing locus defines the zero section of a model.
If no zero section class is known, an error is raised. This attribute is always set simultaneously with zero_section_class.
This information is always available for
Weierstrass and global Tate models, whose zero section classes are known.

```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> foah15_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.190", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Hypersurface model over a concrete base

julia> zero_section_index(foah15_B3)
5
```
"""
function zero_section_index(m::AbstractFTheoryModel)
  @req has_zero_section_class(m) "No zero section class stored for this model"
  return get_attribute(m, :zero_section_index)::Int
end


@doc raw"""
    exceptional_classes(m::AbstractFTheoryModel)

Return the cohomology classes of the exceptional toric divisors of a model as a vector of cohomology classes in the toric ambient space.
This information is only supported for models over a concrete base that is a normal toric variety, but is always available in this case.
After a toric blow up this information is updated.

```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> foah11_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.142", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Hypersurface model over a concrete base

julia> exceptional_classes(foah11_B3)
4-element Vector{CohomologyClass}:
 Cohomology class on a normal toric variety given by e1
 Cohomology class on a normal toric variety given by e2
 Cohomology class on a normal toric variety given by e3
 Cohomology class on a normal toric variety given by e4
```
"""
function exceptional_classes(m::AbstractFTheoryModel)
  @req base_space(m) isa NormalToricVariety "Exceptional divisor classes are only supported for models over a concrete base"
  return get_attribute(m, :exceptional_classes, Vector{CohomologyClass}())
end


@doc raw"""
    exceptional_divisor_indices(m::AbstractFTheoryModel)

Return the indices of the generators of the Cox ring of the ambient space which correspond to exceptional divisors.
This information is only supported for models over a concrete base that is a normal toric variety, but is always available in this case.
After a toric blow up this information is updated.

```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> foah11_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.142", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Hypersurface model over a concrete base

julia> exceptional_divisor_indices(foah11_B3)
4-element Vector{Int64}:
  8
  9
 10
 11
```
"""
@attr Vector{Int} function exceptional_divisor_indices(m::AbstractFTheoryModel)
  @req base_space(m) isa NormalToricVariety "Exceptional divisor indices are only supported for models over a concrete base"
  return get_attribute(m, :exceptional_divisor_indices, Vector{Int}())
end


@doc raw"""
    torsion_sections(m::AbstractFTheoryModel)

Return the torsion sections of the given model.
If no torsion sections are known, an error is raised.

```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> foah15_B3 = literature_model(arxiv_id = "1408.4808", equation = "3.190", type = "hypersurface", base_space = B3, defining_classes = Dict("s7" => Kbar, "s9" => Kbar))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Hypersurface model over a concrete base

julia> length(torsion_sections(foah15_B3))
1
```
"""
function torsion_sections(m::AbstractFTheoryModel)
  @req has_torsion_sections(m) "No torsion sections stored for this model"
  return get_attribute(m, :torsion_sections)
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
over algebraic closure of rational field
```
"""
function gauge_algebra(m::AbstractFTheoryModel)
  @req has_gauge_algebra(m) "No gauge algebra stored for this model"
  return get_attribute(m, :gauge_algebra)
end


@doc raw"""
    global_gauge_group_quotient(m::AbstractFTheoryModel)

Return list of lists of matrices, where each list of matrices corresponds to a gauge factor of the same index given by `gauge_algebra(m)`.
These matrices are elements of the center of the corresponding gauge factor and quotienting by them replicates the action of some discrete group on the center of the lie algebra.
This list combined with `gauge_algebra(m)` completely determines the gauge group of the model.
If no gauge quotients are known, an error is raised.

```jldoctest
julia> t = literature_model(arxiv_id = "1408.4808", equation = "3.190", type = "hypersurface")
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base

julia> global_gauge_group_quotient(t)
5-element Vector{Vector{String}}:
 ["-identity_matrix(C,2)", "-identity_matrix(C,2)"]
 ["-identity_matrix(C,2)"]
 ["-identity_matrix(C,2)"]
 ["-identity_matrix(C,2)", "-identity_matrix(C,2)"]
 ["-identity_matrix(C,1)"]
```
"""
function global_gauge_group_quotient(m::AbstractFTheoryModel)
  @req has_global_gauge_group_quotient(m) "No gauge quotients stored for this model"
  return get_attribute(m, :global_gauge_quotients)
end


@doc raw"""
    chern_class(m::AbstractFTheoryModel, k::Int; check::Bool = true)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute a cohomology class ``h`` on the
toric ambient space ``X_\Sigma``, such that its restriction to ``Y_n`` is the k-th Chern class
``c_k`` of the tangent bundle of ``Y_n``. If those assumptions are satisfied, this method returns
this very cohomology class ``h``, otherwise it raises an error.

The theory guarantees that the implemented algorithm works for toric ambient spaces which are
smooth and complete. The check for completeness can be very time consuming. This check can
be switched off by setting the optional argument `check` to the value `false`, as demonstrated below.

!!!warning
    This method works ONLY for F-theory models which are hypersurfaces in a toric ambient space.
  
!!!warning
    This method represents the Chern classes of said hypersurface by cohomology classes on the toric ambient space.
    These classes counterparts must be restricted to the hypersurface to truly represent the Chern class in question.
    Internally, we integrate those ambient space classes against the class of the hypersurface, which automatically
    executes the restriction to the hypersurface.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> h = chern_class(qsm_model, 4; check = false);

julia> is_trivial(h)
false
```
"""
function chern_class(m::AbstractFTheoryModel, k::Int; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Chern class of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric ambient space"

  # Consistency checks
  @req k >= 0 "Chern class index must be non-negative"
  @req k <= dim(ambient_space(m)) - 1 "Chern class index must not exceed dimension of the space"

  # If thus far, no non-trivial Chern classes have been computed for this toric variety, add an "empty" vector
  if !has_attribute(m, :chern_classes)
    cs = Dict{Int64, CohomologyClass}()
    cs[0] = cohomology_class(ambient_space(m), one(cohomology_ring(ambient_space(m), check = check)), quick = true)
    diff = degree(leading_term(hypersurface_equation(m))) - sum(cox_ring(ambient_space(m)).d)
    cs[1] = cohomology_class(toric_divisor_class(ambient_space(m), diff), quick = true)
    set_attribute!(m, :chern_classes, cs)
    if k == 0
      return cs[0]
    elseif k == 1
      return cs[1]
    end
  end

  # Check if the Chern class in question is known
  cs = get_attribute(m, :chern_classes)::Dict{Int64, CohomologyClass}
  if haskey(cs, k)
    return cs[k]
  end

  # Check if we can compute the Chern classes for the toric ambient space
  if check
    @req is_smooth(ambient_space(m)) && is_complete(ambient_space(m)) "The Chern classes of the tangent bundle of the toric ambient space are only supported if the toric ambient space is smooth and complete"
  end

  # Chern class is not known, so compute and return it...
  cy = cohomology_class(toric_divisor_class(ambient_space(m), degree(leading_term(hypersurface_equation(m)))), quick = true)
  ck_ambient = chern_class(ambient_space(m), k, check = check)
  ckm1 = chern_class(m, k-1, check = check)
  new_poly = lift(polynomial(ck_ambient)) - lift(polynomial(cy)) * lift(polynomial(ckm1))
  coho_R = cohomology_ring(ambient_space(m), check = check)
  cs[k] = cohomology_class(ambient_space(m), coho_R(new_poly), quick = true)
  set_attribute!(m, :chern_classes, cs)
  return cs[k]
end


@doc raw"""
    chern_classes(m::AbstractFTheoryModel; check::Bool = true)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute a cohomology class ``h`` on the
toric ambient space ``X_\Sigma``, such that its restriction to ``Y_n`` is the k-th Chern class
``c_k`` of the tangent bundle of ``Y_n``. If those assumptions are satisfied, this method returns
a vector with the cohomology classes corresponding to all non-trivial Chern classes ``c_k`` of
``Y_n``. Otherwise, this methods raises an error.

As of right now, this method is computationally expensive for involved toric ambient spaces,
such as in the example below.

The theory guarantees that the implemented algorithm works for toric ambient spaces which are
simplicial and complete. The check for completeness can be very time consuming. This check can
be switched off by setting the optional argument `check` to the value `false`, as demonstrated below.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> h = chern_classes(qsm_model; check = false);

julia> is_one(polynomial(h[0]))
true

julia> is_trivial(h[1])
true

julia> is_trivial(h[2])
false
```
"""
function chern_classes(m::AbstractFTheoryModel; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Chern class of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric ambient space"
  for k in 0:dim(ambient_space(m))-1
    chern_class(m, k; check = check)
  end
  return get_attribute(m, :chern_classes)::Dict{Int64, CohomologyClass}
end


@doc raw"""
    euler_characteristic(m::AbstractFTheoryModel; check::Bool = true)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute the Euler characteristic. If this
assumptions is satisfied, this method returns the Euler characteristic, otherwise it raises an
error.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> h = euler_characteristic(qsm_model; check = false)
378

julia> h = euler_characteristic(qsm_model; check = false)
378
```
"""
@attr Int function euler_characteristic(m::AbstractFTheoryModel; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Euler characteristic of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Euler characteristic of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Euler characteristic of F-theory model currently supported only for toric ambient space"

  # Trigger potential short-cut computation of cohomology ring
  cohomology_ring(ambient_space(m); check)

  # Compute the cohomology class corresponding to the hypersurface equation
  cy = cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m))))

  # Compute the Euler characteristic
  return Int(integrate(chern_class(m, 4; check) * cy; check))
end


@doc raw"""
    tunable_sections(m::AbstractFTheoryModel)

Return a vector containing all sections that can be tuned.
This is a list of the names of all parameters appearing in the model.

```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> tunable_sections(m)
5-element Vector{String}:
 "a21"
 "w"
 "a1"
 "a43"
 "a32"
```
"""
@attr Vector{String} tunable_sections(m::AbstractFTheoryModel) = collect(setdiff(keys(explicit_model_sections(m)), keys(model_section_parametrization(m))))


@doc raw"""
    model_sections(m::AbstractFTheoryModel)

Return a vector containing all sections that were used in the definition of the model.
This includes the sections returned by `tunable_sections` and all sections parametrized by
them (the keys of the dictionary returned by `model_section_parametrization`). These are
the keys of the dictionaries returned by of `explicit_model_sections` and `classes_of_model_sections`.

```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> model_sections(m)
9-element Vector{String}:
 "a6"
 "a21"
 "a3"
 "w"
 "a2"
 "a1"
 "a4"
 "a43"
 "a32"
```
"""
model_sections(m::AbstractFTheoryModel) = collect(keys(explicit_model_sections(m)))


@doc raw"""
    classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes(m::AbstractFTheoryModel)

Returns a dictionary giving the classes of all parameters (tunable sections)
in terms of Kbar and the defining classes. Each value gives the divisor
class of the corresponding section/key in this basis. This information is currently only available for literature models.

```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> m = literature_model(arxiv_id = "1212.2949", equation = "3.2", model_parameters = Dict("k" => 5))
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(11) Tate model with parameter values (k = 5) based on arXiv paper 1212.2949 Eq. (3.2)

julia> classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes(m)
Dict{String, Vector{Int64}} with 6 entries:
  "b1" => [1, 0]
  "b3" => [3, -5]
  "b4" => [4, -6]
  "b6" => [6, -11]
  "b2" => [2, -1]
  "ζ0" => [0, 1]
```
"""
@attr Dict{String, Vector{Int}} function classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes(m::AbstractFTheoryModel)
  @req has_attribute(m, :classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes) "No detailed information about tunable sections stored for this model"
  return get_attribute(m, :classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes)
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
  return get_attribute(m, :vertices)::Vector{Vector{QQFieldElem}}
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
  return get_attribute(m, :poly_index)::Int
end


@doc raw"""
    has_quick_triangulation(m::AbstractFTheoryModel)

For a 3-dimensional reflexive polytope in the Kreuzer-Skarke list, the list of
full (sometimes also called fine), regular, star triangulations can be extremely
large. Consequently, one may wonder if the triangulations can be enumerated in a
somewhat reasonable time (say 5 minutes on a personal computer). This method tries
to provide an answer to this. It returns `true` if one should expect a timely response
to the attempt to enumerate all (full, regular, star) triangulations. Otherwise, this
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
  return get_attribute(m, :max_lattice_pts_in_facet)::Int
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
  return get_attribute(m, :estimated_number_of_triangulations)::Int
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
  return get_attribute(m, :Kbar3)::Int
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
  return get_attribute(m, :h11)::Int
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
  return get_attribute(m, :h12)::Int
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
  return get_attribute(m, :h13)::Int
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
  return get_attribute(m, :h22)::Int
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

julia> keys_list = collect(keys(genera_of_ci_curves(qsm_model)));

julia> my_key = only(filter(k -> string(k) == "x7", keys_list))
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

julia> keys_list = collect(keys(degrees_of_kbar_restrictions_to_ci_curves(qsm_model)));

julia> my_key = only(filter(k -> string(k) == "x7", keys_list))
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
  return get_attribute(m, :index_facet_interior_divisors)::Vector{Int}
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
given by Ci = V(xi, s). xi is a homogeneous coordinate of the 3-dimensional toric base space
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
    simplified_dual_graph(m::AbstractFTheoryModel)

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


######################################################################################
### (5) Attributes for flux families (not exported, rather for serialization overhaul)
######################################################################################

@attr QQMatrix function matrix_integral_quant_transverse(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_integral(special_flux_family(m, check = check))
end

@attr QQMatrix function matrix_rational_quant_transverse(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_rational(special_flux_family(m, check = check))
end

@attr Vector{QQFieldElem} function offset_quant_transverse(m::AbstractFTheoryModel; check::Bool = true)
  return offset(special_flux_family(m, check = check))
end

@attr QQMatrix function matrix_integral_quant_transverse_nobreak(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_integral(special_flux_family(m, not_breaking = true; check = check))
end

@attr QQMatrix function matrix_rational_quant_transverse_nobreak(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_rational(special_flux_family(m, not_breaking = true; check = check))
end

@attr Vector{QQFieldElem} function offset_quant_transverse_nobreak(m::AbstractFTheoryModel; check::Bool = true)
  return offset(special_flux_family(m, not_breaking = true; check = check))
end
