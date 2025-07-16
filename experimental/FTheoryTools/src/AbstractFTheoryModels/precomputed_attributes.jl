macro define_model_attribute_getter(arg_expr, doc_example="")
  if !(arg_expr isa Expr && arg_expr.head == :tuple && length(arg_expr.args) == 2)
    error("Expected input like: (function_name, ReturnType)")
  end

  fname_expr = arg_expr.args[1]
  rettype_expr = arg_expr.args[2]
  fname = fname_expr isa Symbol ? fname_expr : eval(fname_expr)
  sym = QuoteNode(fname)
  msg = "No $(replace(string(fname), '_' => ' ')) known for this model"

  default_doc = """
  Returns `$(fname)` of the F-theory model if known, otherwise throws an error.

  See [Literature Models](@ref) for more details.

  $doc_example
  """

  return quote
    @doc $default_doc
    function $(esc(fname))(m::AbstractFTheoryModel)
      @req has_attribute(m, $sym) $msg
      return get_attribute(m, $sym)::$(esc(rettype_expr))
    end
  end
end



##########################################
### (1) Meta data getters
##########################################

meta_data = [
    (:associated_literature_models, Vector{String}),
    (:arxiv_doi, String),
    (:arxiv_id, String),
    (:arxiv_link, String),
    (:arxiv_model_equation_number, String),
    (:arxiv_model_page, String),
    (:arxiv_model_section, String),
    (:arxiv_version, String),
    (:birational_literature_models, Vector{String}),
    (:journal_doi, String),
    (:journal_link, String),
    (:journal_model_equation_number, String),
    (:journal_model_page, String),
    (:journal_model_section, String),
    (:journal_pages, String),
    (:journal_report_numbers, Vector{String}),
    (:journal_volume, String),
    (:journal_name, String),
    (:journal_year, String),
    (:literature_identifier, String),
    (:model_description, String),
    (:model_parameters, Dict{String, Int64}),
    (:paper_authors, Vector{String}),
    (:paper_buzzwords, Vector{String}),
    (:paper_description, String),
    (:paper_title, String)
]

for (name, typ) in meta_data
    @eval @define_model_attribute_getter ($name, $typ)
end




##########################################
### (2) Geometric data getters
##########################################


# Return a list of the known Mordell–Weil generating sections of the given model.  If no generating sections are known, an error is raised.
@doc raw"""
    generating_sections(m::AbstractFTheoryModel)

Returns `generating_sections` of the F-theory model if known, otherwise throws an error.

See [Literature Models](@ref) for more details.    

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
  @req has_attribute(m, :generating_sections) "No generating sections known for this model"
  return get_attribute(m, :generating_sections)::GeneratingSectionsType
end



# Return a list of lists of known Mordell–Weil generating sections for the given model after each known resolution. Each element of the outer list corresponds to a known resolution (in the same order), and each element of the list associated to a given resolution corresponds to a known generating section (in the same order). If no resolution generating sections are known, an error is raised.
@define_model_attribute_getter((resolution_generating_sections, ResolutionGeneratingSectionsType),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> resolution_generating_sections(m)
1-element Vector{Vector{Vector{Vector{QQMPolyRingElem}}}}:
 [[[0, 0, 1], [0, 0, 1], [0, 1], [0, 1], [0, 1], [a32, -a43]]]
```
""")


# Return a list of known Mordell–Weil zero sections for the given model after each blowup of each known resolution. Each element of the list corresponds to a known resolution (in the same order). If no resolution zero sections are known, an error is raised.
@define_model_attribute_getter((resolution_zero_sections, ResolutionZeroSectionsType),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> resolution_zero_sections(m)
1-element Vector{Vector{Vector{QQMPolyRingElem}}}:
 [[1, 1, 0], [1, 1, w], [1, 1], [1, 1], [1, 1], [1, 1]]
```
""")



# Return a list of lists of known Mordell–Weil generating sections for the given model after each known weighted resolution. Each element of the outer list corresponds to a known weighted resolution (in the same order), and each element of the list associated to a given weighted resolution corresponds to a known generating section (in the same order). If no weighted resolution generating sections are known, an error is raised.
@define_model_attribute_getter((weighted_resolution_generating_sections, WeightedResolutionGeneratingSectionsType),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> weighted_resolution_generating_sections(m)
1-element Vector{Vector{Vector{Vector{QQMPolyRingElem}}}}:
 [[[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [a32, -a43]]]
```
""")



# Return a list of known Mordell–Weil zero sections for the given model after each known weighted resolution. Each element of the list corresponds to a known weighted resolution (in the same order). If no weighted resolution zero sections are known, an error is raised.
@define_model_attribute_getter((weighted_resolution_zero_sections, WeightedResolutionZeroSectionsType),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> weighted_resolution_zero_sections(m)
1-element Vector{Vector{Vector{QQMPolyRingElem}}}:
 [[1, 1, 0], [1, 1, w], [1, 1, w], [1, 1, w], [1, 1, w], [1, 1]]
```
""")



# Return the zero section of the given model. If no zero section is known, an error is raised. This information is not typically stored as an attribute for Weierstrass and global Tate models, whose zero sections are known.
@define_model_attribute_getter((zero_section, ZeroSectionType),
"""
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
""")


# Return the torsion sections of the given model. If no torsion sections are known, an error is raised.
@define_model_attribute_getter((torsion_sections, TorsionSectionsType),
"""
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
""")


# Return the list of all known resolutions for the given model. If no resolutions are known, an error is raised.
@define_model_attribute_getter((resolutions, Vector{Tuple{Vector{Vector{String}}, Vector{String}}}),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> resolutions(m)
1-element Vector{Tuple{Vector{Vector{String}}, Vector{String}}}:
 ([["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]], ["e1", "e4", "e2", "e3", "s"])
```
""")


# Return the list of all known weighted resolutions for the given model. If no weighted resolutions are known, an error is raised.
@define_model_attribute_getter((weighted_resolutions, Vector{Tuple{Vector{Tuple{Vector{String}, Vector{Int64}}}, Vector{String}}}),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> weighted_resolutions(m)
1-element Vector{Tuple{Vector{Tuple{Vector{String}, Vector{Int64}}}, Vector{String}}}:
 ([(["x", "y", "w"], [1, 1, 1]), (["x", "y", "w"], [1, 2, 1]), (["x", "y", "w"], [2, 2, 1]), (["x", "y", "w"], [2, 3, 1]), (["x", "y"], [1, 1])], ["e1", "e4", "e2", "e3", "s"])
```
""")


# Return the zero section class of a model as a cohomology class in the toric ambient space. If no zero section class is known, an error is raised. This information is always available for Weierstrass and global Tate models, whose zero section classes are known.
@define_model_attribute_getter((zero_section_class, CohomologyClass),
"""
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> zero_section_class(qsm_model)
Cohomology class on a normal toric variety given by e2 + 2*u + 3*e4 + e1 - w
```
""")


# Return the index of the generator of the Cox ring of the ambient space, whose corresponding vanishing locus defines the zero section of a model. If no zero section class is known, an error is raised. This attribute is always set simultaneously with zero_section_class. This information is always available for Weierstrass and global Tate models, whose zero section classes are known.
@define_model_attribute_getter((zero_section_index, Int),
"""
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
""")


# Return the gauge algebra of the given model. If no gauge algebra is known, an error is raised. This information is typically available for all models, however.
@define_model_attribute_getter((gauge_algebra, DirectSumLieAlgebra{QQBarFieldElem}),
"""
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
""")


# Return list of lists of matrices, where each list of matrices corresponds to a gauge factor of the same index given by `gauge_algebra(m)`. These matrices are elements of the center of the corresponding gauge factor and quotienting by them replicates the action of some discrete group on the center of the lie algebra. This list combined with `gauge_algebra(m)` completely determines the gauge group of the model. If no gauge quotients are known, an error is raised.
@define_model_attribute_getter((global_gauge_group_quotient, Vector{Vector{String}}),
"""
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
""")


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
  return get_attribute(m, :exceptional_divisor_indices, Vector{Int64}())
end
