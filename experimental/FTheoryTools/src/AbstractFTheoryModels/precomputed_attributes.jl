##########################################
### (1) Meta data getters
##########################################

@define_model_attribute_getter((associated_literature_models, Vector{String}), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((arxiv_doi, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((arxiv_id, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((arxiv_link, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((arxiv_model_equation_number, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((arxiv_model_page, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((arxiv_model_section, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((arxiv_version, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((birational_literature_models, Vector{String}), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((journal_doi, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((journal_link, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((journal_model_equation_number, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((journal_model_page, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((journal_model_section, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((journal_pages, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((journal_report_numbers, Vector{String}), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((journal_volume, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((journal_name, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((journal_year, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((literature_identifier, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((model_description, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((model_parameters, Dict{String, Int64}), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((paper_authors, Vector{String}), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((paper_buzzwords, Vector{String}), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((paper_description, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")
@define_model_attribute_getter((paper_title, String), "", "See [Metadata Attributes](@ref meta_data_attributes) for more details.")




##########################################
### (2) Geometric data getters
##########################################


@doc raw"""
    exceptional_classes(m::AbstractFTheoryModel)

Return the cohomology classes of the exceptional toric divisors of a model as a vector of cohomology classes in the toric ambient space.
This information is only supported for models over a concrete base that is a normal toric variety, but is always available in this case.
After a toric blow up this information is updated.

# Examples
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

# Examples
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
""", "See [Advanced Mathematical Attributes](@ref gauge_group_data) for more details.")


@define_model_attribute_getter((generating_sections, GeneratingSectionsType),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> generating_sections(m)
1-element Vector{Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 [0, 0, 1]
```
""", "See [Advanced Mathematical Attributes](@ref mordell_weil_group_data) for more details.")


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
""", "See [Advanced Mathematical Attributes](@ref gauge_group_data) for more details.")


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
""", "See [Advanced Mathematical Attributes](@ref resolving_f_theory_models) for more details.")


@define_model_attribute_getter((resolution_generating_sections, ResolutionGeneratingSectionsType),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> resolution_generating_sections(m)
1-element Vector{Vector{Vector{Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}}:
 [[[0, 0, 1], [0, 0, 1], [0, 1], [0, 1], [0, 1], [a32, a43]]]
```
""", "See [Advanced Mathematical Attributes](@ref resolving_f_theory_models) for more details.")


@define_model_attribute_getter((resolution_zero_sections, ResolutionZeroSectionsType),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> resolution_zero_sections(m)
1-element Vector{Vector{Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}:
 [[1, 1, 0], [1, 1, w], [1, 1], [1, 1], [1, 1], [1, 1]]
```
""", "See [Advanced Mathematical Attributes](@ref resolving_f_theory_models) for more details.")


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
""", "See [Advanced Mathematical Attributes](@ref mordell_weil_group_data) for more details.")


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
""", "See [Advanced Mathematical Attributes](@ref resolving_f_theory_models) for more details.")


@define_model_attribute_getter((weighted_resolution_generating_sections, WeightedResolutionGeneratingSectionsType),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> weighted_resolution_generating_sections(m)
1-element Vector{Vector{Vector{Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}}:
 [[[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [a32, a43]]]
```
""", "See [Advanced Mathematical Attributes](@ref resolving_f_theory_models) for more details.")


@define_model_attribute_getter((weighted_resolution_zero_sections, WeightedResolutionZeroSectionsType),
"""
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> weighted_resolution_zero_sections(m)
1-element Vector{Vector{Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}:
 [[1, 1, 0], [1, 1, w], [1, 1, w], [1, 1, w], [1, 1, w], [1, 1]]
```
""", "See [Advanced Mathematical Attributes](@ref resolving_f_theory_models) for more details.")


@define_model_attribute_getter((zero_section, ZeroSectionType),
"""
```jldoctest
julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5")
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base

julia> zero_section(h)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 0
 1
 0
```
""", "See [Advanced Mathematical Attributes](@ref zero_section_data) for more details.")


@define_model_attribute_getter((zero_section_class, CohomologyClass),
"""
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> zero_section_class(qsm_model)
Cohomology class on a normal toric variety given by e2 + 2*u + 3*e4 + e1 - w
```
""", "See [Advanced Mathematical Attributes](@ref zero_section_data) for more details.")


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
""", "See [Advanced Mathematical Attributes](@ref zero_section_data) for more details.")
