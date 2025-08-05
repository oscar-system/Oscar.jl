##########################################
### (1) Meta Data
##########################################

@doc raw"""
    add_associated_literature_model!(m::AbstractFTheoryModel, addition::String)

Add a new entry to the list of associated literature models for the F-theory model.
If the entry is already present, nothing is changed.

See [Metadata Attributes](@ref meta_data_attributes) for more details.

# Examples
```jldoctest
julia> m = literature_model(arxiv_id = "1408.4808", equation = "3.168", type = "hypersurface")
Hypersurface model over a not fully specified base

julia> add_associated_literature_model!(m, "1408_4808-1")

julia> length(associated_literature_models(m))
15
```
"""
function add_associated_literature_model!(m::AbstractFTheoryModel, addition::String)
  known_values = get_attribute(m, :associated_literature_models, String[])
  if !(addition in known_values)
    set_attribute!(m, :associated_literature_models => vcat(known_values, [addition]))
  end
  return nothing
end


@doc raw"""
    add_birational_literature_model!(m::AbstractFTheoryModel, addition::String)

Add a new entry to the list of birational models for the F-theory model.
If the entry is already present, nothing is changed.

See [Metadata Attributes](@ref meta_data_attributes) for more details.

# Examples
```jldoctest
julia> m = literature_model(arxiv_id = "1408.4808", equation = "3.168", type = "hypersurface")
Hypersurface model over a not fully specified base

julia> add_birational_literature_model!(m, "1408_4808-14-WSF")

julia> length(birational_literature_models(m))
1
```
"""
function add_birational_literature_model!(m::AbstractFTheoryModel, addition::String)
  known_values = get_attribute(m, :birational_literature_models, String[])
  if !(addition in known_values)
    set_attribute!(m, :birational_literature_models => vcat(known_values, [addition]))
  end
  return nothing
end


@doc raw"""
    add_journal_report_number!(m::AbstractFTheoryModel, addition::String)

Add a new entry to the list of journal report numbers for the F-theory model.
If the entry is already present, nothing is changed.

See [Metadata Attributes](@ref meta_data_attributes) for more details.

# Examples
```jldoctest
julia> m = literature_model(arxiv_id = "1408.4808", equation = "3.168", type = "hypersurface")
Hypersurface model over a not fully specified base

julia> add_journal_report_number!(m, "UPR-1264-T")

julia> length(journal_report_numbers(m))
3
```
"""
function add_journal_report_number!(m::AbstractFTheoryModel, addition::String)
  known_values = get_attribute(m, :journal_report_numbers, String[])
  if !(addition in known_values)
    return set_attribute!(m, :journal_report_numbers => vcat(known_values, [addition]))
  end
  return nothing
end


@doc raw"""
    add_paper_author!(m::AbstractFTheoryModel, addition::String)

Add a new entry to the list of paper authors for the F-theory model.
If the entry is already present, nothing is changed.

See [Metadata Attributes](@ref meta_data_attributes) for more details.

# Examples
```jldoctest
julia> m = literature_model(arxiv_id = "1408.4808", equation = "3.168", type = "hypersurface")
Hypersurface model over a not fully specified base

julia> add_paper_author!(m, "Denis Klevers")

julia> length(paper_authors(m))
5
```
"""
function add_paper_author!(m::AbstractFTheoryModel, addition::String)
  known_values = get_attribute(m, :paper_authors, String[])
  if !(addition in known_values)
    return set_attribute!(m, :paper_authors => vcat(known_values, [addition]))
  end
  return nothing
end


@doc raw"""
    add_paper_buzzword!(m::AbstractFTheoryModel, addition::String)

Add a new entry to the list of paper buzzwords for the F-theory model.
If the entry is already present, nothing is changed.

See [Metadata Attributes](@ref meta_data_attributes) for more details.

```jldoctest
julia> m = literature_model(arxiv_id = "1408.4808", equation = "3.168", type = "hypersurface")
Hypersurface model over a not fully specified base

julia> add_paper_buzzword!(m, "Mordell-Weil")

julia> length(paper_buzzwords(m))
4
```
"""
function add_paper_buzzword!(m::AbstractFTheoryModel, addition::String)
  known_values = get_attribute(m, :paper_buzzwords, String[])
  if !(addition in known_values)
    return set_attribute!(m, :paper_buzzwords => vcat(known_values, [addition]))
  end
  return nothing
end



##########################################
### (2) Zero Sections
##########################################



##########################################
### (3) Mordell-Weil Group
##########################################

@doc raw"""
    add_generating_section!(m::AbstractFTheoryModel, addition::Vector{String})

Add a generating section for a model.

!!! warning
    The newly added section must be specified in terms of the base space coordinates only.

See [Mordell–Weil Group](@ref mordell_weil_group_data) for more details.

# Examples
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> add_generating_section!(m, ["0","0","1"])

julia> length(generating_sections(m))
1
```
"""
function add_generating_section!(m::AbstractFTheoryModel, addition::Vector{String})
  R = parent(first(values(explicit_model_sections(m))))
  new_entry = deepmap(s -> eval_poly(s, R), addition)
  known_values =  get_attribute(m, :generating_sections, GeneratingSectionsType())
  if !(new_entry in known_values)
    set_attribute!(m, :generating_sections => vcat(known_values, [new_entry]))
  end
  return nothing
end


@doc raw"""
    add_torsion_section!(m::AbstractFTheoryModel, addition::Vector{String})

Add a torsion section for a model.

!!! warning
    The newly added section must be specified in terms of the base space coordinates only.

See [Mordell–Weil Group](@ref mordell_weil_group_data) for more details.

# Examples
```jldoctest
julia> m = literature_model(arxiv_id = "1408.4808", equation = "3.190", type = "hypersurface")
Hypersurface model over a not fully specified base

julia> add_torsion_section!(m, ["1", "s5", "-s2", "1", "1", "1", "1", "0"])

julia> length(torsion_sections(m))
1
```
"""
function add_torsion_section!(m::AbstractFTheoryModel, addition::Vector{String})
  R = parent(first(values(explicit_model_sections(m))))
  new_entry = deepmap(s -> eval_poly(s, R), addition)
  known_values = get_attribute(m, :torsion_sections, TorsionSectionsType())
  if !(new_entry in known_values)
    set_attribute!(m, :torsion_sections => vcat(known_values, [new_entry]))
  end
  return nothing
end



##########################################
### (4) Gauge Group
##########################################



##########################################
### (5) (Weighted) Resolutions
##########################################

@doc raw"""
    add_resolution!(m::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})

Add a resolution for a model.

See [Registering And Extracting Known Resolution Sequences](@ref working_with_resolution_sequences) for more details.

# Examples
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> add_resolution!(m, [["x", "y", "w"], ["y", "e1"], ["x", "e4"], ["y", "e2"], ["x", "y"]], ["e1", "e4", "e2", "e3", "s"])

julia> length(resolutions(m))
1
```
"""
function add_resolution!(m::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})
  @req length(exceptionals) == length(centers) "Number of exceptionals must match number of centers"
  new_entry = (centers, exceptionals)
  known_resolutions = get_attribute(m, :resolutions, Tuple{Vector{Vector{String}}, Vector{String}}[])
  if !(new_entry in known_resolutions)
    set_attribute!(m, :resolutions => vcat(known_resolutions, [new_entry]))
  end
  return nothing
end


@doc raw"""
    add_weighted_resolution!(m::AbstractFTheoryModel, centers::Vector{Tuple{Vector{String}, Vector{Int64}}}, exceptionals::Vector{String})

Add a weighted resolution for a model.

See [Registering And Extracting Known Resolution Sequences](@ref working_with_resolution_sequences) for more details.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, model_sections = Dict("w" => w), completeness_check = false)
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> centers = [(["x", "y", "w"], [1, 1, 1]), (["x", "y", "w"], [1, 2, 1]), (["x", "y", "w"], [2, 2, 1]), (["x", "y", "w"], [2, 3, 1]), (["x", "y"], [1, 1])];

julia> exceptionals = ["e1", "e4", "e2", "e3", "s"];

julia> add_weighted_resolution!(m, centers, exceptionals)

julia> length(weighted_resolutions(m))
1
```
"""
function add_weighted_resolution!(m::AbstractFTheoryModel, centers::Vector{Tuple{Vector{String}, Vector{Int64}}}, exceptionals::Vector{String})
  @req length(exceptionals) == length(centers) "Number of exceptionals must match number of centers"
  new_entry = (centers, exceptionals)
  known_values = get_attribute(m, :weighted_resolutions, Tuple{Vector{Tuple{Vector{String}, Vector{Int64}}}, Vector{String}}[])
  if !(new_entry in known_values)
    set_attribute!(m, :weighted_resolutions => vcat(known_values, [new_entry]))
  end
  return nothing
end



##########################################
### (6) Resolution Metadata Functions
##########################################

@doc raw"""
    add_resolution_zero_section!(m::AbstractFTheoryModel, addition::Vector{Vector{String}})

Add a resolution zero section for a model.

!!! warning
    The newly added section must be specified in terms of the base space coordinates only.

See [Resolution Metadata Functions](@ref resolution_meta_data) for more details.

# Examples
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> add_resolution_zero_section!(m, [["1", "1", "0"], ["1", "1", "w"], ["1", "1"], ["1", "1"], ["1", "1"], ["1", "1"]])

julia> length(resolution_zero_sections(m))
1
```
"""
function add_resolution_zero_section!(m::AbstractFTheoryModel, addition::Vector{Vector{String}})
  R = parent(first(values(explicit_model_sections(m))))
  new_entry = deepmap(s -> eval_poly(s, R), addition)
  known_values = get_attribute(m, :resolution_zero_sections, ResolutionZeroSectionsType())
  if !(new_entry in known_values)
    set_attribute!(m, :resolution_zero_sections => vcat(known_values, [new_entry]))
  end
  return nothing
end


@doc raw"""
    add_resolution_generating_section!(m::AbstractFTheoryModel, addition::Vector{Vector{Vector{String}}})

Add a resolution generating section for a model.

!!! warning
    The newly added section must be specified in terms of the base space coordinates only.

See [Resolution Metadata Functions](@ref resolution_meta_data) for more details.

# Examples
```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> add_resolution_generating_section!(m, [[["0", "0", "1"], ["0", "0", "1"], ["0", "1"], ["0", "1"], ["0", "1"], ["a32", "a43"]]])

julia> length(resolution_generating_sections(m))
1
```
"""
function add_resolution_generating_section!(m::AbstractFTheoryModel, addition::Vector{Vector{Vector{String}}})
  R = parent(first(values(explicit_model_sections(m))))
  new_entry = deepmap(s -> eval_poly(s, R), addition)
  known_values = get_attribute(m, :resolution_generating_sections, ResolutionGeneratingSectionsType())
  if !(new_entry in known_values)
    set_attribute!(m, :resolution_generating_sections => vcat(known_values, [new_entry]))
  end
  return nothing
end


@doc raw"""
    add_weighted_resolution_zero_section!(m::AbstractFTheoryModel, addition::Vector{Vector{String}})

Add a weighted resolution zero section for a model.

!!! warning
    The newly added section must be specified in terms of the base space coordinates only.

See [Resolution Metadata Functions](@ref resolution_meta_data) for more details.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, model_sections = Dict("w" => w), completeness_check = false)
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> add_weighted_resolution_zero_section!(m, [["1", "1", "0"], ["1", "1", "x1"], ["1", "1", "x1"], ["1", "1", "x1"], ["1", "1", "x1"], ["1", "1"]])

julia> length(weighted_resolution_zero_sections(m))
1
```
"""
function add_weighted_resolution_zero_section!(m::AbstractFTheoryModel, addition::Vector{Vector{String}})
  R = parent(first(values(explicit_model_sections(m))))
  new_entry = deepmap(s -> eval_poly(s, R), addition)
  known_values = get_attribute(m, :weighted_resolution_zero_sections, WeightedResolutionZeroSectionsType())
  if !(new_entry in known_values)
    set_attribute!(m, :weighted_resolution_zero_sections => vcat(known_values, [new_entry]))
  end
  return nothing
end


@doc raw"""
    add_weighted_resolution_generating_section!(m::AbstractFTheoryModel, addition::Vector{Vector{Vector{String}}})

Add a weighted resolution generating section for a model.

!!! warning
    The newly added section must be specified in terms of the base space coordinates only.

See [Resolution Metadata Functions](@ref resolution_meta_data) for more details.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, model_sections = Dict("w" => w), completeness_check = false)
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> addition = [[["0", "0", "1"], ["0", "0", "1"], ["0", "0", "1"], ["0", "0", "1"], ["0", "0", "1"], ["1980*x1^10", "0", "0"]]];

julia> add_weighted_resolution_generating_section!(m, addition)

julia> length(weighted_resolution_generating_sections(m))
2
```
"""
function add_weighted_resolution_generating_section!(m::AbstractFTheoryModel, addition::Vector{Vector{Vector{String}}})
  R = parent(first(values(explicit_model_sections(m))))
  new_entry = deepmap(s -> eval_poly(s, R), addition)
  known_values = get_attribute(m, :weighted_resolution_generating_sections, WeightedResolutionGeneratingSectionsType())
  if !(new_entry in known_values)
    set_attribute!(m, :weighted_resolution_generating_sections => vcat(known_values, [new_entry]))
  end
  return nothing
end



##########################################
### (7) Exceptional Divisors
##########################################
