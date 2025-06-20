##########################################
### (1) General properties
##########################################

@doc raw"""
    is_base_space_fully_specified(m::AbstractFTheoryModel)

Return `true` if the F-theory model has a concrete base space and `false` otherwise.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> is_base_space_fully_specified(t)
false
```
"""
is_base_space_fully_specified(m::AbstractFTheoryModel) = !(m.base_space isa FamilyOfSpaces)


@doc raw"""
    is_partially_resolved(m::AbstractFTheoryModel)

Return `true` if resolution techniques were applied to the F-theory model,
thereby potentially resolving its singularities. Otherwise, return `false`.

```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> is_partially_resolved(t)
false

julia> t2 = blow_up(t, ["x", "y", "x1"]; coordinate_name = "e1")
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> is_partially_resolved(t2)
true
```
"""
is_partially_resolved(m::AbstractFTheoryModel) = get_attribute(m, :partially_resolved)::Bool



##########################################
### (2) Meta-data properties
##########################################

has_arxiv_id(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_id)
has_arxiv_doi(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_doi)
has_arxiv_link(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_link)
has_arxiv_model_equation_number(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_model_equation_number)
has_arxiv_model_page(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_model_page)
has_arxiv_model_section(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_model_section)
has_arxiv_version(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_version)
has_associated_literature_models(m::AbstractFTheoryModel) = has_attribute(m, :associated_literature_models)
has_journal_doi(m::AbstractFTheoryModel) = has_attribute(m, :journal_doi)
has_journal_link(m::AbstractFTheoryModel) = has_attribute(m, :journal_link)
has_journal_model_equation_number(m::AbstractFTheoryModel) = has_attribute(m, :journal_model_equation_number)
has_journal_model_page(m::AbstractFTheoryModel) = has_attribute(m, :journal_model_page)
has_journal_model_section(m::AbstractFTheoryModel) = has_attribute(m, :journal_model_section)
has_journal_name(m::AbstractFTheoryModel) = has_attribute(m, :journal_name)
has_journal_pages(m::AbstractFTheoryModel) = has_attribute(m, :journal_pages)
has_journal_report_numbers(m::AbstractFTheoryModel) = has_attribute(m, :journal_report_numbers)
has_journal_volume(m::AbstractFTheoryModel) = has_attribute(m, :journal_volume)
has_journal_year(m::AbstractFTheoryModel) = has_attribute(m, :journal_year)
has_literature_identifier(m::AbstractFTheoryModel) = has_attribute(m, :literature_identifier)
has_model_description(m::AbstractFTheoryModel) = has_attribute(m, :model_description)
has_model_parameters(m::AbstractFTheoryModel) = has_attribute(m, :model_parameters)
has_paper_authors(m::AbstractFTheoryModel) = has_attribute(m, :paper_authors)
has_paper_buzzwords(m::AbstractFTheoryModel) = has_attribute(m, :paper_buzzwords)
has_paper_description(m::AbstractFTheoryModel) = has_attribute(m, :paper_description)
has_paper_title(m::AbstractFTheoryModel) = has_attribute(m, :paper_title)
has_birational_literature_models(m::AbstractFTheoryModel) = has_attribute(m, :birational_literature_models)



##########################################
### (3) Specialized model properties
##########################################

has_generating_sections(m::AbstractFTheoryModel) = has_attribute(m, :generating_sections)
has_resolutions(m::AbstractFTheoryModel) = has_attribute(m, :resolutions)
has_resolution_generating_sections(m::AbstractFTheoryModel) = has_attribute(m, :resolution_generating_sections)
has_resolution_zero_sections(m::AbstractFTheoryModel) = has_attribute(m, :resolution_zero_sections)
has_weighted_resolutions(m::AbstractFTheoryModel) = has_attribute(m, :weighted_resolutions)
has_weighted_resolution_generating_sections(m::AbstractFTheoryModel) = has_attribute(m, :weighted_resolution_generating_sections)
has_weighted_resolution_zero_sections(m::AbstractFTheoryModel) = has_attribute(m, :weighted_resolution_zero_sections)
has_zero_section(m::AbstractFTheoryModel) = has_attribute(m, :zero_section)
has_zero_section_class(m::AbstractFTheoryModel) = has_attribute(m, :zero_section_class)
has_torsion_sections(m::AbstractFTheoryModel) = has_attribute(m, :torsion_sections)
has_gauge_algebra(m::AbstractFTheoryModel) = has_attribute(m, :gauge_algebra)
has_global_gauge_group_quotient(m::AbstractFTheoryModel) = has_attribute(m, :global_gauge_quotients)



##########################################
### (4) Consistency checks
##########################################

@doc raw"""
    verify_euler_characteristic_from_hodge_numbers(m::AbstractFTheoryModel; check::Bool = true)

Verify if the Euler characteristic, as computed from integrating the 4-th Chern class,
agrees with the results obtained from using the alternating sum of the Hodge numbers.
If so, this method returns `true`. However, should information be missing, (e.g. some
Hodge numbers), or the dimension of the F-theory model differ form 4, then this method
raises an error.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> verify_euler_characteristic_from_hodge_numbers(qsm_model; check = false)
true
```
"""
@attr Bool function verify_euler_characteristic_from_hodge_numbers(m::AbstractFTheoryModel; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Verification of Euler characteristic of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Verification of Euler characteristic of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Verification of Euler characteristic of F-theory model currently supported only for toric ambient space"
  @req dim(base_space(m)) == 3 "Verification of Euler characteristic of F-theory model currently supported only for toric base spaces of dimension 3"
  @req dim(ambient_space(m)) == 5 "Verification of Euler characteristic of F-theory model currently supported only for toric ambient spaces of dimension 5"
  @req has_attribute(m, :h11) "Verification of Euler characteristic of F-theory model requires h11"
  @req has_attribute(m, :h12) "Verification of Euler characteristic of F-theory model requires h12"
  @req has_attribute(m, :h13) "Verification of Euler characteristic of F-theory model requires h13"
  @req has_attribute(m, :h22) "Verification of Euler characteristic of F-theory model requires h22"

  # Computer Euler characteristic from integrating c4
  ec = euler_characteristic(m; check)

  # Compute Euler characteristic from adding Hodge numbers
  ec2 = 4 + 2 * hodge_h11(m) - 4 * hodge_h12(m) + 2 * hodge_h13(m) + hodge_h22(m)

  # Compute result of verification
  return ec == ec2
end


@doc raw"""
    is_calabi_yau(m::AbstractFTheoryModel; check::Bool = true)

Verify if the first Chern class of the tangent bundle of the F-theory geometry
$Y_n$ vanishes. If so, this confirms that this geometry is indeed Calabi-Yau,
as required by the reasoning of F-theory.

The implemented algorithm works for hypersurface, Weierstrass and global Tate models,
which are defined in a toric ambient space. It expresses $c_1(Y_n)$ as the restriction
of a cohomology class $h$ on the toric ambient space. This in turn requires that the
toric ambient space is simplicial and complete. We provide a switch to turn off 
these computationally very demanding checks. This is demonstrated in the example below.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> is_calabi_yau(qsm_model, check = false)
true
```
"""
@attr Bool function is_calabi_yau(m::AbstractFTheoryModel; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Verification of Euler characteristic of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Verification of Euler characteristic of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Verification of Euler characteristic of F-theory model currently supported only for toric ambient space"
  return is_trivial(chern_class(m, 1; check))
end
