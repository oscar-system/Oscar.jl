################
# 1: Constructor
################

@doc raw"""
    family_of_g4_fluxes(m::AbstractFTheoryModel, mat_int::QQMatrix, mat_rat::QQMatrix; check::Bool = true)

Given an F-theory model with a toric ambient space, we can 
identify ambient space candidates of G4-fluxes. In terms of these
candidates, we can define a family of G4-fluxes as:
- $\mathbb{Z}$-linear combinations, provided by a matrix $mat_int$,
- $\mathbb{Q}$-linear combinations, provided by a matrix $mat_rat$.
An example is in order.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> family_of_g4_fluxes(qsm_model, mat_int, mat_rat, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
```
"""
function family_of_g4_fluxes(m::AbstractFTheoryModel, mat_int::QQMatrix, mat_rat::QQMatrix; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Family of G4-fluxes only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Family of G4-flux currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Family of G4-flux currently supported only for toric ambient space"
  @req nrows(mat_int) == nrows(mat_rat) "Number of rows in both matrices must coincide"
  n_gens = length(chosen_g4_flux_basis(m, check = check))
  @req nrows(mat_int) == n_gens "Number of rows in both matrices must agree with the number of ambient space models of G4-fluxes"
  return FamilyOfG4Fluxes(m, mat_int, mat_rat)
end


################################################
# 2: Equality and hash
################################################

function Base.:(==)(gf1::FamilyOfG4Fluxes, gf2::FamilyOfG4Fluxes)
  return model(gf1) === model(gf2) && mat_int(gf1) == mat_int(gf2) && mat_rat(gf1) == mat_rat(gf2)
end

function Base.hash(gf::FamilyOfG4Fluxes, h::UInt)
  b = 0x92bd6ac4f87d834e % UInt
  h = hash(model(gf), h)
  h = hash(mat_int(gf), h)
  h = hash(mat_rat(gf), h)
  return xor(h, b)
end


################################################
# 3: Display
################################################

function Base.show(io::IO, gf::FamilyOfG4Fluxes)
  properties_string = ["A family of G4 fluxes:"]

  # Check for elementary quantization checks
  if has_attribute(gf, :is_well_quantized) && get_attribute(gf, :is_well_quantized) !== nothing
    if is_well_quantized(gf)
      push!(properties_string, "  - Elementary quantization checks: satisfied")
    else
      push!(properties_string, "  - Elementary quantization checks: failed")
    end
  else
    push!(properties_string, "  - Elementary quantization checks: not executed")
  end

  # Check for transversality checks
  if has_attribute(gf, :passes_transversality_checks) && get_attribute(gf, :passes_transversality_checks) !== nothing
    if passes_transversality_checks(gf)
      push!(properties_string, "  - Transversality checks: satisfied")
    else
      push!(properties_string, "  - Transversality checks: failed")
    end
  else
    push!(properties_string, "  - Transversality checks: not executed")
  end

  # Check for non-abelian gauge group breaking
  if has_attribute(gf, :breaks_non_abelian_gauge_group) && get_attribute(gf, :breaks_non_abelian_gauge_group) !== nothing
    if breaks_non_abelian_gauge_group(gf)
      push!(properties_string, "  - Non-abelian gauge group: broken")
    else
      push!(properties_string, "  - Non-abelian gauge group: not broken")
    end
  else
    push!(properties_string, "  - Non-abelian gauge group: breaking pattern not analyzed")
  end

  # Is the tadpole constrained worked out as polynomial?
  #=
  if has_attribute(gf, :d3_tadpole_constraint) && get_attribute(gf, :d3_tadpole_constraint) !== nothing
    push!(properties_string, "  - Tadpole constraint: evaluated")
  else
    push!(properties_string, "  - Tadpole constraint: not analyzed")
  end
  =#

  # Print each line separately, to avoid extra line break at the end
  for (i, line) in enumerate(properties_string)
    if i == length(properties_string)
      print(io, line) # Last line without extra newline
    else
      println(io, line) # Print all other lines with line break
    end
  end

end
