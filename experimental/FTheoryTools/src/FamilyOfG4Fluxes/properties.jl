#####################################################
# 1 Basic properties
#####################################################

@doc raw"""
    is_well_quantized(fgs::FamilyOfG4Fluxes)

In case it is known if the family of G4-fluxes is well-quantized,
this method returns this boolean value -- true if well-quantized
and false if not well-quantized. In case it is not known if the
family of G4-fluxes is well-quantized, an error is raised.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gf = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: failed
  - Non-abelian gauge group: broken

julia> is_well_quantized(gf)
true
```
"""
function is_well_quantized(fgs::FamilyOfG4Fluxes)
  @req has_attribute(fgs, :is_well_quantized) "Cannot (yet) tell if this family of G4-fluxes is well-quantized"
  return get_attribute(fgs, :is_well_quantized)
end


@doc raw"""
    is_vertical(fgs::FamilyOfG4Fluxes)

In case it is known if the family of G4-fluxes is vertical,
this method returns this boolean value -- true if vertical
and false if not. In case it is not known if the family of
G4-fluxes is vertical, an error is raised.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gf = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: failed
  - Non-abelian gauge group: broken

julia> is_vertical(gf)
false

julia> gf2 = special_flux_family(qsm_model, vert = true, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: satisfied
  - Non-abelian gauge group: broken

julia> is_vertical(gf2)
true
```
"""
function is_vertical(fgs::FamilyOfG4Fluxes)
  @req has_attribute(fgs, :is_vertical) "Cannot (yet) tell if this family of G4-fluxes is vertical"
  return get_attribute(fgs, :is_vertical)
end


@doc raw"""
    breaks_non_abelian_gauge_group(fgs::FamilyOfG4Fluxes)

In case it is known if the family of G4-fluxes does break
the non-abelian gauge group, this method returns this boolean
value -- true if it break the non-abelian gauge group and false
if not. In case it is not known if the family of G4-fluxes does
break the non-abelian gauge group, an error is raised.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gf = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: failed
  - Non-abelian gauge group: broken

julia> breaks_non_abelian_gauge_group(gf)
true

julia> gf2 = special_flux_family(qsm_model, vert = true, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: satisfied
  - Non-abelian gauge group: broken

julia> breaks_non_abelian_gauge_group(gf2)
true

julia> gf3 = special_flux_family(qsm_model, vert = true, not_breaking = true, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: satisfied
  - Non-abelian gauge group: not broken

julia> breaks_non_abelian_gauge_group(gf3)
false
```
"""
function breaks_non_abelian_gauge_group(fgs::FamilyOfG4Fluxes)
  @req has_attribute(fgs, :breaks_non_abelian_gauge_group) "Cannot (yet) tell if this family of G4-fluxes breaks the non-abelian gauge group"
  return get_attribute(fgs, :breaks_non_abelian_gauge_group)
end
