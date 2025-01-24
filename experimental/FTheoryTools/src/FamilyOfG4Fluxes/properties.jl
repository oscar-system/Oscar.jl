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
  - Tadpole constraint: not analyzed

julia> is_well_quantized(gf)
true
```
"""
function is_well_quantized(fgs::FamilyOfG4Fluxes)
  @req has_attribute(fgs, :is_well_quantized) "Cannot (yet) tell if this family of G4-fluxes is well-quantized"
  return get_attribute(fgs, :is_well_quantized)
end


@doc raw"""
    is_vertical(fgs::FamilyOfG4Fluxes; check::Bool = true)

Checks if the given family of $G_4$-fluxes is vertical.
If so, this method returns `true` and otherwise `false`.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gf = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: failed
  - Non-abelian gauge group: broken
  - Tadpole constraint: not analyzed

julia> is_vertical(gf, check = false)
false

julia> gf2 = special_flux_family(qsm_model, vert = true, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: satisfied
  - Non-abelian gauge group: broken
  - Tadpole constraint: not analyzed

julia> is_vertical(gf2, check = false)
true

julia> m1 = matrix_integral(gf2);

julia> m2 = matrix_rational(gf2);

julia> gf3 = family_of_g4_fluxes(qsm_model, m1, m2, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole constraint: not analyzed

julia> is_vertical(gf3)
true
```
"""
function is_vertical(fgs::FamilyOfG4Fluxes; check::Bool = true)
  # Entry checks
  m = model(fgs)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Verticality check only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Verticality check currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Verticality check currently supported only for toric ambient space"
  
  # Is the result known?
  if has_attribute(fgs, :is_vertical)
    return get_attribute(fgs, :is_vertical)
  end

  # Extract ambient space model of g4-fluxes, in terms of which we express the generators of the flux family
  mb = ambient_space_models_of_g4_fluxes(model(fgs), check = check)
  nmb = length(mb)

  # Verify that each generator of the flux family is vertical
  my_mat = matrix_integral(fgs)
  for k in 1:ncols(my_mat)
    class = sum(my_mat[l,k] * mb[l] for l in 1:nmb)
    gen_k = g4_flux(model(fgs), class)
    if !passes_verticality_checks(gen_k)
      set_attribute!(fgs, :is_vertical, false)
      return false
    end
  end
  my_mat = matrix_rational(fgs)
  for k in 1:ncols(my_mat)
    class = sum(my_mat[l,k] * mb[l] for l in 1:nmb)
    gen_k = g4_flux(model(fgs), class)
    if !passes_verticality_checks(gen_k)
      set_attribute!(fgs, :is_vertical, false)
      return false
    end
  end
  set_attribute!(fgs, :is_vertical, true)
  return true
end


@doc raw"""
    breaks_non_abelian_gauge_group(fgs::FamilyOfG4Fluxes; check::Bool = true)

Checks if a family of G4-fluxes breaks the non-abelian
gauge group. If so, this method returns `true` and
otherwise `false`.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gf = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: failed
  - Non-abelian gauge group: broken
  - Tadpole constraint: not analyzed

julia> breaks_non_abelian_gauge_group(gf)
true

julia> gf2 = special_flux_family(qsm_model, vert = true, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: satisfied
  - Non-abelian gauge group: broken
  - Tadpole constraint: not analyzed

julia> breaks_non_abelian_gauge_group(gf2)
true

julia> gf3 = special_flux_family(qsm_model, vert = true, not_breaking = true, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: satisfied
  - Non-abelian gauge group: not broken
  - Tadpole constraint: not analyzed

julia> breaks_non_abelian_gauge_group(gf3)
false

julia> m1 = matrix_integral(gf3);

julia> m2 = matrix_rational(gf3);

julia> gf4 = family_of_g4_fluxes(qsm_model, m1, m2, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole constraint: not analyzed

julia> breaks_non_abelian_gauge_group(gf4, check = false)
false
```
"""
function breaks_non_abelian_gauge_group(fgs::FamilyOfG4Fluxes; check::Bool = true)
  # Entry checks
  m = model(fgs)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Gauge group breaking check only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Gauge group breaking check currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Gauge group breaking check currently supported only for toric ambient space"
  
  # Is the result known?
  if has_attribute(fgs, :breaks_non_abelian_gauge_group)
    return get_attribute(fgs, :breaks_non_abelian_gauge_group)
  end

  # Extract ambient space model of g4-fluxes, in terms of which we express the generators of the flux family
  mb = ambient_space_models_of_g4_fluxes(model(fgs), check = check)
  nmb = length(mb)

  # Verify that each generator of the flux family does not break the non-abelian gauge group
  my_mat = matrix_integral(fgs)
  for k in 1:ncols(my_mat)
    class = sum(my_mat[l,k] * mb[l] for l in 1:nmb)
    gen_k = g4_flux(model(fgs), class)
    if breaks_non_abelian_gauge_group(gen_k)
      set_attribute!(fgs, :breaks_non_abelian_gauge_group, true)
      return true
    end
  end
  my_mat = matrix_rational(fgs)
  for k in 1:ncols(my_mat)
    class = sum(my_mat[l,k] * mb[l] for l in 1:nmb)
    gen_k = g4_flux(model(fgs), class)
    if breaks_non_abelian_gauge_group(gen_k)
      set_attribute!(fgs, :breaks_non_abelian_gauge_group, true)
      return true
    end
  end
  set_attribute!(fgs, :breaks_non_abelian_gauge_group, false)
  return false
end
