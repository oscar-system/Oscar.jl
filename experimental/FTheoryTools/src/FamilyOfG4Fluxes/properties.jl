#####################################################
# 1 Basic properties
#####################################################

@doc raw"""
    is_well_quantized(fgs::FamilyOfG4Fluxes; check::Bool = true)

Executes elementary tests (they are necessary but not sufficient)
to tell if a family of $G_4$-fluxes is well-quantized. In case
any of these tests fails, we know that this family of $G_4$-fluxes
is definitely not well-quantized. This method then returns `false`.

In the opposite case that all elementary tests pass, this method
returns `true`. Note however that this does not imply that the
family is well-quantized. At least, according to the current
theoretical understanding it does not. Rather, this means that
to the best of our current means, this family of fluxes appears
to be well-quantized.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gf = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: broken

julia> is_well_quantized(gf)
true

julia> m1 = matrix_integral(gf);

julia> m2 = matrix_rational(gf);

julia> gf2 = family_of_g4_fluxes(qsm_model, m1, m2, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> is_well_quantized(gf2, check = false)
true
```
"""
@attr Bool function is_well_quantized(fgs::FamilyOfG4Fluxes; check::Bool = true)
  # Entry checks
  m = model(fgs)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Elementary quantization check only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Elementary quantization checks currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Elementary quantization checks currently supported only for toric ambient space"
  @req all(==(0), offset(fgs)) "Currently, is_well_quantized is only supported for flux families with trivial offset"
  # TODO: Remove this limitation, i.e. support this functionality for all flux families!

  # Extract ambient space model of g4-fluxes, in terms of which we express the generators of the flux family
  mb = chosen_g4_flux_basis(model(fgs), check = check)
  nmb = length(mb)

  # Verify that each integral generator is well-quantized
  my_mat = matrix_integral(fgs)
  for k in 1:ncols(my_mat)
    gen_k = sum(my_mat[l,k] * mb[l] for l in 1:nmb)
    if !is_well_quantized(gen_k)
      return false
    end
  end

  # Verify that each rational generator is well-quantized, in that all relevant integrals vanish.
  cy = polynomial(cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))))
  c_ds = [polynomial(cohomology_class(d)) for d in torusinvariant_prime_divisors(ambient_space(m))]
  my_mat = matrix_rational(fgs)
  for k in 1:ncols(my_mat)
    gen_k = sum(my_mat[l,k] * mb[l] for l in 1:nmb)
    twist_g4 = polynomial(gen_k + 1//2 * chern_class(m, 2; check = false))
    for i in 1:length(c_ds)
      for j in i:length(c_ds)
        numb = integrate(cohomology_class(ambient_space(m), twist_g4 * c_ds[i] * c_ds[j] * cy); check = false)
        if !is_zero(numb)
          return false    
        end
      end
    end
  end

  # All other tests passed, so must be well-quantized according to elementary tests.
  return true
end


@doc raw"""
    passes_transversality_checks(fgs::FamilyOfG4Fluxes; check::Bool = true)

Check if the given family of $G_4$-fluxes passes the transversality checks.
If so, this method returns `true` and otherwise `false`.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gf = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: broken

julia> passes_transversality_checks(gf, check = false)
true

julia> m1 = matrix_integral(gf);

julia> m2 = matrix_rational(gf);

julia> gf3 = family_of_g4_fluxes(qsm_model, m1, m2, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> passes_transversality_checks(gf3)
true
```
"""
@attr Bool function passes_transversality_checks(fgs::FamilyOfG4Fluxes; check::Bool = true)
  # Entry checks
  m = model(fgs)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Transversality checks supported only for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Transversality checks supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Transversality checks supported only for toric ambient space"
  @req all(==(0), offset(fgs)) "Currently, the transversality check is only supported for flux families with trivial offset"
  # TODO: Remove this limitation, i.e. support this functionality for all flux families!

  # Extract ambient space model of g4-fluxes, in terms of which we express the generators of the flux family
  mb = chosen_g4_flux_basis(model(fgs), check = check)
  nmb = length(mb)

  # Verify that each generator of the flux family is vertical
  my_mat = matrix_integral(fgs)
  for k in 1:ncols(my_mat)
    gen_k = sum(my_mat[l,k] * mb[l] for l in 1:nmb)
    if !passes_transversality_checks(gen_k)
      return false
    end
  end
  my_mat = matrix_rational(fgs)
  for k in 1:ncols(my_mat)
    gen_k = sum(my_mat[l,k] * mb[l] for l in 1:nmb)
    if !passes_transversality_checks(gen_k)
      return false
    end
  end
  return true
end


@doc raw"""
    breaks_non_abelian_gauge_group(fgs::FamilyOfG4Fluxes; check::Bool = true)

Check if a family of G4-fluxes breaks the non-abelian
gauge group. If so, this method returns `true` and
otherwise `false`.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gf = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: broken

julia> breaks_non_abelian_gauge_group(gf)
true

julia> gf3 = special_flux_family(qsm_model, not_breaking = true, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: not broken

julia> breaks_non_abelian_gauge_group(gf3)
false

julia> m1 = matrix_integral(gf3);

julia> m2 = matrix_rational(gf3);

julia> gf4 = family_of_g4_fluxes(qsm_model, m1, m2, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> breaks_non_abelian_gauge_group(gf4, check = false)
false
```
"""
@attr Bool function breaks_non_abelian_gauge_group(fgs::FamilyOfG4Fluxes; check::Bool = true)
  # Entry checks
  m = model(fgs)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Gauge group breaking check only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Gauge group breaking check currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Gauge group breaking check currently supported only for toric ambient space"
  @req all(==(0), offset(fgs)) "Currently, the check for breaking the non-abelian gauge group is only supported for flux families with trivial offset"
  # TODO: Remove this limitation, i.e. support this functionality for all flux families!
  
  # Extract ambient space model of g4-fluxes, in terms of which we express the generators of the flux family
  mb = chosen_g4_flux_basis(model(fgs), check = check)
  nmb = length(mb)

  # Verify that each generator of the flux family does not break the non-abelian gauge group
  my_mat = matrix_integral(fgs)
  for k in 1:ncols(my_mat)
    gen_k = sum(my_mat[l,k] * mb[l] for l in 1:nmb)
    if breaks_non_abelian_gauge_group(gen_k)
      return true
    end
  end
  my_mat = matrix_rational(fgs)
  for k in 1:ncols(my_mat)
    gen_k = sum(my_mat[l,k] * mb[l] for l in 1:nmb)
    if breaks_non_abelian_gauge_group(gen_k)
      return true
    end
  end
  return false
end
