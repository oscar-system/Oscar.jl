##########################################
### (1) Fluxes within family of G4-fluxes
##########################################


@doc raw"""
    flux_instance(fgs::FamilyOfG4Fluxes, int_combination::ZZMatrix, rat_combination::QQMatrix; check::Bool = true)

Create an element of a family of G4-fluxes.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> fgs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> int_combination = matrix(ZZ, [[3]])
[3]

julia> rat_combination = matrix(QQ, [[5//2]])
[5//2]

julia> flux_instance(fgs, int_combination, rat_combination, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Tadpole cancellation check: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
```
"""
function flux_instance(fgs::FamilyOfG4Fluxes, int_combination::ZZMatrix, rat_combination::QQMatrix; check::Bool = true)
  @req ncols(int_combination) == 1 "int_combination is expected to be a single column vector"
  @req ncols(rat_combination) == 1 "rat_combination is expected to be a single column vector"
  @req nrows(int_combination) == ncols(matrix_integral(fgs)) "Number of specified integers must match the number of integral combinations in G4-flux family"
  @req nrows(rat_combination) == ncols(matrix_rational(fgs)) "Number of specified rationals must match the number of integral combinations in G4-flux family"
  m1 = matrix_integral(fgs) * int_combination
  m2 = matrix_rational(fgs) * rat_combination
  gens = ambient_space_models_of_g4_fluxes(model(fgs), check = check)
  c1 = [m1[k,1] * gens[k] for k in 1:length(gens)]
  c2 = [m2[k,1] * gens[k] for k in 1:length(gens)]
  flux = g4_flux(model(fgs), sum(c1+c2), check = check)
  if has_attribute(fgs, :is_well_quantized)
    set_attribute!(flux, :passes_elementary_quantization_checks, is_well_quantized(fgs))
  end
  if has_attribute(fgs, :is_vertical)
    set_attribute!(flux, :passes_verticality_checks, is_well_quantized(fgs))
  end
  if has_attribute(fgs, :breaks_non_abelian_gauge_group)
    set_attribute!(flux, :breaks_non_abelian_gauge_group, breaks_non_abelian_gauge_group(fgs))
  end
  return flux
end


@doc raw"""
    random_flux_instance(fgs::FamilyOfG4Fluxes; check::Bool = true)

Create a random element of a family of G4-fluxes.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> fgs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> random_flux_instance(fgs, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Tadpole cancellation check: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
```
"""
function random_flux_instance(fgs::FamilyOfG4Fluxes; check::Bool = true)
  int_combination = zero_matrix(ZZ, ncols(matrix_integral(fgs)), 1)
  rat_combination = zero_matrix(QQ, ncols(matrix_rational(fgs)), 1)
  for i in 1:ncols(matrix_integral(fgs))
    int_combination[i] = rand(-100:100)
  end
  for i in 1:ncols(matrix_rational(fgs))
    numerator = rand(-100:100)
    denominator = rand(1:100)
    rat_combination[i] = numerator // denominator
  end
  return flux_instance(fgs, int_combination, rat_combination, check = check)
end



##########################################
### (2) Special random flux on model
##########################################

@doc raw"""
    random_flux(m::AbstractFTheoryModel; vert::Bool = false, not_breaking::Bool = false, check::Bool = true)

Create a random $G_4$-flux on a given F-theory model.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> rf = random_flux(qsm_model, vert = true, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Tadpole cancellation check: not executed
  - Verticality checks: satisfied
  - Non-abelian gauge group: broken
```
"""
function random_flux(m::AbstractFTheoryModel; vert::Bool = false, not_breaking::Bool = false, check::Bool = true)
  family = special_flux_family(m, vert = vert, not_breaking = not_breaking, check = check)
  return random_flux_instance(family, check = check)
end
