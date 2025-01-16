#####################################################
# 1 Basic attributes
#####################################################

@doc raw"""
    model(gf::FamilyOfG4Fluxes)

Return the F-theory model for which this family of $G_4$-flux candidates is defined.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> f_gs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> model(f_gs) == qsm_model
true
```
"""
model(gf::FamilyOfG4Fluxes) = gf.model


@doc raw"""
    matrix_integral(gf::FamilyOfG4Fluxes)

Return the matrix whose columns specify those combinations of ambient space G4-flux
candidates, of which integral linear combinations are contained in this family
of G4-fluxes.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> f_gs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> matrix_integral(f_gs) == mat_int
true
```
"""
matrix_integral(gf::FamilyOfG4Fluxes) = gf.mat_int


@doc raw"""
    matrix_rational(gf::FamilyOfG4Fluxes)

Return the matrix whose columns specify those combinations of ambient space G4-flux
candidates, of which rational linear combinations are contained in this family
of G4-fluxes.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> f_gs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> matrix_rational(f_gs) == mat_rat
true
```
"""
matrix_rational(gf::FamilyOfG4Fluxes) = gf.mat_rat
