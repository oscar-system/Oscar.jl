#####################################################
# 1 Basic attributes
#####################################################

@doc raw"""
    model(gf::G4Flux)

Return the F-theory model for which this $G_4$-flux candidate is defined.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> cohomology_ring(ambient_space(qsm_model), check = false);

julia> g4_class = cohomology_class(anticanonical_divisor_class(ambient_space(qsm_model)))^2;

julia> g4f = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Tadpole cancellation check: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed  

julia> model(g4f)
Hypersurface model over a concrete base
```
"""
model(gf::G4Flux) = gf.model


@doc raw"""
    cohomology_class(gf::G4Flux)

Return the cohomology class which defines the $G_4$-flux candidate.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> cohomology_ring(ambient_space(qsm_model), check = false);

julia> g4_class = cohomology_class(anticanonical_divisor_class(ambient_space(qsm_model)))^2;

julia> g4f = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Tadpole cancellation check: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> cohomology_class(g4f) == g4_class
true
```
"""
cohomology_class(gf::G4Flux) = gf.class
