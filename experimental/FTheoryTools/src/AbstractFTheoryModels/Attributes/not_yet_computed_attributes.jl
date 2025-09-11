@define_model_attribute_getter((hodge_h11, Int),
"""
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> hodge_h11(qsm_model)
31
```
""", "See [Advanced Mathematical Attributes](@ref non_yet_algorithmic_advanced_attributes) for more details.", h11)


@define_model_attribute_getter((hodge_h12, Int),
"""
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> hodge_h12(qsm_model)
10
```
""", "See [Advanced Mathematical Attributes](@ref non_yet_algorithmic_advanced_attributes) for more details.", h12)


@define_model_attribute_getter((hodge_h13, Int),
"""
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> hodge_h13(qsm_model)
34
```
""", "See [Advanced Mathematical Attributes](@ref non_yet_algorithmic_advanced_attributes) for more details.", h13)


@define_model_attribute_getter((hodge_h22, Int),
"""
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> hodge_h22(qsm_model)
284
```
""", "See [Advanced Mathematical Attributes](@ref non_yet_algorithmic_advanced_attributes) for more details.", h22)

@define_model_attribute_getter((kbar3, Int),
"""
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> kbar3(qsm_model)
6
```
""", "See [Topological Data of a `QSM`](@ref base_top_data) for more details.", Kbar3)
