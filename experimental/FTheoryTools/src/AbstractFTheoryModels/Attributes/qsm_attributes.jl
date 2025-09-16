######################################################################
# (1) Attributes regarding the polytope in the Kreuzer-Skarke database
######################################################################

@define_model_attribute_getter((vertices, Vector{Vector{QQFieldElem}}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> vertices(qsm_model)
  4-element Vector{Vector{QQFieldElem}}:
   [-1, -1, -1]
   [1, -1//2, -1//2]
   [-1, 2, -1]
   [-1, -1, 5]
  ```
  """, "See [Underlying Polytope](@ref qsm_polytope) for more details.")

@define_model_attribute_getter((polytope_index, Int),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> polytope_index(qsm_model)
  4
  ```
  """, "See [Underlying Polytope](@ref qsm_polytope) for more details.", poly_index)

@define_model_attribute_getter((has_quick_triangulation, Bool),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> has_quick_triangulation(qsm_model)
  true
  ```
  """, "See [Underlying Polytope](@ref qsm_polytope) for more details.", triang_quick)

@define_model_attribute_getter((max_lattice_pts_in_facet, Int),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> max_lattice_pts_in_facet(qsm_model)
  16
  ```
  """, "See [Underlying Polytope](@ref qsm_polytope) for more details.")

@define_model_attribute_getter((estimated_number_of_triangulations, Int),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> estimated_number_of_triangulations(qsm_model)
  212533333333
  ```
  """, "See [Underlying Polytope](@ref qsm_polytope) for more details.")

######################################################################
# (2) Attributes regarding the Ci-curves
######################################################################

@define_model_attribute_getter((genera_of_ci_curves, Dict{MPolyDecRingElem,Int64}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> typeof(genera_of_ci_curves(qsm_model))
  Dict{MPolyDecRingElem, Int64}
  ```
  """, "See [The Nodal Curve](@ref qsm_nodal_curve) for more details.", genus_ci)

@define_model_attribute_getter(
  (degrees_of_kbar_restrictions_to_ci_curves, Dict{MPolyDecRingElem,Int64}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> typeof(degrees_of_kbar_restrictions_to_ci_curves(qsm_model))
  Dict{MPolyDecRingElem, Int64}
  ```
  """, "See [The Nodal Curve](@ref qsm_nodal_curve) for more details.",
  degree_of_Kbar_of_tv_restricted_to_ci)

@define_model_attribute_getter(
  (topological_intersection_numbers_among_ci_curves, Matrix{Int64}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> size(topological_intersection_numbers_among_ci_curves(qsm_model))
  (29, 29)
  ```
  """, "See [The Nodal Curve](@ref qsm_nodal_curve) for more details.",
  intersection_number_among_ci_cj)

@define_model_attribute_getter((indices_of_trivial_ci_curves, Vector{Int64}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> indices_of_trivial_ci_curves(qsm_model)
  10-element Vector{Int64}:
   23
   22
   18
   19
   20
   26
   10
   11
   12
   15
  ```
  """, "See [The Nodal Curve](@ref qsm_nodal_curve) for more details.",
  index_facet_interior_divisors)

@define_model_attribute_getter(
  (topological_intersection_numbers_among_nontrivial_ci_curves, Matrix{Int64}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> size(topological_intersection_numbers_among_nontrivial_ci_curves(qsm_model))
  (19, 19)
  ```
  """, "See [The Nodal Curve](@ref qsm_nodal_curve) for more details.",
  intersection_number_among_nontrivial_ci_cj)

######################################################################
# (3) Attributes regarding the dual graph
######################################################################

@define_model_attribute_getter((dual_graph, Graph{Undirected}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> dual_graph(qsm_model)
  Undirected graph with 21 nodes and the following edges:
  (5, 1)(6, 5)(7, 6)(8, 7)(9, 4)(9, 8)(10, 1)(11, 4)(12, 3)(12, 10)(13, 3)(13, 11)(14, 1)(15, 4)(16, 3)(17, 3)(18, 2)(18, 14)(19, 2)(19, 15)(20, 2)(20, 16)(21, 2)(21, 17)
  ```
  """, "See [The Dual Graph](@ref qsm_dual_graph) for more details.")

@define_model_attribute_getter((components_of_dual_graph, Vector{String}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> length(components_of_dual_graph(qsm_model))
  21
  ```
  """, "See [The Dual Graph](@ref qsm_dual_graph) for more details.")

@define_model_attribute_getter(
  (degrees_of_kbar_restrictions_to_components_of_dual_graph, Dict{String,Int64}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> typeof(degrees_of_kbar_restrictions_to_components_of_dual_graph(qsm_model))
  Dict{String, Int64}
  ```
  """, "See [The Dual Graph](@ref qsm_dual_graph) for more details.",
  degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph)

@define_model_attribute_getter((genera_of_components_of_dual_graph, Dict{String,Int64}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> typeof(genera_of_components_of_dual_graph(qsm_model))
  Dict{String, Int64}
  ```
  """, "See [The Dual Graph](@ref qsm_dual_graph) for more details.",
  genus_of_components_of_dual_graph)

######################################################################
# (4) Attributes regarding the simplified dual graph
######################################################################

@define_model_attribute_getter((simplified_dual_graph, Graph{Undirected}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> simplified_dual_graph(qsm_model)
  Undirected graph with 4 nodes and the following edges:
  (2, 1)(3, 1)(3, 2)(4, 1)(4, 2)(4, 3)
  ```
  """, "See [The Simplified Dual Graph](@ref qsm_simple_dual_graph) for more details.")

@define_model_attribute_getter((components_of_simplified_dual_graph, Vector{String}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> length(components_of_simplified_dual_graph(qsm_model))
  4
  ```
  """, "See [The Simplified Dual Graph](@ref qsm_simple_dual_graph) for more details.")

@define_model_attribute_getter(
  (degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph, Dict{String,Int64}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> typeof(degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph(qsm_model))
  Dict{String, Int64}
  ```
  """, "See [The Simplified Dual Graph](@ref qsm_simple_dual_graph) for more details.",
  degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph)

@define_model_attribute_getter(
  (genera_of_components_of_simplified_dual_graph, Dict{String,Int64}),
  """
  ```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
  julia> using Random;

  julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
  Hypersurface model over a concrete base

  julia> typeof(genera_of_components_of_simplified_dual_graph(qsm_model))
  Dict{String, Int64}
  ```
  """, "See [The Simplified Dual Graph](@ref qsm_simple_dual_graph) for more details.",
  genus_of_components_of_simplified_dual_graph)
